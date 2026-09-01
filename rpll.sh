#!/bin/bash

rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
          inp1l_*.txt synth_base_*.txt tmp_round*.txt PIS*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt bres.tmp \
          out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt bextxzx.txt suminpxzx.txt data_middle1.txt idk.txt \
          expfy_samples_done.tmp .expfy_samples_*.tmp
rm -rf baseline .generate_outinp_* .run_round.*

set -e

# Node-dependent work is calibrated to a 1,000-node graph.  EFFECTIVE_N is
# clamped so small graphs still get enough trials for stable estimates while
# large graphs cannot make every candidate evaluation arbitrarily expensive.
REFERENCE_N=1000
MIN_EFFECTIVE_N=500
MAX_EFFECTIVE_N=5000
MAX_APPLY_SAMPLES=100000

# The final optimization round always uses these tightest constraints. hdgb
# is expressed at the 1,000-node reference size and is scaled after the target
# graph has been cleaned and EFFECTIVE_N is known. TIGHTEST_DEGB is the integer
# code passed to expfy; expfy divides it by 140, so 14 means a log-degree
# tolerance of 0.1.
TIGHTEST_EVB=0.5
TIGHTEST_DEGB=14
TIGHTEST_HDGB_REFERENCE=40
ROUND_BACKWARD_FACTOR=1.25

# After the requested RSEAS rounds, run one four-round equal-weight RMSE
# polishing block.  The command-line round count still means RSEAS rounds.
EXTRA_RMSE_ROUNDS=4

clamp_int() {
    local VALUE="$1"
    local LOWER="$2"
    local UPPER="$3"

    if (( VALUE < LOWER )); then VALUE="$LOWER"; fi
    if (( VALUE > UPPER )); then VALUE="$UPPER"; fi
    printf '%d\n' "$VALUE"
}

# cleanup renumbers the largest connected component as 0,1,...,n-1, which is
# also exactly how expfy determines n (maximum endpoint plus one).
graph_node_count() {
    awk '
        NF >= 2 {
            if ($1 !~ /^[0-9]+$/ || $2 !~ /^[0-9]+$/) bad = 1
            u = $1 + 0
            v = $2 + 0
            if (!seen || u > largest) largest = u
            if (!seen || v > largest) largest = v
            seen = 1
        }
        END {
            if (bad || !seen) exit 1
            print largest + 1
        }
    ' "$1"
}

scale_from_reference() {
    local BASE="$1"
    # Add half the denominator so integer arithmetic rounds instead of floors.
    printf '%d\n' $(( (BASE * EFFECTIVE_N + REFERENCE_N / 2) / REFERENCE_N ))
}

float_below_by() {
    awk -v A="$1" -v B="$2" -v D="$3" 'BEGIN { exit !(A < B - D) }'
}

float_above_by() {
    awk -v A="$1" -v B="$2" -v D="$3" 'BEGIN { exit !(A > B + D) }'
}


# Graphlet-frequency logarithms use the same floor in Bash/Python and in the
# four phase-specific C++ optimizers.
LOG_FLOOR=-15.0

# Round-level RSEAS state.  The six positions correspond to graphlets
# 5, 6, 7, 8, 9, and 10, in the order emitted by bnb 4.
RSEAS_E=0.0
RSEAS_E_INITIALIZED=0
RSEAS_WEIGHTS=(1 1 4 30 1 2)
PREVIOUS_WORST_GRAPHLET=""

# Log displacement: log(A)-log(B), using LOG_FLOOR for a zero frequency.
displ() {
    python3 - "$1" "$2" "$3" "$LOG_FLOOR" << 'EOF_PY'
import math
import sys

path_a, path_b, output_path, log_floor_text = sys.argv[1:5]
log_floor = float(log_floor_text)

A = [float(x) for x in open(path_a, encoding="utf-8").read().split()]
B = [float(x) for x in open(path_b, encoding="utf-8").read().split()]
if len(A) < 6 or len(B) < 6:
    raise SystemExit("Expected at least six graphlet-frequency values")

def safe_log(value):
    return math.log(value) if value > 0.0 else log_floor

diff = [safe_log(A[i]) - safe_log(B[i]) for i in range(6)]
open(output_path, "w", encoding="utf-8").write(
    " ".join(format(value, ".17g") for value in diff) + "\n"
)
EOF_PY
}

# Score two six-coordinate displacement vectors with one fixed RSEAS state.
# The caller supplies e followed by six positive weights.
rseas_with_config() {
    local FILE_A="$1"
    local FILE_B="$2"
    local ERROR_RANGE="$3"
    shift 3

    if (( $# != 6 )); then
        echo "rseas_with_config requires exactly six weights." >&2
        return 1
    fi

    python3 - "$FILE_A" "$FILE_B" "$ERROR_RANGE" "$@" << 'EOF_PY'
import math
import sys

path_a, path_b = sys.argv[1:3]
error_range = float(sys.argv[3])
weights = [float(x) for x in sys.argv[4:10]]

A = [float(x) for x in open(path_a, encoding="utf-8").read().split()]
B = [float(x) for x in open(path_b, encoding="utf-8").read().split()]
if len(A) < 6 or len(B) < 6:
    raise SystemExit("Expected at least six values in each RSEAS vector")
if not math.isfinite(error_range) or error_range < 0.0:
    raise SystemExit("RSEAS error range must be a finite nonnegative number")
if len(weights) != 6 or any((not math.isfinite(w) or w <= 0.0) for w in weights):
    raise SystemExit("RSEAS requires six positive finite weights")

distances = [abs(A[i] - B[i]) for i in range(6)]
maximum_distance = max(distances)

if any(distance > error_range for distance in distances):
    result = sum(
        weights[i] * max(0.0, distances[i] - error_range) ** 2
        for i in range(6)
    )
else:
    result = maximum_distance - error_range

print(format(result, ".17g"))
EOF_PY
}

rseas() {
    rseas_with_config "$1" "$2" "$RSEAS_E" "${RSEAS_WEIGHTS[@]}"
}

# Equal-weight RMSE used in the old pipeline.  With six unit weights this is
# sqrt(sum_i (A_i-B_i)^2).  The missing division by six is retained deliberately
# so reported values remain comparable with earlier runs; it does not affect
# candidate ordering or the optimizer's chosen scales.
unit_rmse() {
    python3 - "$1" "$2" << 'EOF_PY'
import math
import sys

path_a, path_b = sys.argv[1:3]
A = [float(x) for x in open(path_a, encoding="utf-8").read().split()]
B = [float(x) for x in open(path_b, encoding="utf-8").read().split()]
if len(A) < 6 or len(B) < 6:
    raise SystemExit("Expected at least six values in each RMSE vector")

value = sum((A[i] - B[i]) ** 2 for i in range(6))
print(format(math.sqrt(value), ".17g"))
EOF_PY
}

# run_round sets ROUND_OBJECTIVE to RSEAS or RMSE.  Bash's dynamic scoping and
# process inheritance make that fixed phase choice visible inside its workers.
objective_score() {
    case "$ROUND_OBJECTIVE" in
        RSEAS) rseas "$1" "$2" ;;
        RMSE)  unit_rmse "$1" "$2" ;;
        *)
            echo "Unknown round objective '$ROUND_OBJECTIVE'." >&2
            return 1
            ;;
    esac
}

# Print: worst_graphlet max_log_error candidate_error_range. The candidate
# range is 0.75*max_log_error; set_round_rseas_state then caps it at the
# preceding round's e. Ties are broken by taking the first coordinate, so
# graphlet 5 wins a complete tie, then 6, and so on.
log_error_summary() {
    python3 - "$1" "$2" "$LOG_FLOOR" << 'EOF_PY'
import math
import sys

current_path, target_path, log_floor_text = sys.argv[1:4]
log_floor = float(log_floor_text)

current = [float(x) for x in open(current_path, encoding="utf-8").read().split()]
target = [float(x) for x in open(target_path, encoding="utf-8").read().split()]
if len(current) < 6 or len(target) < 6:
    raise SystemExit("Expected six graphlet frequencies when deriving RSEAS state")

def safe_log(value):
    return math.log(value) if value > 0.0 else log_floor

distances = [
    abs(safe_log(current[i]) - safe_log(target[i]))
    for i in range(6)
]
worst_index = max(range(6), key=lambda index: distances[index])
maximum_distance = distances[worst_index]
error_range = 0.75 * maximum_distance

print(
    worst_index + 5,
    format(maximum_distance, ".17g"),
    format(error_range, ".17g"),
)
EOF_PY
}

set_round_rseas_state() {
    local ROUND_INDEX="$1"
    local CURRENT_VEC="$2"
    local TARGET_VEC="$3"
    local WORST_GRAPHLET MAX_ERROR NEW_ERROR_RANGE WEIGHT_INDEX
    local PREVIOUS_ERROR_RANGE=""
    local REPEATED_TEXT=""
    local ERROR_RANGE_TEXT=""

    read -r WORST_GRAPHLET MAX_ERROR NEW_ERROR_RANGE \
        < <(log_error_summary "$CURRENT_VEC" "$TARGET_VEC")

    WEIGHT_INDEX=$((WORST_GRAPHLET - 5))
    if [[ -n "$PREVIOUS_WORST_GRAPHLET" && \
          "$WORST_GRAPHLET" == "$PREVIOUS_WORST_GRAPHLET" ]]; then
        RSEAS_WEIGHTS[$WEIGHT_INDEX]=$((RSEAS_WEIGHTS[$WEIGHT_INDEX] + 1))
        REPEATED_TEXT="; repeated worst graphlet, increased its weight"
    fi

    PREVIOUS_WORST_GRAPHLET="$WORST_GRAPHLET"

    # The first RSEAS round has no preceding e. After that, e is monotone
    # nonincreasing:
    #     e = min(0.75 * maxLogDist, previous_e).
    if (( RSEAS_E_INITIALIZED == 0 )); then
        RSEAS_E="$NEW_ERROR_RANGE"
        RSEAS_E_INITIALIZED=1
    else
        PREVIOUS_ERROR_RANGE="$RSEAS_E"
        RSEAS_E=$(awk \
            -v CANDIDATE="$NEW_ERROR_RANGE" \
            -v PREVIOUS="$PREVIOUS_ERROR_RANGE" \
            'BEGIN {
                printf "%.17g\n", (CANDIDATE < PREVIOUS ? CANDIDATE : PREVIOUS)
            }')

        if awk -v CANDIDATE="$NEW_ERROR_RANGE" -v PREVIOUS="$PREVIOUS_ERROR_RANGE" \
            'BEGIN { exit !(CANDIDATE > PREVIOUS) }'; then
            ERROR_RANGE_TEXT="; candidate e=$NEW_ERROR_RANGE capped at previous e=$PREVIOUS_ERROR_RANGE"
        fi
    fi

    echo "Round $ROUND_INDEX/$ROUNDS RSEAS state: worst graphlet=$WORST_GRAPHLET, max log error=$MAX_ERROR, e=$RSEAS_E, weights=${RSEAS_WEIGHTS[*]}$REPEATED_TEXT$ERROR_RANGE_TEXT"
}

configure_node_scaled_parameters() {
    if ! NODE_COUNT=$(graph_node_count "$TARGET"); then
        echo "Could not determine the node count from $TARGET." >&2
        exit 1
    fi
    if (( NODE_COUNT < 4 )); then
        echo "The cleaned target must contain at least 4 nodes." >&2
        exit 1
    fi

    EFFECTIVE_N=$(clamp_int "$NODE_COUNT" "$MIN_EFFECTIVE_N" "$MAX_EFFECTIVE_N")

    # Initial synthesis passes.  These are direct expfy success targets, so
    # they scale with n.  The first, unusually heavy pass gets an extra cap.
    INIT_PRIMARY_SAMPLES=$(clamp_int "$(scale_from_reference 50000)" 10000 100000)
    INIT_MAIN_SAMPLES=$(clamp_int "$(scale_from_reference 10000)" 2500 50000)
    INIT_SECONDARY_SAMPLES=$(clamp_int "$(scale_from_reference 1000)" 250 5000)

    # Preliminary eigenvector rewiring and the per-round EVC probes.
    PIS1_SAMPLES=$(scale_from_reference 3000)
    PIS2_SAMPLES=$(scale_from_reference 1000)
    PIS3_SAMPLES=$(scale_from_reference 400)
    OUTPEV_SMALL_SAMPLES=$(scale_from_reference 100)
    OUTPEV_LARGE_SAMPLES=$(scale_from_reference 200)

    # hdgb is an absolute count tolerance in a degree-histogram bin, so it
    # should scale with node count. evb is dimensionless. degb is an encoded
    # log-degree tolerance; expfy divides it by 140 and computes its separate
    # n-dependent absolute fallback internally, so rpll does not rescale it.
    INIT_HDGB_1=$(scale_from_reference 50)
    INIT_HDGB_2=$(scale_from_reference 40)
    INIT_HDGB_3=$(scale_from_reference 30)
    INIT_HDGB_4=$(scale_from_reference 26)
    INIT_HDGB_5=$(scale_from_reference 24)
    INIT_HDGB_6=$(scale_from_reference 23)
    INIT_HDGB_7=$(scale_from_reference 22)
    # The last optimization round uses this value exactly.  Earlier rounds
    # are generated backward from it by multiplying by ROUND_BACKWARD_FACTOR.
    TIGHTEST_HDGB=$(scale_from_reference "$TIGHTEST_HDGB_REFERENCE")

    echo "Graph size: $NODE_COUNT nodes (workload scaling uses $EFFECTIVE_N; clamp $MIN_EFFECTIVE_N..$MAX_EFFECTIVE_N)."
    echo "PIS samples: $PIS1_SAMPLES, $PIS2_SAMPLES, $PIS3_SAMPLES; outpev samples: $OUTPEV_SMALL_SAMPLES, $OUTPEV_LARGE_SAMPLES."
}

# Build the constraint schedule backward from the final round. For evb the
# previous value is exactly 1.25 times the next value. The encoded degb value
# and hdgb are positive integers, so each backward multiplication is rounded
# to the nearest integer. Consequently, ROUNDS=1 uses the tightest values
# immediately, while
# every larger run reaches those same values on its final iteration.
declare -a ROUND_EVB_SCHEDULE
declare -a ROUND_DEGB_SCHEDULE
declare -a ROUND_HDGB_SCHEDULE

build_round_constraint_schedule() {
    local ROUND_INDEX NEXT_INDEX NEXT_EVB NEXT_DEGB NEXT_HDGB
    local PREVIOUS_EVB PREVIOUS_DEGB PREVIOUS_HDGB

    ROUND_EVB_SCHEDULE=()
    ROUND_DEGB_SCHEDULE=()
    ROUND_HDGB_SCHEDULE=()

    ROUND_EVB_SCHEDULE[$ROUNDS]="$TIGHTEST_EVB"
    ROUND_DEGB_SCHEDULE[$ROUNDS]="$TIGHTEST_DEGB"
    ROUND_HDGB_SCHEDULE[$ROUNDS]="$TIGHTEST_HDGB"

    for ((ROUND_INDEX=ROUNDS-1; ROUND_INDEX>=1; ROUND_INDEX--)); do
        NEXT_INDEX=$((ROUND_INDEX + 1))
        NEXT_EVB="${ROUND_EVB_SCHEDULE[$NEXT_INDEX]}"
        NEXT_DEGB="${ROUND_DEGB_SCHEDULE[$NEXT_INDEX]}"
        NEXT_HDGB="${ROUND_HDGB_SCHEDULE[$NEXT_INDEX]}"

        read -r PREVIOUS_EVB PREVIOUS_DEGB PREVIOUS_HDGB \
            < <(awk \
                -v EVB="$NEXT_EVB" \
                -v DEGB="$NEXT_DEGB" \
                -v HDGB="$NEXT_HDGB" \
                -v FACTOR="$ROUND_BACKWARD_FACTOR" \
                'BEGIN {
                    printf "%.12g %d %d\n", \
                           EVB * FACTOR, \
                           int(DEGB * FACTOR + 0.5), \
                           int(HDGB * FACTOR + 0.5)
                }')

        ROUND_EVB_SCHEDULE[$ROUND_INDEX]="$PREVIOUS_EVB"
        ROUND_DEGB_SCHEDULE[$ROUND_INDEX]="$PREVIOUS_DEGB"
        ROUND_HDGB_SCHEDULE[$ROUND_INDEX]="$PREVIOUS_HDGB"
    done
}

set_round_constraints() {
    local ROUND_INDEX="$1"

    evb="${ROUND_EVB_SCHEDULE[$ROUND_INDEX]}"
    degb="${ROUND_DEGB_SCHEDULE[$ROUND_INDEX]}"
    hdgb="${ROUND_HDGB_SCHEDULE[$ROUND_INDEX]}"

    echo "RSEAS round $ROUND_INDEX/$ROUNDS constraints: evb=$evb, degb=$degb, hdgb=$hdgb"
}

set_rmse_round_constraints() {
    local ROUND_INDEX="$1"

    # The RSEAS schedule already reaches these tightest values on its last
    # round.  Keep them fixed during all four appended polishing rounds.
    evb="$TIGHTEST_EVB"
    degb="$TIGHTEST_DEGB"
    hdgb="$TIGHTEST_HDGB"

    echo "RMSE polishing round $ROUND_INDEX/$EXTRA_RMSE_ROUNDS constraints: evb=$evb, degb=$degb, hdgb=$hdgb"
}

INITIAL_GRAPH_GIVEN=0

if [ $# -eq 3 ]; then
    TARGET0="$1"
    FINAL="$2"
    ROUNDS="$3"
    PRE_INIT_SYNTH=""
elif [ $# -eq 4 ]; then
    TARGET0="$1"
    FINAL="$2"
    PRE_INIT_SYNTH="$3"
    ROUNDS="$4"
    INITIAL_GRAPH_GIVEN=1
else
    echo "Usage: $0 <Target> <Output> <Rounds> | $0 <Target> <Output> <Base_Synth> <Rounds>"
    exit 1
fi

if ! [[ "$ROUNDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Rounds must be a positive integer, not '$ROUNDS'." >&2
    exit 1
fi

TARGET="target_c.txt"
./cleanup "$TARGET0" "$TARGET"
configure_node_scaled_parameters
build_round_constraint_schedule

if [ $# -eq 3 ]; then
    ./gen_deg2 "$TARGET" synth_base_1.txt

    ./bnb 4 synth_base_1.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_1.txt "$INIT_PRIMARY_SAMPLES" "$INIT_PRIMARY_SAMPLES" 1 8 2.0 20 "$INIT_HDGB_1" > synth_base_2.txt
    ./bnb 4 synth_base_2.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_2.txt "$INIT_MAIN_SAMPLES" "$INIT_MAIN_SAMPLES" 1 8 1.7 18 "$INIT_HDGB_2" > synth_base_3.txt
    ./bnb 4 synth_base_3.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_3.txt "$INIT_MAIN_SAMPLES" "$INIT_MAIN_SAMPLES" 1 8 1.5 18 "$INIT_HDGB_3" > synth_base_4.txt
    ./bnb 4 synth_base_4.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_4.txt "$INIT_MAIN_SAMPLES" "$INIT_MAIN_SAMPLES" 1 8 1.3 17 "$INIT_HDGB_4" > synth_base_5.txt
    ./bnb 4 synth_base_5.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_5.txt "$INIT_MAIN_SAMPLES" "$INIT_SECONDARY_SAMPLES" 1 17 1.2 16 "$INIT_HDGB_5" > synth_base_6.txt
    ./bnb 4 synth_base_6.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_6.txt "$INIT_MAIN_SAMPLES" "$INIT_SECONDARY_SAMPLES" 1 17 1.1 15 "$INIT_HDGB_6" > synth_base_7.txt
    ./bnb 4 synth_base_7.txt > bres.tmp
    cat data_middle.txt > data_middle1.txt
    ./expfy "$TARGET" synth_base_7.txt "$INIT_MAIN_SAMPLES" "$INIT_SECONDARY_SAMPLES" 1 17 1 14 "$INIT_HDGB_7" > synth_base_8.txt

    PRE_INIT_SYNTH=synth_base_8.txt
else
    if ! SYNTH_NODE_COUNT=$(graph_node_count "$PRE_INIT_SYNTH"); then
        echo "Could not determine the node count from $PRE_INIT_SYNTH." >&2
        exit 1
    fi
    if (( SYNTH_NODE_COUNT != NODE_COUNT )); then
        echo "Target has $NODE_COUNT nodes but base synth has $SYNTH_NODE_COUNT; expfy requires equal node counts." >&2
        exit 1
    fi
fi

# When a base synth was supplied by the user, each complete PIS candidate is
# accepted only when its RSEAS is no greater than that of the last accepted
# graph. PIS is outside the optimization rounds, so it uses six unit weights
# and one error range fixed from the graph entering PIS.
PIS_TARGET_VEC="pis_target_vec.tmp"
PIS_CURRENT_RSEAS=""
PIS_RSEAS_E=0.0
PIS_RSEAS_WEIGHTS=(1 1 1 1 1 1)

pis_graph_rseas() {
    local GRAPH="$1"
    local VEC_FILE DISPL_FILE SCORE STATUS

    VEC_FILE=$(mktemp ".pis_vec.XXXXXX")
    DISPL_FILE=$(mktemp ".pis_displ.XXXXXX")
    if ! ./bnb 4 "$GRAPH" > "$VEC_FILE"; then
        rm -f "$VEC_FILE" "$DISPL_FILE"
        return 1
    fi

    displ "$VEC_FILE" "$PIS_TARGET_VEC" "$DISPL_FILE"
    SCORE=$(rseas_with_config f0.txt "$DISPL_FILE" "$PIS_RSEAS_E" \
        "${PIS_RSEAS_WEIGHTS[@]}") || {
        STATUS=$?
        rm -f "$VEC_FILE" "$DISPL_FILE"
        return "$STATUS"
    }

    rm -f "$VEC_FILE" "$DISPL_FILE"
    printf '%s\n' "$SCORE"
}

accept_pis_candidate() {
    local BEFORE="$1"
    local CANDIDATE="$2"
    local OUTPUT="$3"
    local LABEL="$4"

    if (( ! INITIAL_GRAPH_GIVEN )); then
        mv "$CANDIDATE" "$OUTPUT"
        return
    fi

    # The threshold branch can deliberately produce an exact copy. Do not
    # re-sample the same graph and accidentally change the remembered score.
    if cmp -s "$BEFORE" "$CANDIDATE"; then
        mv "$CANDIDATE" "$OUTPUT"
        echo "$LABEL made no change; RSEAS remains $PIS_CURRENT_RSEAS."
        return
    fi

    local CANDIDATE_RSEAS
    CANDIDATE_RSEAS=$(pis_graph_rseas "$CANDIDATE")

    echo "$LABEL RSEAS: before=$PIS_CURRENT_RSEAS, candidate=$CANDIDATE_RSEAS"

    if awk -v NEW="$CANDIDATE_RSEAS" -v OLD="$PIS_CURRENT_RSEAS" \
        'BEGIN { exit !(NEW <= OLD) }'; then
        mv "$CANDIDATE" "$OUTPUT"
        PIS_CURRENT_RSEAS="$CANDIDATE_RSEAS"
        echo "$LABEL accepted."
    else
        cp "$BEFORE" "$OUTPUT"
        rm -f "$CANDIDATE"
        echo "$LABEL rejected; keeping the previous graph."
    fi
}

if (( INITIAL_GRAPH_GIVEN )); then
    local_pis_initial_vec=$(mktemp ".pis_initial_vec.XXXXXX")
    local_pis_initial_displ=$(mktemp ".pis_initial_displ.XXXXXX")

    ./bnb 4 "$TARGET" > "$PIS_TARGET_VEC"
    ./bnb 4 "$PRE_INIT_SYNTH" > "$local_pis_initial_vec"

    read -r PIS_WORST_GRAPHLET PIS_MAX_ERROR PIS_RSEAS_E \
        < <(log_error_summary "$local_pis_initial_vec" "$PIS_TARGET_VEC")

    displ "$local_pis_initial_vec" "$PIS_TARGET_VEC" "$local_pis_initial_displ"
    PIS_CURRENT_RSEAS=$(rseas_with_config f0.txt "$local_pis_initial_displ" \
        "$PIS_RSEAS_E" "${PIS_RSEAS_WEIGHTS[@]}")

    rm -f "$local_pis_initial_vec" "$local_pis_initial_displ"

    echo "PIS RSEAS state: worst graphlet=$PIS_WORST_GRAPHLET, max log error=$PIS_MAX_ERROR, e=$PIS_RSEAS_E, weights=${PIS_RSEAS_WEIGHTS[*]}"
    echo "Initial graph RSEAS before PIS: $PIS_CURRENT_RSEAS"
fi

PIS1_CANDIDATE="PIS1_candidate.txt"
./evs "$PRE_INIT_SYNTH" > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
if float_below_by "$EV1" "$EV2" 14; then
    ./evc1 "$PRE_INIT_SYNTH" "$PIS1_SAMPLES" > "$PIS1_CANDIDATE"
elif float_above_by "$EV1" "$EV2" 14; then
    ./evc2 "$PRE_INIT_SYNTH" "$PIS1_SAMPLES" > "$PIS1_CANDIDATE"
else
    cp "$PRE_INIT_SYNTH" "$PIS1_CANDIDATE"
fi
accept_pis_candidate "$PRE_INIT_SYNTH" "$PIS1_CANDIDATE" PIS1.txt PIS1

./evs PIS1.txt > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
PIS2_CANDIDATE="PIS2_candidate.txt"
if float_below_by "$EV1" "$EV2" 7; then
    ./evc1 PIS1.txt "$PIS2_SAMPLES" > "$PIS2_CANDIDATE"
elif float_above_by "$EV1" "$EV2" 7; then
    ./evc2 PIS1.txt "$PIS2_SAMPLES" > "$PIS2_CANDIDATE"
else
    cp PIS1.txt "$PIS2_CANDIDATE"
fi
accept_pis_candidate PIS1.txt "$PIS2_CANDIDATE" PIS2.txt PIS2

./evs PIS2.txt > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
PIS3_CANDIDATE="PIS3_candidate.txt"
if float_below_by "$EV1" "$EV2" 2; then
    ./evc1 PIS2.txt "$PIS3_SAMPLES" > "$PIS3_CANDIDATE"
elif float_above_by "$EV1" "$EV2" 2; then
    ./evc2 PIS2.txt "$PIS3_SAMPLES" > "$PIS3_CANDIDATE"
else
    cp PIS2.txt "$PIS3_CANDIDATE"
fi
accept_pis_candidate PIS2.txt "$PIS3_CANDIDATE" PIS3.txt PIS3

rm -f "$PIS_TARGET_VEC"
INIT_SYNTH=PIS3.txt

# Initialize the globals used by apply_exp/apply_exp2.  The main loop resets
# them from the precomputed schedule at the start of every round.
evb="${ROUND_EVB_SCHEDULE[1]}"
degb="${ROUND_DEGB_SCHEDULE[1]}"
hdgb="${ROUND_HDGB_SCHEDULE[1]}"
CANDS=($(seq 1 60 | grep -v -E '^(25|26)$'))
ACTIVE_CANDS=("${CANDS[@]}")

# Map each transformation number to the graphlet ID required by expfy.
# Transformations 25 and 26 are unused.  Transformation 60 starts from a
# five-node graphlet, so it cannot be classified from bnb 4's columns.
declare -a initial_graphlet
for i in {1..9} 29 34 35 37 38 52; do initial_graphlet[$i]=6; done
for i in {10..17} 28 33 41 42 51 56; do initial_graphlet[$i]=7; done
for i in 18 19 20 30 36 39 53 54; do initial_graphlet[$i]=5; done
for i in 21 22 23 27 32 40 50 57 58; do initial_graphlet[$i]=8; done
for i in 24 31 48 49 55; do initial_graphlet[$i]=9; done
for i in {43..47} 59; do initial_graphlet[$i]=10; done
initial_graphlet[25]=-1
initial_graphlet[26]=-1
initial_graphlet[60]=-1

# data_middle1.txt column index -> graphlet ID.  Several columns can describe
# different node roles in the same graphlet.
column_graphlet=(0 0 0 0 6 6 5 5 8 7 7 7 9 9 10)

declare -a graphlet_node_count
declare -a rare_graphlet

# Number of simultaneous workers used by generate_outinp.
# Override when launching, for example: PARALLEL_JOBS=4 ./rpll_parallel_round.sh ...
PARALLEL_JOBS="${PARALLEL_JOBS:-2}"

# Number of simultaneous workers used inside run_round.  By default it uses
# the same limit as generate_outinp, but it can be tuned independently.
ROUND_JOBS="${ROUND_JOBS:-$PARALLEL_JOBS}"

declare -a ec
for i in {1..60}; do
    ec[$i]=0
done
ec[27]=2;
ec[28]=2;
ec[29]=3;
ec[30]=3;
ec[31]=1;
ec[32]=1;
ec[33]=1;
ec[34]=2;
ec[35]=2;
ec[36]=2;
ec[37]=1;
ec[38]=1;
ec[39]=1;
ec[40]=-1;
ec[41]=-1;
ec[42]=-1;
ec[43]=-1;
ec[44]=-2;
ec[45]=-2;
ec[46]=-3;
ec[47]=-3;
ec[48]=-1;
ec[49]=-1;
ec[50]=-1;
ec[51]=-1;
ec[52]=-1;
ec[53]=-1;
ec[54]=0;
ec[55]=0;
ec[56]=0;
ec[57]=0;
ec[58]=-2;
ec[59]=+3;
ec[60]=+1;

BASELINE_DIR="baseline"
mkdir -p "$BASELINE_DIR"

# Create an isolated current directory for a parallel worker.  Several of the
# project executables create fixed-name scratch files such as data_middle.txt,
# so merely changing the shell output filename is not enough to prevent races.
prepare_parallel_worker() {
    local WORK_DIR="$1"
    local ROOT_DIR="$2"
    local SHARED_FILE TOOL

    mkdir -p "$WORK_DIR"

    for SHARED_FILE in data_middle.txt data_middle1.txt f0.txt; do
        if [ -f "$ROOT_DIR/$SHARED_FILE" ]; then
            cp "$ROOT_DIR/$SHARED_FILE" "$WORK_DIR/$SHARED_FILE"
        fi
    done

    # bnb and expfy may invoke other executables through paths such as ./blant.
    # Linking every project executable keeps those calls working inside the
    # private directory without sharing their generated scratch files.
    for TOOL in "$ROOT_DIR"/*; do
        if [ ! -f "$TOOL" ]; then
            continue
        fi

        if [ -x "$TOOL" ]; then
            ln -s "$TOOL" "$WORK_DIR/${TOOL##*/}"
        else
            case "$TOOL" in
                *.py|*.awk|*.pl)
                    ln -s "$TOOL" "$WORK_DIR/${TOOL##*/}"
                    ;;
            esac
        fi
    done
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

float_less() {
    awk -v A="$1" -v B="$2" 'BEGIN { exit !(A < B) }'
}

float_greater() {
    awk -v A="$1" -v B="$2" 'BEGIN { exit !(A > B) }'
}

integer_part() {
    awk -v X="$1" 'BEGIN { print int(X) }'
}

# Count how many nodes were seen in each graphlet represented by the columns of
# data_middle1.txt.  If several columns map to the same graphlet, a node is
# counted only once when any of those columns is positive.  A graphlet is rare
# when it appears at at most 0.1% of the nodes, i.e. count * 1000 <= n.
refresh_rare_graphlets() {
    local DATA_FILE="$1"
    local SUMMARY
    local COLUMN_MAP="${column_graphlet[*]}"

    if [ ! -r "$DATA_FILE" ]; then
        echo "Cannot read graphlet-node data from $DATA_FILE." >&2
        return 1
    fi

    if ! SUMMARY=$(awk -v column_map="$COLUMN_MAP" '
        BEGIN {
            expected_columns = split(column_map, graphlet_for_column, " ")
            for (i = 1; i <= expected_columns; i++) {
                graphlet_exists[graphlet_for_column[i]] = 1
            }
        }
        NF == 0 { next }
        NF != expected_columns {
            printf "Inconsistent column count in %s at line %d: expected %d, got %d.\n", \
                   FILENAME, NR, expected_columns, NF > "/dev/stderr"
            bad = 1
            next
        }
        {
            n++
            for (graphlet in graphlet_exists) present[graphlet] = 0

            for (i = 1; i <= expected_columns; i++) {
                if (($i + 0) > 0) present[graphlet_for_column[i]] = 1
            }

            for (graphlet in graphlet_exists) {
                if (present[graphlet]) seen[graphlet]++
            }
        }
        END {
            if (bad || n == 0) exit 1
            printf "%d", n
            for (graphlet = 0; graphlet <= 10; graphlet++) {
                if (graphlet in graphlet_exists) {
                    printf " %d:%d", graphlet, seen[graphlet] + 0
                }
            }
            printf "\n"
        }
    ' "$DATA_FILE"); then
        return 1
    fi

    local -a STATS=()
    read -r -a STATS <<< "$SUMMARY"

    local NODE_COUNT="${STATS[0]}"
    local ENTRY GRAPHLET COUNT VAL INITIAL INDEX
    local -a RARE_IDS=()
    local -a SKIPPED=()

    graphlet_node_count=()
    rare_graphlet=()

    echo "Graphlet node coverage (rare means at most 0.1% of $NODE_COUNT nodes):"
    for ((INDEX=1; INDEX<${#STATS[@]}; INDEX++)); do
        ENTRY="${STATS[$INDEX]}"
        GRAPHLET="${ENTRY%%:*}"
        COUNT="${ENTRY#*:}"
        graphlet_node_count[$GRAPHLET]="$COUNT"

        if (( COUNT * 1000 <= NODE_COUNT )); then
            rare_graphlet[$GRAPHLET]=1
            RARE_IDS+=("$GRAPHLET")
            echo "  graphlet $GRAPHLET: $COUNT/$NODE_COUNT nodes (rare)"
        else
            rare_graphlet[$GRAPHLET]=0
            echo "  graphlet $GRAPHLET: $COUNT/$NODE_COUNT nodes"
        fi
    done

    ACTIVE_CANDS=()
    for VAL in "${CANDS[@]}"; do
        if [[ ! -v "initial_graphlet[$VAL]" ]]; then
            echo "No initial-graphlet mapping for transformation $VAL." >&2
            return 1
        fi

        INITIAL="${initial_graphlet[$VAL]}"
        if (( INITIAL < 0 )); then
            # Transformation 60 uses a five-node initial graphlet, which is not
            # represented in data_middle1.txt produced by bnb 4.
            ACTIVE_CANDS+=("$VAL")
        elif [[ ! -v "graphlet_node_count[$INITIAL]" ]]; then
            echo "Transformation $VAL maps to graphlet $INITIAL, which has no data_middle1.txt columns." >&2
            return 1
        elif (( rare_graphlet[INITIAL] )); then
            SKIPPED+=("$VAL")
            rm -f "outinp_val_${VAL}.txt"
        else
            ACTIVE_CANDS+=("$VAL")
        fi
    done

    if (( ${#RARE_IDS[@]} == 0 )); then
        echo "Rare graphlets: none"
    else
        echo "Rare graphlets: ${RARE_IDS[*]}"
    fi

    if (( ${#SKIPPED[@]} == 0 )); then
        echo "Transformations skipped for rarity: none"
    else
        echo "Transformations skipped for rarity: ${SKIPPED[*]}"
    fi
    echo "Active transformations: ${ACTIVE_CANDS[*]}"

    if (( ${#ACTIVE_CANDS[@]} == 0 )); then
        echo "No transformations remain after rare-graphlet filtering." >&2
        return 1
    fi
}

# Convert the magnitude used by bxk4f (calibrated at 1,000 nodes) to an
# expfy success target.  For ec==0 this is
#
#     requested * EFFECTIVE_N / (1000 * 600),
#
# exactly the n/1000 * 1/600 rule.  Transformations with nonzero ec retain the
# old ten-to-one reduction because each success also changes the edge count.
scale_apply_samples() {
    local REQUESTED="$1"
    local DIVISOR="$2"
    local DENOMINATOR
    local SCALED

    if (( REQUESTED <= 0 )); then
        printf '0\n'
        return
    fi

    DENOMINATOR=$(( REFERENCE_N * DIVISOR ))
    SCALED=$(( (REQUESTED * EFFECTIVE_N + DENOMINATOR / 2) / DENOMINATOR ))

    # A positive optimizer request must perform at least one modification, but
    # an extreme magnitude must not turn a single candidate into an unbounded
    # run.  expfy also has its own attempted-sample limit.
    clamp_int "$SCALED" 1 "$MAX_APPLY_SAMPLES"
}

# Convert a completed expfy success count back into bxk4f scale units.
# If optimizer scale X was converted to REQUESTED_SAMPLES=f(X), but expfy
# completed only DONE_SAMPLES, the effective tested scale is
#
#     X * DONE_SAMPLES / f(X).
#
# Integer arithmetic rounds to the nearest scale unit.
completed_scale_from_samples() {
    local REQUESTED_SCALE="$1"
    local REQUESTED_SAMPLES="$2"
    local DONE_SAMPLES="$3"

    if (( REQUESTED_SCALE <= 0 || REQUESTED_SAMPLES <= 0 || DONE_SAMPLES <= 0 )); then
        printf '0\n'
        return
    fi

    if (( DONE_SAMPLES > REQUESTED_SAMPLES )); then
        DONE_SAMPLES="$REQUESTED_SAMPLES"
    fi

    printf '%d\n' $((
        (REQUESTED_SCALE * DONE_SAMPLES + REQUESTED_SAMPLES / 2) /
        REQUESTED_SAMPLES
    ))
}

read_completed_sample_file() {
    local SAMPLE_FILE="$1"
    local RESULT_VAR1="$2"
    local RESULT_VAR2="$3"
    local VALUE1 VALUE2 EXTRA

    if ! read -r VALUE1 VALUE2 EXTRA < "$SAMPLE_FILE"; then
        echo "Could not read completed sample counts from $SAMPLE_FILE." >&2
		cat "$SAMPLE_FILE"
        return 1
    fi

    if ! [[ "$VALUE1" =~ ^[0-9]+$ && "$VALUE2" =~ ^[0-9]+$ ]] || [ -n "$EXTRA" ]; then
        echo "Invalid completed sample data in $SAMPLE_FILE: expected two nonnegative integers." >&2
        return 1
    fi

    printf -v "$RESULT_VAR1" '%s' "$VALUE1"
    printf -v "$RESULT_VAR2" '%s' "$VALUE2"
}

# Apply expfy with one type.  The fifth argument is the name of a caller
# variable that receives the effective scale actually completed by expfy.
apply_exp() {
    local INPUT_FILE="$1"
    local REQUESTED_SCALE="$2"
    local VAL="$3"
    local OUTPUT_FILE="$4"
    local TESTED_SCALE_VAR="$5"
    local DIVISOR REQUESTED_SAMPLES SAMPLE_FILE DONE1 DONE2 TESTED_SCALE

    if [ "${ec[$VAL]}" -eq 0 ]; then
        DIVISOR=60
    else
        DIVISOR=600
    fi
    REQUESTED_SAMPLES=$(scale_apply_samples "$REQUESTED_SCALE" "$DIVISOR")

    SAMPLE_FILE=$(mktemp ".expfy_samples.XXXXXX.tmp")
    if ! ./expfy "$TARGET" "$INPUT_FILE" "$REQUESTED_SAMPLES" 0 \
        "$VAL" 1 "$evb" "$degb" "$hdgb" "$SAMPLE_FILE" \
        > "$OUTPUT_FILE"; then
        rm -f "$SAMPLE_FILE" "$OUTPUT_FILE"
        return 1
    fi

    if ! read_completed_sample_file "$SAMPLE_FILE" DONE1 DONE2; then
        rm -f "$SAMPLE_FILE" "$OUTPUT_FILE"
        return 1
    fi
    rm -f "$SAMPLE_FILE"

    TESTED_SCALE=$(completed_scale_from_samples \
        "$REQUESTED_SCALE" "$REQUESTED_SAMPLES" "$DONE1")
    printf -v "$TESTED_SCALE_VAR" '%s' "$TESTED_SCALE"
}

# Apply expfy with two types.  The final two arguments are caller-variable
# names that receive the two effective scales actually completed by expfy.
apply_exp2() {
    local INPUT_FILE="$1"
    local REQUESTED_SCALE1="$2"
    local REQUESTED_SCALE2="$3"
    local VAL="$4"
    local VAL2="$5"
    local OUTPUT_FILE="$6"
    local TESTED_SCALE_VAR1="$7"
    local TESTED_SCALE_VAR2="$8"
    local DIVISOR1 DIVISOR2 REQUESTED_SAMPLES1 REQUESTED_SAMPLES2
    local SAMPLE_FILE DONE1 DONE2 TESTED_SCALE1 TESTED_SCALE2

    if [ "${ec[$VAL]}" -eq 0 ]; then
        DIVISOR1=60
    else
        DIVISOR1=600
    fi
    if [ "${ec[$VAL2]}" -eq 0 ]; then
        DIVISOR2=60
    else
        DIVISOR2=600
    fi

    REQUESTED_SAMPLES1=$(scale_apply_samples "$REQUESTED_SCALE1" "$DIVISOR1")
    REQUESTED_SAMPLES2=$(scale_apply_samples "$REQUESTED_SCALE2" "$DIVISOR2")

    SAMPLE_FILE=$(mktemp ".expfy_samples.XXXXXX.tmp")
    if ! ./expfy "$TARGET" "$INPUT_FILE" \
        "$REQUESTED_SAMPLES1" "$REQUESTED_SAMPLES2" "$VAL" "$VAL2" \
        "$evb" "$degb" "$hdgb" "$SAMPLE_FILE" > "$OUTPUT_FILE"; then
        rm -f "$SAMPLE_FILE" "$OUTPUT_FILE"
        return 1
    fi

    if ! read_completed_sample_file "$SAMPLE_FILE" DONE1 DONE2; then
        rm -f "$SAMPLE_FILE" "$OUTPUT_FILE"
        return 1
    fi
    rm -f "$SAMPLE_FILE"

    TESTED_SCALE1=$(completed_scale_from_samples \
        "$REQUESTED_SCALE1" "$REQUESTED_SAMPLES1" "$DONE1")
    TESTED_SCALE2=$(completed_scale_from_samples \
        "$REQUESTED_SCALE2" "$REQUESTED_SAMPLES2" "$DONE2")

    printf -v "$TESTED_SCALE_VAR1" '%s' "$TESTED_SCALE1"
    printf -v "$TESTED_SCALE_VAR2" '%s' "$TESTED_SCALE2"
}


# optimization round
run_round() {
    # Keep the cleanup trap and worker-local variables inside a subshell.  The
    # next optimization round still starts only after this one has produced its
    # output, because rounds depend on one another and cannot safely overlap.
    (
        set -e

        local INPUT_SYNTH="$1"
        local OUTPUT_SYNTH="$2"
        local ROUND_OBJECTIVE="$3"

        case "$ROUND_OBJECTIVE" in
            RSEAS|RMSE) ;;
            *)
                echo "run_round objective must be RSEAS or RMSE, not '$ROUND_OBJECTIVE'." >&2
                exit 1
                ;;
        esac
        local ROOT_DIR
        ROOT_DIR=$(pwd -P)

        local MAX_JOBS="$ROUND_JOBS"
        if ! is_positive_integer "$MAX_JOBS"; then
            echo "ROUND_JOBS must be a positive integer, not '$MAX_JOBS'." >&2
            exit 1
        fi

        local INPUT_ABS OUTPUT_ABS TARGET_ABS VEC_TAR_ABS F0_ABS
        local BNB_ABS BXK4F_ABS BXK4ONE_ABS
        INPUT_ABS=$(realpath "$INPUT_SYNTH")
        OUTPUT_ABS=$(realpath -m "$OUTPUT_SYNTH")
        TARGET_ABS=$(realpath "$TARGET")
        VEC_TAR_ABS=$(realpath vecTar.txt)
        F0_ABS=$(realpath f0.txt)
        BNB_ABS=$(realpath ./bnb)
        if [ "$ROUND_OBJECTIVE" = "RSEAS" ]; then
            BXK4F_ABS=$(realpath ./bxk4f_rseas)
            BXK4ONE_ABS=$(realpath ./bxk4one_rseas)
        else
            BXK4F_ABS=$(realpath ./bxk4f_rmse)
            BXK4ONE_ABS=$(realpath ./bxk4one_rmse)
        fi

        local ROUND_DIR
        ROUND_DIR=$(mktemp -d "$ROOT_DIR/.run_round.XXXXXX")
        trap 'rm -rf -- "$ROUND_DIR"' EXIT

        local VEC_INS="$ROOT_DIR/vecInS.txt"
        local INPL="$ROOT_DIR/inpl.txt"
        local INP="$ROOT_DIR/inp.txt"

        "$BNB_ABS" 4 "$INPUT_ABS" > "$VEC_INS"
        cp "$ROOT_DIR/data_middle.txt" "$ROOT_DIR/data_middle1.txt"

        displ "$VEC_INS" "$VEC_TAR_ABS" "$INPL"
        cat "$VEC_INS" "$VEC_TAR_ABS" > "$INP"

        local had_score
        had_score=$(objective_score "$F0_ABS" "$INPL")

        echo "New $ROUND_OBJECTIVE round with input $INPUT_SYNTH"
        echo "Running independent round work with up to $MAX_JOBS parallel jobs"

        # ---------------------------------------------------------
        # 1. Generate and evaluate the four EVC candidates in parallel.
        # ---------------------------------------------------------
        local EVC_DIR="$ROUND_DIR/evc"
        mkdir -p "$EVC_DIR"

        local -a EVC_LABELS=(outpev11 outpev12 outpev21 outpev22)
        local -a EVC_TOOLS=(evc1 evc1 evc2 evc2)
        local -a EVC_AMOUNTS=(
            "$OUTPEV_SMALL_SAMPLES" "$OUTPEV_LARGE_SAMPLES"
            "$OUTPEV_SMALL_SAMPLES" "$OUTPEV_LARGE_SAMPLES"
        )
        local -a EVC_PIDS=()
        local EVC_FAILED=0
        local IDX LABEL TOOL AMOUNT PID

        for IDX in "${!EVC_LABELS[@]}"; do
            LABEL="${EVC_LABELS[$IDX]}"
            TOOL="${EVC_TOOLS[$IDX]}"
            AMOUNT="${EVC_AMOUNTS[$IDX]}"

            (
                set -e
                local WORK_DIR="$EVC_DIR/work_$LABEL"
                prepare_parallel_worker "$WORK_DIR" "$ROOT_DIR"
                trap 'rm -rf -- "$WORK_DIR"' EXIT
                cd "$WORK_DIR"

                "./$TOOL" "$INPUT_ABS" "$AMOUNT" > candidate.graph
                ./bnb 4 candidate.graph > candidate.vec
                displ "$VEC_INS" candidate.vec candidate.displ

                local SCORE
                SCORE=$(objective_score candidate.displ "$INPL")

                mv candidate.graph "$EVC_DIR/$LABEL.graph"
                printf '%s|%s\n' "$SCORE" "$EVC_DIR/$LABEL.graph" \
                    > "$EVC_DIR/$LABEL.result.tmp"
                mv "$EVC_DIR/$LABEL.result.tmp" "$EVC_DIR/$LABEL.result"
            ) &
            EVC_PIDS+=("$!")

            if (( ${#EVC_PIDS[@]} >= MAX_JOBS )); then
                if ! wait "${EVC_PIDS[0]}"; then
                    EVC_FAILED=1
                fi
                EVC_PIDS=("${EVC_PIDS[@]:1}")
            fi
        done

        for PID in "${EVC_PIDS[@]}"; do
            if ! wait "$PID"; then
                EVC_FAILED=1
            fi
        done

        if (( EVC_FAILED )); then
            echo "At least one parallel EVC evaluation failed." >&2
            exit 1
        fi

        local BEST_EVC_GRAPH="$INPUT_ABS"
        local BEST_EVC_SCORE="$had_score"
        local CUR_SCORE CUR_GRAPH

        # Read in the original deterministic order so equal scores are handled
        # exactly as they were in the serial version.
        for LABEL in "${EVC_LABELS[@]}"; do
            IFS='|' read -r CUR_SCORE CUR_GRAPH < "$EVC_DIR/$LABEL.result"
            echo "$LABEL $ROUND_OBJECTIVE = $CUR_SCORE"

            if float_less "$CUR_SCORE" "$BEST_EVC_SCORE"; then
                BEST_EVC_SCORE="$CUR_SCORE"
                BEST_EVC_GRAPH="$CUR_GRAPH"
            fi
        done

        if [ "$BEST_EVC_GRAPH" != "$INPUT_ABS" ]; then
            echo "Replacing $INPUT_SYNTH with the best EVC candidate ($ROUND_OBJECTIVE $BEST_EVC_SCORE)"
            cp "$BEST_EVC_GRAPH" "$INPUT_ABS"
            had_score="$BEST_EVC_SCORE"
        else
            echo "No EVC candidate improved the input graph."
        fi

        # Recompute the current graph vector after the possible EVC replacement.
        "$BNB_ABS" 4 "$INPUT_ABS" > "$VEC_INS"
        cp "$ROOT_DIR/data_middle.txt" "$ROOT_DIR/data_middle1.txt"
        displ "$VEC_INS" "$VEC_TAR_ABS" "$INPL"
        cat "$VEC_INS" "$VEC_TAR_ABS" > "$INP"

        # ---------------------------------------------------------
        # 2. Score all candidate pairs using a fixed worker pool.
        # ---------------------------------------------------------
        local PAIR_DIR="$ROUND_DIR/pairs"
        local RESULTS_FILE="$ROUND_DIR/round_results.tmp"
        mkdir -p "$PAIR_DIR"
        : > "$RESULTS_FILE"

        local CAND_COUNT=${#ACTIVE_CANDS[@]}
        local PAIR_WORKERS="$MAX_JOBS"
        if (( PAIR_WORKERS > CAND_COUNT )); then
            PAIR_WORKERS=$CAND_COUNT
        fi

        local -a PAIR_PIDS=()
        local PAIR_FAILED=0
        local WORKER

        for ((WORKER=0; WORKER<PAIR_WORKERS; WORKER++)); do
            (
                set -e
                local WORK_DIR="$PAIR_DIR/work_$WORKER"
                prepare_parallel_worker "$WORK_DIR" "$ROOT_DIR"
                trap 'rm -rf -- "$WORK_DIR"' EXIT
                cd "$WORK_DIR"

                local LOCAL_RESULTS="pair_results.tmp"
                : > "$LOCAL_RESULTS"

                local VI VAL VAL2 MAG MAG2 EXPECTED
                local OUT_INP OUT_INP2

                # Round-robin assignment keeps the triangular pair workload
                # reasonably balanced among the workers.
                for ((VI=WORKER; VI<CAND_COUNT; VI+=PAIR_WORKERS)); do
                    VAL="${ACTIVE_CANDS[$VI]}"

                    for VAL2 in "${ACTIVE_CANDS[@]}"; do
                        if (( VAL2 > VAL )); then
                            continue
                        fi

                        OUT_INP="$ROOT_DIR/outinp_val_${VAL}.txt"
                        OUT_INP2="$ROOT_DIR/outinp_val_${VAL2}.txt"

                        if [ ! -s "$OUT_INP" ] || [ ! -s "$OUT_INP2" ]; then
                            echo "Missing transformation data for pair $VAL,$VAL2." >&2
                            exit 1
                        fi

                        # bxk4f input order is:
                        #   current[6], target[6],
                        #   first.ka first.kb first.A[6] first.B[6],
                        #   second.ka second.kb second.A[6] second.B[6].
                        # Each outinp file begins with its actual tested ka/kb.
                        cat "$INP" "$OUT_INP" "$OUT_INP2" > pair_input.tmp

                        local PAIR_INPUT_WORDS
                        PAIR_INPUT_WORDS=$(wc -w < pair_input.tmp)
                        if (( PAIR_INPUT_WORDS != 40 )); then
                            echo "Malformed bxk4f input for pair $VAL,$VAL2: expected 40 values, got $PAIR_INPUT_WORDS." >&2
                            exit 1
                        fi

                        if [ "$ROUND_OBJECTIVE" = "RSEAS" ]; then
                            "$BXK4F_ABS" "${ec[$VAL]}" "${ec[$VAL2]}" "$RSEAS_E" \
                                "${RSEAS_WEIGHTS[@]}" \
                                < pair_input.tmp > pair_output.tmp
                        else
                            "$BXK4F_ABS" "${ec[$VAL]}" "${ec[$VAL2]}" \
                                < pair_input.tmp > pair_output.tmp
                        fi

                        read -r MAG MAG2 EXPECTED < pair_output.tmp
                        printf '%s|%s|%s|%s|%s\n' \
                            "$EXPECTED" "$VAL" "$VAL2" "$MAG" "$MAG2" \
                            >> "$LOCAL_RESULTS"
                    done
                done

                mv "$LOCAL_RESULTS" "$PAIR_DIR/results_$WORKER.txt"
            ) &
            PAIR_PIDS+=("$!")
        done

        for PID in "${PAIR_PIDS[@]}"; do
            if ! wait "$PID"; then
                PAIR_FAILED=1
            fi
        done

        if (( PAIR_FAILED )); then
            echo "At least one parallel pair-scoring worker failed." >&2
            exit 1
        fi

        for ((WORKER=0; WORKER<PAIR_WORKERS; WORKER++)); do
            cat "$PAIR_DIR/results_$WORKER.txt" >> "$RESULTS_FILE"
        done

        had_score=$(objective_score "$F0_ABS" "$INPL")
        echo "$ROUND_OBJECTIVE before modification is $had_score"

        # ---------------------------------------------------------
        # 3. Evaluate the ten most promising pairs in parallel.
        # ---------------------------------------------------------
        local TOP_DIR="$ROUND_DIR/top"
        mkdir -p "$TOP_DIR"

        local -a TOP_LINES=()
        mapfile -t TOP_LINES < <(sort -n "$RESULTS_FILE" | head -n 10)

        if (( ${#TOP_LINES[@]} == 0 )); then
            echo "No pair-scoring results were produced." >&2
            exit 1
        fi

        local -a TOP_PIDS=()
        local TOP_FAILED=0
        local RANK LINE
        local sorted_score sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2
        local SCALE SCALE2

        for IDX in "${!TOP_LINES[@]}"; do
            RANK=$((IDX + 1))
            LINE="${TOP_LINES[$IDX]}"
            IFS='|' read -r sorted_score sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2 <<< "$LINE"

            SCALE=$(integer_part "$sorted_MAG")
            SCALE2=$(integer_part "$sorted_MAG2")

            if (( SCALE <= 0 && SCALE2 <= 0 && RANK > 1 )); then
                echo "Skipping rank $RANK pair $sorted_VAL,$sorted_VAL2 because both scales are non-positive."
                continue
            fi

            if (( SCALE <= 0 && SCALE2 <= 0 && RANK == 1 )); then
                echo "First candidate has non-positive scales. Forcing evaluation."
            fi

            (
                set -e
                local WORK_DIR="$TOP_DIR/work_$RANK"
                prepare_parallel_worker "$WORK_DIR" "$ROOT_DIR"
                trap 'rm -rf -- "$WORK_DIR"' EXIT
                cd "$WORK_DIR"

                # apply_exp2 refers to TARGET dynamically and invokes ./expfy.
                local TARGET="$TARGET_ABS"
                local STAGE1="candidate_stage1.graph"
                local FINAL_CANDIDATE="candidate_final.graph"
                local MAG_FACTOR GOT_SCORE
                local ADJUSTED_SCALE="$SCALE"
                local ADJUSTED_SCALE2="$SCALE2"
                local STAGE1_TESTED_SCALE STAGE1_TESTED_SCALE2
                local FINAL_TESTED_SCALE FINAL_TESTED_SCALE2

                apply_exp2 "$INPUT_ABS" "$ADJUSTED_SCALE" "$ADJUSTED_SCALE2" \
                    "$sorted_VAL" "$sorted_VAL2" "$STAGE1" \
                    STAGE1_TESTED_SCALE STAGE1_TESTED_SCALE2

                ./bnb 4 "$STAGE1" > stage1.vec
                displ "$VEC_INS" stage1.vec stage1.displ
                cat "$INP" stage1.displ > scale_input.tmp
                if [ "$ROUND_OBJECTIVE" = "RSEAS" ]; then
                    "$BXK4ONE_ABS" "$RSEAS_E" "${RSEAS_WEIGHTS[@]}" \
                        < scale_input.tmp > scale_output.tmp
                else
                    "$BXK4ONE_ABS" < scale_input.tmp > scale_output.tmp
                fi

                read -r MAG_FACTOR < scale_output.tmp

                # bxk4one fits a multiplier to the displacement actually seen in
                # STAGE1.  Multiply the effective completed scales, not the
                # larger scales originally requested from expfy.
                ADJUSTED_SCALE=$(awk -v S="$STAGE1_TESTED_SCALE" -v M="$MAG_FACTOR" \
                    'BEGIN { print int(S * M) }')
                ADJUSTED_SCALE2=$(awk -v S="$STAGE1_TESTED_SCALE2" -v M="$MAG_FACTOR" \
                    'BEGIN { print int(S * M) }')

                apply_exp2 "$INPUT_ABS" "$ADJUSTED_SCALE" "$ADJUSTED_SCALE2" \
                    "$sorted_VAL" "$sorted_VAL2" "$FINAL_CANDIDATE" \
                    FINAL_TESTED_SCALE FINAL_TESTED_SCALE2

                ./bnb 4 "$FINAL_CANDIDATE" > final.vec
                displ "$VEC_INS" final.vec final.displ
                GOT_SCORE=$(objective_score final.displ "$INPL")

                # Preserve the exact graph whose objective score was measured.  expfy is
                # randomized, so regenerating later would be slower and could
                # produce a different graph from the candidate selected here.
                local SAVED_CANDIDATE
                SAVED_CANDIDATE="$TOP_DIR/candidate_$(printf '%02d' "$RANK").graph"
                mv "$FINAL_CANDIDATE" "$SAVED_CANDIDATE"

                # Store both the adjusted parameters and the saved graph path.
                # The result file is published only after the graph is safely
                # in TOP_DIR, so the parent never selects a missing candidate.
                printf '%s|%s|%s|%s|%s|%s|%s|%s\n' \
                    "$GOT_SCORE" "$sorted_VAL" "$sorted_VAL2" \
                    "$FINAL_TESTED_SCALE" "$FINAL_TESTED_SCALE2" "$sorted_score" "$RANK" \
                    "$SAVED_CANDIDATE" \
                    > "$TOP_DIR/result_$(printf '%02d' "$RANK").tmp"
                mv "$TOP_DIR/result_$(printf '%02d' "$RANK").tmp" \
                   "$TOP_DIR/result_$(printf '%02d' "$RANK").txt"
            ) &
            TOP_PIDS+=("$!")

            if (( ${#TOP_PIDS[@]} >= MAX_JOBS )); then
                if ! wait "${TOP_PIDS[0]}"; then
                    TOP_FAILED=1
                fi
                TOP_PIDS=("${TOP_PIDS[@]:1}")
            fi
        done

        for PID in "${TOP_PIDS[@]}"; do
            if ! wait "$PID"; then
                TOP_FAILED=1
            fi
        done

        if (( TOP_FAILED )); then
            echo "At least one parallel top-pair evaluation failed." >&2
            exit 1
        fi

        local BEST_VAL=""
        local BEST_VAL2=""
        local BEST_MAG=""
        local BEST_MAG2=""
        local BEST_SCORE=""
        local BEST_GRAPH=""
        local GOT_SCORE EXPECTED_SCORE RESULT_FILE CANDIDATE_GRAPH

        # Consume results by original rank, preserving the serial tie-breaking.
        for ((RANK=1; RANK<=${#TOP_LINES[@]}; RANK++)); do
            RESULT_FILE="$TOP_DIR/result_$(printf '%02d' "$RANK").txt"
            if [ ! -f "$RESULT_FILE" ]; then
                continue
            fi

            IFS='|' read -r GOT_SCORE sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2 \
                EXPECTED_SCORE _ CANDIDATE_GRAPH < "$RESULT_FILE"

            if [ ! -s "$CANDIDATE_GRAPH" ]; then
                echo "Saved graph for rank $RANK is missing or empty: $CANDIDATE_GRAPH" >&2
                exit 1
            fi

            echo "VAL=$sorted_VAL,$sorted_VAL2 MAG=$sorted_MAG,MAG2=$sorted_MAG2. expected $ROUND_OBJECTIVE=$EXPECTED_SCORE, got $ROUND_OBJECTIVE=$GOT_SCORE."

            if [ -z "$BEST_SCORE" ] || float_less "$GOT_SCORE" "$BEST_SCORE"; then
                BEST_SCORE="$GOT_SCORE"
                BEST_VAL="$sorted_VAL"
                BEST_VAL2="$sorted_VAL2"
                BEST_MAG="$sorted_MAG"
                BEST_MAG2="$sorted_MAG2"
                BEST_GRAPH="$CANDIDATE_GRAPH"
            fi
        done

        if [ -z "$BEST_SCORE" ]; then
            echo "No top-pair candidate completed successfully." >&2
            exit 1
        fi

        echo "Best pair was $BEST_VAL and $BEST_VAL2 with $ROUND_OBJECTIVE $BEST_SCORE and parameters $BEST_MAG,$BEST_MAG2"

        if float_greater "$BEST_SCORE" "$had_score"; then
            echo "The best evaluated pair was worse than the input ($ROUND_OBJECTIVE $had_score); keeping the input graph."
            cp "$INPUT_ABS" "$OUTPUT_ABS"
        else
            echo "Using the saved best candidate directly; expfy will not be rerun."
            cp "$BEST_GRAPH" "$OUTPUT_ABS"
        fi

        echo "Achieved graphlet frequency:"
        "$BNB_ABS" 4 "$OUTPUT_ABS"
    )
}

generate_outinp() {
        local InS="$1"
        local ROOT_DIR
        ROOT_DIR=$(pwd -P)

        local INS_ABS TARGET_ABS VEC_INS_ABS EXPFY_ABS BNB_ABS
        INS_ABS=$(realpath "$InS")
        TARGET_ABS=$(realpath "$TARGET")
        VEC_INS_ABS=$(realpath vecInS.txt)
        EXPFY_ABS=$(realpath ./expfy)
        BNB_ABS=$(realpath ./bnb)

        local MAX_JOBS="$PARALLEL_JOBS"
        if ! [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]; then
                echo "PARALLEL_JOBS must be a positive integer, not '$MAX_JOBS'." >&2
                return 1
        fi

        echo
        echo "Re-getting transformation results with $MAX_JOBS parallel jobs"
        echo

        local -a PIDS=()
        local FAILED=0
        local VAL

        for VAL in "${ACTIVE_CANDS[@]}"; do
                (
                        set -e

                        local WORK_DIR
                        WORK_DIR=$(mktemp -d "$ROOT_DIR/.generate_outinp_${VAL}.XXXXXX")
                        trap 'rm -rf -- "$WORK_DIR"' EXIT

                        # expfy/bnb use files in the current directory.  Each candidate gets
                        # a private directory so their temporary files cannot overwrite one another.
                        for SHARED_FILE in data_middle.txt data_middle1.txt; do
                                if [ -f "$ROOT_DIR/$SHARED_FILE" ]; then
                                        cp "$ROOT_DIR/$SHARED_FILE" "$WORK_DIR/$SHARED_FILE"
                                fi
                        done

                        # Also expose every executable from the project directory.  This covers
                        # wrappers such as bnb that may launch another local executable (for
                        # example ./blant) while keeping their generated data files job-local.
                        for TOOL in "$ROOT_DIR"/*; do
                                if [ -f "$TOOL" ] && [ -x "$TOOL" ]; then
                                        ln -s "$TOOL" "$WORK_DIR/${TOOL##*/}"
                                fi
                        done

                        # These two must exist even if their executable bits are unusual.
                        [ -e "$WORK_DIR/expfy" ] || ln -s "$EXPFY_ABS" "$WORK_DIR/expfy"
                        [ -e "$WORK_DIR/bnb" ] || ln -s "$BNB_ABS" "$WORK_DIR/bnb"

                        cd "$WORK_DIR"

                        # Bash uses dynamic scoping, so apply_exp sees this absolute TARGET.
                        local TARGET="$TARGET_ABS"
                        local BASE_SYNTH1="baseline1_val_${VAL}.txt"
                        local BASE_SYNTH2="baseline2_val_${VAL}.txt"
                        local OUT_INP="outinp_val_${VAL}.txt"
                        local SCALE TESTED_SCALE1 TESTED_SCALE2

                        echo "Getting $VAL result"

                        SCALE=30000
                        apply_exp "$INS_ABS" "$SCALE" "$VAL" "$BASE_SYNTH1" \
                            TESTED_SCALE1
                        ./bnb 4 "$BASE_SYNTH1" > "vecBS1_val_${VAL}.txt"
                        displ "$VEC_INS_ABS" "vecBS1_val_${VAL}.txt" "trfBS1_val_${VAL}.txt"

                        SCALE=300000
                        apply_exp "$INS_ABS" "$SCALE" "$VAL" "$BASE_SYNTH2" \
                            TESTED_SCALE2
                        ./bnb 4 "$BASE_SYNTH2" > "vecBS2_val_${VAL}.txt"
                        displ "$VEC_INS_ABS" "vecBS2_val_${VAL}.txt" "trfBS2_val_${VAL}.txt"

                        # Each transformation file is self-contained:
                        #   line 1: effective tested scales ka kb
                        #   line 2: displacement measured at ka
                        #   line 3: displacement measured at kb
                        printf '%s %s\n' "$TESTED_SCALE1" "$TESTED_SCALE2" > "$OUT_INP"
                        cat "trfBS1_val_${VAL}.txt" "trfBS2_val_${VAL}.txt" >> "$OUT_INP"

                        # bxk4f expects exactly 14 values per transformation:
                        #   ka kb + six values at ka + six values at kb.
                        # Catch malformed probe data here instead of letting bxk4f
                        # fail later with only "incomplete input".
                        local OUT_INP_WORDS
                        OUT_INP_WORDS=$(wc -w < "$OUT_INP")
                        if (( OUT_INP_WORDS != 14 )); then
                                echo "Malformed $OUT_INP for transformation $VAL: expected 14 values (ka kb A[6] B[6]), got $OUT_INP_WORDS." >&2
                                exit 1
                        fi

                        # mv is atomic on the same filesystem, so run_round never sees a half-written file.
                        mv "$OUT_INP" "$ROOT_DIR/outinp_val_${VAL}.txt"
                        echo "Finished $VAL result (tested scales: $TESTED_SCALE1, $TESTED_SCALE2)"
                ) &

                PIDS+=("$!")

                # Keep at most MAX_JOBS candidates alive at once.
                if (( ${#PIDS[@]} >= MAX_JOBS )); then
                        if ! wait "${PIDS[0]}"; then
                                FAILED=1
                        fi
                        PIDS=("${PIDS[@]:1}")
                fi
        done

        # Wait for the final partial batch.
        local PID
        for PID in "${PIDS[@]}"; do
                if ! wait "$PID"; then
                        FAILED=1
                fi
        done

        if (( FAILED )); then
                echo "At least one parallel generate_outinp job failed." >&2
                return 1
        fi
}

# Main loop
CURRENT="$INIT_SYNTH"
./bnb 4 "$TARGET" > vecTar.txt

# Phase 1: the requested number of dynamic-weight RSEAS rounds.
for ((i=1; i<=ROUNDS; i++)); do
    set_round_constraints "$i"
    ./bnb 4 "$CURRENT" > vecInS.txt
    set_round_rseas_state "$i" vecInS.txt vecTar.txt

    if (( (i-1) % 4 == 0 )); then
        ./bnb 4 "$CURRENT" > bres.tmp
        cat data_middle.txt > data_middle1.txt
        refresh_rare_graphlets data_middle1.txt
        generate_outinp "$CURRENT"
    fi

    NEXT="tmp_round${i}.txt"
    run_round "$CURRENT" "$NEXT" RSEAS
    CURRENT="$NEXT"
done

# Phase 2: regenerate transformation probes from the final RSEAS graph, then
# perform exactly four equal-weight RMSE polishing rounds.  These rounds do not
# change or consult the accumulated RSEAS weights/error range.
echo
echo "Completed $ROUNDS RSEAS rounds. Starting $EXTRA_RMSE_ROUNDS equal-weight RMSE polishing rounds."
evb="$TIGHTEST_EVB"
degb="$TIGHTEST_DEGB"
hdgb="$TIGHTEST_HDGB"
./bnb 4 "$CURRENT" > vecInS.txt
cat data_middle.txt > data_middle1.txt
refresh_rare_graphlets data_middle1.txt
generate_outinp "$CURRENT"

for ((j=1; j<=EXTRA_RMSE_ROUNDS; j++)); do
    set_rmse_round_constraints "$j"
    TOTAL_INDEX=$((ROUNDS + j))
    NEXT="tmp_round${TOTAL_INDEX}.txt"
    run_round "$CURRENT" "$NEXT" RMSE
    CURRENT="$NEXT"
done

cp "$CURRENT" "$FINAL"
echo
echo "Final output stored in $FINAL after $ROUNDS RSEAS rounds and $EXTRA_RMSE_ROUNDS RMSE rounds."
rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
          inp1l_*.txt synth_base_*.txt tmp_round*.txt PIS*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt bres.tmp \
          out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt bextxzx.txt suminpxzx.txt data_middle1.txt idk.txt \
          expfy_samples_done.tmp .expfy_samples_*.tmp
rm -rf baseline .generate_outinp_* .run_round.*
