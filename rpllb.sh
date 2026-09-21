#!/bin/bash

# expfyc implements modes 1..95; expfyb currently stops at 83 (25/26 unused).
# Resolve and query the actual executable before cleanup or entering workers.
EXPFY_BIN=$(realpath -- "${EXPFY_BIN:-./expfy}") || exit 1
if [ ! -x "$EXPFY_BIN" ]; then
    echo "Transformation backend is not executable: $EXPFY_BIN. Run ./build.sh first." >&2
    exit 1
fi
if ! EXPFY_MODE_LIMIT=$("$EXPFY_BIN" --max-transformation); then
    echo "Cannot query transformation support from $EXPFY_BIN. Rebuild with ./build.sh." >&2
    exit 1
fi
if ! [[ "$EXPFY_MODE_LIMIT" =~ ^[1-9][0-9]{0,8}$ ]]; then
    echo "Invalid transformation limit from $EXPFY_BIN: '$EXPFY_MODE_LIMIT'." >&2
    exit 1
fi

rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
          inp1l_*.txt synth_base_*.txt tmp_round*.txt PIS*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt bres.tmp evh.txt \
          out_inplxzx.txt data1.txt vecInS.txt vecTar.txt bextxzx.txt suminpxzx.txt data_middle1.txt idk.txt \
          expfy_samples_done.tmp .expfy_samples_*.tmp
rm -rf baseline .generate_outinp_* .run_round.*

set -e

# Set graphlet frequency checker
g++ -O2 bnb.cpp -o bnb

# Launch from the project directory containing bnb, blant, canon_maps, and
# orca_jesse_blant_table. Candidate work directories must not become BLANT's cwd.
RPLLB_ROOT_DIR=$(pwd -P)

# Dense-graph setup: compute node orbits once on the initial synthetic graph.
# Keep this data_middle1.txt unchanged for all later selection and filtering.
# Use a headered temporary copy rather than modifying the input graph.
initialize_node_weights() (
    set -e
    local GRAPH_ABS="$1"
    local ORCA_INPUT ORCA_OUTPUT ORBIT_NODE_COUNT
    ORCA_INPUT=$(mktemp .rpllb_orca_input.XXXXXX) || exit 1
    ORCA_OUTPUT=""
    trap 'rm -f -- "$ORCA_INPUT" ${ORCA_OUTPUT:+"$ORCA_OUTPUT"}' EXIT
    ORCA_OUTPUT=$(mktemp .rpllb_orca_output.XXXXXX) || exit 1

    ORBIT_NODE_COUNT=$(python3 - "$GRAPH_ABS" "$ORCA_INPUT" << 'EOF_ORCA_INPUT'
import sys

source, destination = sys.argv[1:]
edges = set()
with open(source, encoding="utf-8") as graph:
    for line in graph:
        if not line.strip():
            continue
        u, v = map(int, line.split())
        if u < 0 or v < 0:
            raise SystemExit("ORCA input must have nonnegative vertex IDs")
        if u != v:
            edges.add((min(u, v), max(u, v)))
if not edges:
    raise SystemExit("ORCA input graph has no edges")
n = max(max(edge) for edge in edges) + 1
with open(destination, "w", encoding="utf-8") as graph:
    graph.write(f"{n} {len(edges)}\n")
    for u, v in sorted(edges):
        graph.write(f"{u} {v}\n")
print(n)
EOF_ORCA_INPUT
    ) || exit 1

    "$RPLLB_ROOT_DIR/orca" 4 "$ORCA_INPUT" "$ORCA_OUTPUT" > smth.txt || exit 1
    # Some ORCA failures return zero: reject empty, partial or malformed output.
    if ! awk -v N="$ORBIT_NODE_COUNT" '
        NF != 15 { bad = 1 }
        { for (i = 1; i <= NF; ++i) if ($i !~ /^[0-9]+$/) bad = 1 }
        END { exit (bad || NR != N) }
    ' "$ORCA_OUTPUT"; then
        echo "ORCA did not produce $ORBIT_NODE_COUNT complete 15-column rows for $GRAPH_ABS." >&2
        exit 1
    fi
    mv -- "$ORCA_OUTPUT" "$RPLLB_ROOT_DIR/data_middle1.txt" || exit 1
)

run_bnb() {
    local K="$1"
    local GRAPH_ABS GRAPH_QUOTED
    GRAPH_ABS=$(realpath -- "$2") || return 1
    # bnb.cpp embeds argv[2] in a /bin/sh command. Pass a POSIX-quoted path
    # so spaces, apostrophes and shell metacharacters retain their literal meaning.
    GRAPH_QUOTED="'${GRAPH_ABS//\'/\'\\\'\'}'"
    (
        cd -- "$RPLLB_ROOT_DIR" || exit 1
        ./bnb "$K" "$GRAPH_QUOTED"
    ) || return 1
}

run_expfy() {
    # Measure the current graph with BLANT while retaining the initial ORCA
    # weights. Candidate directories receive a copy in prepare_serial_worker.
    run_bnb 4 "$2" > bres.tmp || return 1
    if [ ! -r data_middle1.txt ]; then
        echo "Initial node-selection weights are missing: data_middle1.txt." >&2
        return 1
    fi
    "$EXPFY_BIN" "$@"
}

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
TIGHTEST_EVB=50
TIGHTEST_DEGB=14
TIGHTEST_HDGB_REFERENCE=40
ROUND_BACKWARD_FACTOR=1.05

# RMSE_ROUNDS and RBGO are command-line parameters. After a round accepts its
# best pair, refine that same pair this many additional times.
PAIR_REFINEMENT_PASSES=3

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

# Graphlet mode is configured after command-line parsing. All scoring vectors
# have GRAPHLET_DIM coordinates: 6 for k=4 and 21 for k=5. The transformation
# sampler itself remains based on four-node graphlets.
GRAPHLET_K=""
GRAPHLET_DIM=0
declare -a GRAPHLET_IDS=()
PAIR_INPUT_WORDS_EXPECTED=0
OUTINP_WORDS_EXPECTED=0
SCALE_INPUT_WORDS_EXPECTED=0
RSEAS_WEIGHT_INCREMENT=1

# Round-level RSEAS state.
RSEAS_E=0.0
RSEAS_E_INITIALIZED=0
declare -a RSEAS_WEIGHTS=()
PREVIOUS_WORST_INDEX=""

configure_graphlet_mode() {
    case "$GRAPHLET_K" in
        4)
            GRAPHLET_DIM=6
            GRAPHLET_IDS=(5 6 7 8 9 10)
            # Preserve the current k=4 starting weights.
            RSEAS_WEIGHTS=(1 1 1 1 1 1)
            RSEAS_WEIGHT_INCREMENT=1
            ;;
        5)
            GRAPHLET_DIM=21
            # bnb 5 emits these 21 graphlet IDs in this order.
            GRAPHLET_IDS=(4 10 11 14 15 16 17 18 19 22 23 24 25 26 27 28 29 30 31 32 33)
            RSEAS_WEIGHTS=()
            for ((GRAPHLET_INDEX=0; GRAPHLET_INDEX<GRAPHLET_DIM; GRAPHLET_INDEX++)); do
                RSEAS_WEIGHTS+=(1)
            done
            RSEAS_WEIGHT_INCREMENT=3
            ;;
        *)
            echo "Graphlet size k must be 4 or 5, not '$GRAPHLET_K'." >&2
            exit 1
            ;;
    esac

    PAIR_INPUT_WORDS_EXPECTED=$((6 * GRAPHLET_DIM + 4))
    OUTINP_WORDS_EXPECTED=$((2 * GRAPHLET_DIM + 2))
    SCALE_INPUT_WORDS_EXPECTED=$((3 * GRAPHLET_DIM))
}

# Log displacement: log(A)-log(B), using LOG_FLOOR for a zero frequency.
displ() {
    python3 - "$1" "$2" "$3" "$LOG_FLOOR" "$GRAPHLET_DIM" << 'EOF_PY'
import math
import sys

path_a, path_b, output_path, log_floor_text, dimension_text = sys.argv[1:6]
log_floor = float(log_floor_text)
dimension = int(dimension_text)

A = [float(x) for x in open(path_a, encoding="utf-8").read().split()]
B = [float(x) for x in open(path_b, encoding="utf-8").read().split()]
if len(A) != dimension or len(B) != dimension:
    raise SystemExit(
        f"Expected exactly {dimension} graphlet-frequency values; got {len(A)} and {len(B)}"
    )

def safe_log(value):
    return math.log(value) if value > 0.0 else log_floor

diff = [safe_log(A[i]) - safe_log(B[i]) for i in range(dimension)]
open(output_path, "w", encoding="utf-8").write(
    " ".join(format(value, ".17g") for value in diff) + "\n"
)
EOF_PY
}

# Score two GRAPHLET_DIM-coordinate displacement vectors with one fixed RSEAS
# state. The caller supplies e followed by one positive weight per coordinate.
rseas_with_config() {
    local FILE_A="$1"
    local FILE_B="$2"
    local ERROR_RANGE="$3"
    shift 3

    if (( $# != GRAPHLET_DIM )); then
        echo "rseas_with_config requires exactly $GRAPHLET_DIM weights for k=$GRAPHLET_K." >&2
        return 1
    fi

    python3 - "$FILE_A" "$FILE_B" "$ERROR_RANGE" "$GRAPHLET_DIM" "$@" << 'EOF_PY'
import math
import sys

path_a, path_b = sys.argv[1:3]
error_range = float(sys.argv[3])
dimension = int(sys.argv[4])
weights = [float(x) for x in sys.argv[5:]]

A = [float(x) for x in open(path_a, encoding="utf-8").read().split()]
B = [float(x) for x in open(path_b, encoding="utf-8").read().split()]
if len(A) < dimension or len(B) < dimension:
    raise SystemExit(
        f"Expected at least {dimension} values in each RSEAS vector; got {len(A)} and {len(B)}"
    )
if not math.isfinite(error_range) or error_range < 0.0:
    raise SystemExit("RSEAS error range must be a finite nonnegative number")
if len(weights) != dimension or any((not math.isfinite(w) or w <= 0.0) for w in weights):
    raise SystemExit(f"RSEAS requires {dimension} positive finite weights")

distances = [abs(A[i] - B[i]) for i in range(dimension)]
maximum_distance = max(distances)

if any(distance > error_range for distance in distances):
    result = sum(
        weights[i] * max(0.0, distances[i] - error_range) ** 2
        for i in range(dimension)
    )
else:
    result = maximum_distance - error_range

print(format(result, ".17g"))
EOF_PY
}

rseas() {
    rseas_with_config "$1" "$2" "$RSEAS_E" "${RSEAS_WEIGHTS[@]}"
}

# Equal-weight RMSE used in the polishing phase. This preserves the existing
# convention sqrt(sum_i (A_i-B_i)^2), now over GRAPHLET_DIM coordinates.
unit_rmse() {
    python3 - "$1" "$2" "$GRAPHLET_DIM" << 'EOF_PY'
import math
import sys

path_a, path_b = sys.argv[1:3]
dimension = int(sys.argv[3])
A = [float(x) for x in open(path_a, encoding="utf-8").read().split()]
B = [float(x) for x in open(path_b, encoding="utf-8").read().split()]
if len(A) < dimension or len(B) < dimension:
    raise SystemExit(
        f"Expected at least {dimension} values in each RMSE vector; got {len(A)} and {len(B)}"
    )

value = sum((A[i] - B[i]) ** 2 for i in range(dimension))
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

# Print: worst_index worst_graphlet_id max_log_error candidate_error_range.
# The candidate range is 0.75*max_log_error; set_round_rseas_state then caps it
# at the preceding round's e. Ties are broken by the first emitted coordinate.
log_error_summary() {
    local GRAPHLET_ID_CSV
    GRAPHLET_ID_CSV=$(IFS=,; echo "${GRAPHLET_IDS[*]}")

    python3 - "$1" "$2" "$LOG_FLOOR" "$GRAPHLET_DIM" "$GRAPHLET_ID_CSV" << 'EOF_PY'
import math
import sys

current_path, target_path, log_floor_text, dimension_text, graphlet_ids_text = sys.argv[1:6]
log_floor = float(log_floor_text)
dimension = int(dimension_text)
graphlet_ids = [int(x) for x in graphlet_ids_text.split(",") if x]

current = [float(x) for x in open(current_path, encoding="utf-8").read().split()]
target = [float(x) for x in open(target_path, encoding="utf-8").read().split()]
if len(current) != dimension or len(target) != dimension:
    raise SystemExit(
        f"Expected exactly {dimension} graphlet frequencies; got {len(current)} and {len(target)}"
    )
if len(graphlet_ids) != dimension:
    raise SystemExit("Internal graphlet-ID list has the wrong length")

def safe_log(value):
    return math.log(value) if value > 0.0 else log_floor

distances = [
    abs(safe_log(current[i]) - safe_log(target[i]))
    for i in range(dimension)
]
worst_index = max(range(dimension), key=lambda index: distances[index])
maximum_distance = distances[worst_index]
error_range = 0.75 * maximum_distance

print(
    worst_index,
    graphlet_ids[worst_index],
    format(maximum_distance, ".17g"),
    format(error_range, ".17g"),
)
EOF_PY
}

set_round_rseas_state() {
    local ROUND_INDEX="$1"
    local CURRENT_VEC="$2"
    local TARGET_VEC="$3"
    local WORST_INDEX WORST_GRAPHLET MAX_ERROR NEW_ERROR_RANGE WEIGHT_INDEX
    local PREVIOUS_ERROR_RANGE=""
    local REPEATED_TEXT=""
    local ERROR_RANGE_TEXT=""

    read -r WORST_INDEX WORST_GRAPHLET MAX_ERROR NEW_ERROR_RANGE \
        < <(log_error_summary "$CURRENT_VEC" "$TARGET_VEC")

    WEIGHT_INDEX="$WORST_INDEX"
    if [[ -n "$PREVIOUS_WORST_INDEX" && \
          "$WORST_INDEX" == "$PREVIOUS_WORST_INDEX" ]]; then
        RSEAS_WEIGHTS[$WEIGHT_INDEX]=$((
            RSEAS_WEIGHTS[$WEIGHT_INDEX] + RSEAS_WEIGHT_INCREMENT
        ))
        REPEATED_TEXT="; repeated worst graphlet, increased its weight by $RSEAS_WEIGHT_INCREMENT"
    fi

    PREVIOUS_WORST_INDEX="$WORST_INDEX"

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
    if (( NODE_COUNT < GRAPHLET_K )); then
        echo "The cleaned target must contain at least $GRAPHLET_K nodes for k=$GRAPHLET_K." >&2
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
    # round. Keep them fixed during every appended polishing round.
    evb="$TIGHTEST_EVB"
    degb="$TIGHTEST_DEGB"
    hdgb="$TIGHTEST_HDGB"

    echo "RMSE polishing round $ROUND_INDEX/$RMSE_ROUNDS constraints: evb=$evb, degb=$degb, hdgb=$hdgb"
}

INITIAL_GRAPH_GIVEN=0

if [ $# -eq 6 ]; then
    GRAPHLET_K="$1"
    TARGET0="$2"
    FINAL="$3"
    ROUNDS="$4"
    RMSE_ROUNDS="$5"
    RBGO="$6"
    PRE_INIT_SYNTH=""
elif [ $# -eq 7 ]; then
    GRAPHLET_K="$1"
    TARGET0="$2"
    FINAL="$3"
    PRE_INIT_SYNTH="$4"
    ROUNDS="$5"
    RMSE_ROUNDS="$6"
    RBGO="$7"
    INITIAL_GRAPH_GIVEN=1
else
    echo "Usage: $0 <k:4|5> <Target> <Output> <RSEAS_Rounds> <RMSE_Rounds> <RBGO>"
    echo "   or: $0 <k:4|5> <Target> <Output> <Base_Synth> <RSEAS_Rounds> <RMSE_Rounds> <RBGO>"
    exit 1
fi
echo "$TARGET0" >> wcic.txt

if ! [[ "$GRAPHLET_K" =~ ^[45]$ ]]; then
    echo "Graphlet size k must be 4 or 5, not '$GRAPHLET_K'." >&2
    exit 1
fi
if ! [[ "$ROUNDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "RSEAS rounds must be a positive integer, not '$ROUNDS'." >&2
    exit 1
fi
if ! [[ "$RMSE_ROUNDS" =~ ^[0-9]+$ ]]; then
    echo "RMSE rounds must be a nonnegative integer, not '$RMSE_ROUNDS'." >&2
    exit 1
fi
if ! [[ "$RBGO" =~ ^[1-9][0-9]*$ ]]; then
    echo "RBGO must be a positive integer, not '$RBGO'." >&2
    exit 1
fi

configure_graphlet_mode

echo "Optimizing k=$GRAPHLET_K graphlets with $GRAPHLET_DIM frequency coordinates."
echo "Schedule: $ROUNDS RSEAS rounds, $RMSE_ROUNDS RMSE rounds, regenerate transformation probes every $RBGO round(s) within each phase."

TARGET="target_c.txt"
./cleanup "$TARGET0" "$TARGET"
configure_node_scaled_parameters
build_round_constraint_schedule

if (( ! INITIAL_GRAPH_GIVEN )); then
    ./gen_deg2 "$TARGET" synth_base_1.txt
    initialize_node_weights synth_base_1.txt || exit 1

    run_expfy "$TARGET" synth_base_1.txt "$INIT_PRIMARY_SAMPLES" "$INIT_PRIMARY_SAMPLES" 1 8 2.0 20 "$INIT_HDGB_1" > synth_base_2.txt
    run_expfy "$TARGET" synth_base_2.txt "$INIT_MAIN_SAMPLES" "$INIT_MAIN_SAMPLES" 1 8 1.7 18 "$INIT_HDGB_2" > synth_base_3.txt
    run_expfy "$TARGET" synth_base_3.txt "$INIT_MAIN_SAMPLES" "$INIT_MAIN_SAMPLES" 1 8 1.5 18 "$INIT_HDGB_3" > synth_base_4.txt
    run_expfy "$TARGET" synth_base_4.txt "$INIT_MAIN_SAMPLES" "$INIT_MAIN_SAMPLES" 1 8 1.3 17 "$INIT_HDGB_4" > synth_base_5.txt
    run_expfy "$TARGET" synth_base_5.txt "$INIT_MAIN_SAMPLES" "$INIT_SECONDARY_SAMPLES" 1 17 1.2 16 "$INIT_HDGB_5" > synth_base_6.txt
    run_expfy "$TARGET" synth_base_6.txt "$INIT_MAIN_SAMPLES" "$INIT_SECONDARY_SAMPLES" 1 17 1.1 15 "$INIT_HDGB_6" > synth_base_7.txt
    run_expfy "$TARGET" synth_base_7.txt "$INIT_MAIN_SAMPLES" "$INIT_SECONDARY_SAMPLES" 1 17 1 14 "$INIT_HDGB_7" > synth_base_8.txt

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
    initialize_node_weights "$PRE_INIT_SYNTH" || exit 1
fi

# PIS uses one target vector, unit weights, and one error range fixed from the
# graph entering PIS.  The state is initialized for both generated and supplied
# base graphs because evc1 and evc3 must be compared under the same objective.
PIS_TARGET_VEC="pis_target_vec.tmp"
PIS_CURRENT_RSEAS=""
PIS_RSEAS_E=0.0
declare -a PIS_RSEAS_WEIGHTS=()
for ((GRAPHLET_INDEX=0; GRAPHLET_INDEX<GRAPHLET_DIM; GRAPHLET_INDEX++)); do
    PIS_RSEAS_WEIGHTS+=(1)
done

pis_graph_rseas() {
    local GRAPH="$1"
    local VEC_FILE DISPL_FILE SCORE STATUS

    VEC_FILE=$(mktemp ".pis_vec.XXXXXX")
    DISPL_FILE=$(mktemp ".pis_displ.XXXXXX")
    if ! run_bnb "$GRAPHLET_K" "$GRAPH" > "$VEC_FILE"; then
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

    # The threshold branch can deliberately produce an exact copy. Do not
    # resample it and accidentally perturb the remembered score.
    if cmp -s "$BEFORE" "$CANDIDATE"; then
        mv "$CANDIDATE" "$OUTPUT"
        echo "$LABEL made no change; RSEAS remains $PIS_CURRENT_RSEAS."
        return
    fi

    local CANDIDATE_RSEAS
    CANDIDATE_RSEAS=$(pis_graph_rseas "$CANDIDATE")

    echo "$LABEL RSEAS: before=$PIS_CURRENT_RSEAS, candidate=$CANDIDATE_RSEAS"

    # Preserve the existing behavior for an automatically generated base: the
    # PIS transformation is applied directly.  For a user-supplied base, reject
    # a candidate that worsens the fixed PIS objective.
    if (( ! INITIAL_GRAPH_GIVEN )); then
        mv "$CANDIDATE" "$OUTPUT"
        PIS_CURRENT_RSEAS="$CANDIDATE_RSEAS"
        echo "$LABEL accepted."
    elif awk -v NEW="$CANDIDATE_RSEAS" -v OLD="$PIS_CURRENT_RSEAS" \
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

# When the current eigenvalue is too small, both evc1 and evc2 move it in the
# required direction. Generate both from the same input using the same sample
# count, then pass the lower-RSEAS candidate to the normal PIS acceptance rule.
make_pis_increase_candidate() {
    local INPUT_GRAPH="$1"
    local NUM_SAMPLES="$2"
    local LABEL="$3"
    local OUTPUT_CANDIDATE="$4"

    local EVC1_CANDIDATE="${LABEL}_evc1_candidate.txt"
    local EVC2_CANDIDATE="${LABEL}_evc2_candidate.txt"
    local EVC1_SCORE EVC2_SCORE

    ./evc1 "$INPUT_GRAPH" "$NUM_SAMPLES" > "$EVC1_CANDIDATE"
    ./evc2 "$INPUT_GRAPH" "$NUM_SAMPLES" > "$EVC2_CANDIDATE"

    EVC1_SCORE=$(pis_graph_rseas "$EVC1_CANDIDATE")
    EVC2_SCORE=$(pis_graph_rseas "$EVC2_CANDIDATE")

    echo "$LABEL eigenvalue-increase candidates: evc1 RSEAS=$EVC1_SCORE, evc2 RSEAS=$EVC2_SCORE"

    # Prefer evc1 on an exact tie so the old PIS behavior remains the stable
    # tie-breaker. evc3 is selected only when its measured RSEAS is lower.
    if awk -v EVC2="$EVC2_SCORE" -v EVC1="$EVC1_SCORE" \
        'BEGIN { exit !(EVC2 < EVC1) }'; then
        mv "$EVC2_CANDIDATE" "$OUTPUT_CANDIDATE"
        rm -f "$EVC1_CANDIDATE"
        echo "$LABEL selected evc2."
    else
        mv "$EVC1_CANDIDATE" "$OUTPUT_CANDIDATE"
        rm -f "$EVC2_CANDIDATE"
        echo "$LABEL selected evc1."
    fi
}

# Initialize the fixed PIS objective for every run. This is needed even for a
# generated base so evc1 and evc3 can be compared fairly.
local_pis_initial_vec=$(mktemp ".pis_initial_vec.XXXXXX")
local_pis_initial_displ=$(mktemp ".pis_initial_displ.XXXXXX")

run_bnb "$GRAPHLET_K" "$TARGET" > "$PIS_TARGET_VEC"
run_bnb "$GRAPHLET_K" "$PRE_INIT_SYNTH" > "$local_pis_initial_vec"

read -r PIS_WORST_INDEX PIS_WORST_GRAPHLET PIS_MAX_ERROR PIS_RSEAS_E \
    < <(log_error_summary "$local_pis_initial_vec" "$PIS_TARGET_VEC")

displ "$local_pis_initial_vec" "$PIS_TARGET_VEC" "$local_pis_initial_displ"
PIS_CURRENT_RSEAS=$(rseas_with_config f0.txt "$local_pis_initial_displ" \
    "$PIS_RSEAS_E" "${PIS_RSEAS_WEIGHTS[@]}")

rm -f "$local_pis_initial_vec" "$local_pis_initial_displ"

echo "PIS RSEAS state: worst graphlet=$PIS_WORST_GRAPHLET, max log error=$PIS_MAX_ERROR, e=$PIS_RSEAS_E, weights=${PIS_RSEAS_WEIGHTS[*]}"
echo "Initial graph RSEAS before PIS: $PIS_CURRENT_RSEAS"

PIS1_CANDIDATE="PIS1_candidate.txt"
./evs "$PRE_INIT_SYNTH" > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
if float_below_by "$EV1" "$EV2" 14; then
    make_pis_increase_candidate "$PRE_INIT_SYNTH" "$PIS1_SAMPLES" PIS1 "$PIS1_CANDIDATE"
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
    make_pis_increase_candidate PIS1.txt "$PIS2_SAMPLES" PIS2 "$PIS2_CANDIDATE"
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
    make_pis_increase_candidate PIS2.txt "$PIS3_SAMPLES" PIS3 "$PIS3_CANDIDATE"
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
MAX_TRANSFORMATION=95
CANDIDATE_LIMIT="$MAX_TRANSFORMATION"
if (( ${EXPFY_MODE_LIMIT:-MAX_TRANSFORMATION} < MAX_TRANSFORMATION )); then
    CANDIDATE_LIMIT="$EXPFY_MODE_LIMIT"
    echo "Selected backend supports modes through $CANDIDATE_LIMIT; later modes are excluded." >&2
fi
CANDS=($(seq 1 "$CANDIDATE_LIMIT" | grep -v -E '^(25|26)$'))
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
# New G6/G7-directed moves. Mode 64 starts from a G6 leaf slide and also
# requires a disjoint pendant triangle, checked internally by expfy.
initial_graphlet[61]=6
initial_graphlet[62]=7
initial_graphlet[63]=7
initial_graphlet[64]=6
# Edge-neutral G10 reduction: all three rewrites start from a clique.
initial_graphlet[65]=10
initial_graphlet[66]=10
initial_graphlet[67]=10
initial_graphlet[68]=6  # Triangle-guided P4 leaf relocation.
initial_graphlet[69]=6  # Degree-histogram-preserving P4 reduction.
initial_graphlet[70]=6  # Degree-preserving leaf/branch switch creating boxes.
initial_graphlet[71]=7  # Open an isolated triangle into boxes using a leaf swap.
initial_graphlet[72]=7  # Open a triangle with a branch exchange, keeping hub degrees fixed.
initial_graphlet[73]=6  # Inverse T20: path -> star.
initial_graphlet[74]=5  # Inverse T9: star -> path.
initial_graphlet[75]=6  # Inverse T4: path -> path.
initial_graphlet[76]=7  # Inverse T22: paw -> cycle.
initial_graphlet[77]=7  # Inverse T23: paw -> cycle.
initial_graphlet[78]=6  # Inverse T1: path -> path.
initial_graphlet[79]=5  # Inverse T8: star -> path.
initial_graphlet[80]=8  # Inverse T16: cycle -> paw.
initial_graphlet[81]=7  # Inverse T13: paw -> paw.
initial_graphlet[82]=7  # Inverse T10: paw -> paw.
initial_graphlet[83]=7  # Inverse T14: paw -> paw.
initial_graphlet[84]=6  # Inverse T2: path -> path.
initial_graphlet[85]=6  # Inverse T3: path -> path.
initial_graphlet[86]=6  # Inverse T5: path -> path.
initial_graphlet[87]=6  # Inverse T6: path -> path.
initial_graphlet[88]=5  # Inverse T7: star -> path.
initial_graphlet[89]=7  # Inverse T11: paw -> paw.
initial_graphlet[90]=7  # Inverse T12: paw -> paw.
initial_graphlet[91]=7  # Inverse T15: paw -> paw.
initial_graphlet[92]=8  # Inverse T17: cycle -> paw.
initial_graphlet[93]=5  # Inverse T18: star -> star.
initial_graphlet[94]=6  # Inverse T19: path -> star.
initial_graphlet[95]=8  # Inverse T21: cycle -> cycle.

# data_middle1.txt column index -> graphlet ID.  Several columns can describe
# different node roles in the same graphlet.
column_graphlet=(0 0 0 0 6 6 5 5 8 7 7 7 9 9 10)

declare -a graphlet_node_count
declare -a rare_graphlet

# Transformation probes and round candidates are evaluated sequentially.

declare -a ec
for ((i=1; i<=MAX_TRANSFORMATION; i++)); do
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
ec[61]=0;   # Guarded P4 leaf slide.
ec[62]=0;   # Guarded paw-tail relocation.
ec[63]=-1;  # Pendant triangle opening.
ec[64]=-1;  # Guarded P4 leaf slide plus pendant triangle opening.
ec[65]=0;   # Exchange a clique edge across a shared triangle.
ec[66]=0;   # Switch clique/box edges; every degree stays unchanged.
ec[67]=0;   # Move a clique edge into an adjoining box diagonal.
ec[68]=0;   # Reattach a true leaf; one edge removed and one added.
ec[69]=0;   # Reattach a true leaf while exchanging its roots' degrees.
ec[70]=0;   # Swap two edges, preserving every node's degree.
ec[71]=0;   # Open a paw into a box and move a leaf; every degree stays fixed.
ec[72]=0;   # Exchange two edges between a triangle and a branch.
# Modes 73..95 retain the default ec=0: all learned inverses preserve edge count.

BASELINE_DIR="baseline"
mkdir -p "$BASELINE_DIR"

# Keep candidate scratch files in private directories even though work is serial.
# BLANT-backed bnb calls are routed to the launch directory by run_bnb.
prepare_serial_worker() {
    local WORK_DIR="$1"
    local ROOT_DIR="$2"
    local SHARED_FILE TOOL

    mkdir -p "$WORK_DIR"

    for SHARED_FILE in data_middle1.txt f0.txt; do
        if [ -f "$ROOT_DIR/$SHARED_FILE" ]; then
            cp "$ROOT_DIR/$SHARED_FILE" "$WORK_DIR/$SHARED_FILE"
        fi
    done

    # Other tools can still use project-local executables in this directory.
    # run_bnb explicitly returns to RPLLB_ROOT_DIR for BLANT and its resources.
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

# Read the output of `evs graph > file`. The gate relies on evs writing one
# finite floating-point value and nothing else to stdout.
read_single_eigenvalue() {
    local FILE="$1"
    local LABEL="$2"

    python3 - "$FILE" "$LABEL" << 'EOF_PY'
import math
import sys

path, label = sys.argv[1:3]
values = open(path, "r", encoding="utf-8").read().split()
if len(values) != 1:
    raise SystemExit(
        f"{label} eigenvalue output must contain exactly one value; "
        f"found {len(values)} in {path}"
    )

try:
    value = float(values[0])
except ValueError as error:
    raise SystemExit(f"Invalid {label} eigenvalue {values[0]!r}: {error}")

if not math.isfinite(value):
    raise SystemExit(f"{label} eigenvalue must be finite, got {value}")

print(format(value, ".17g"))
EOF_PY
}

# Print:
#   status|after_log_distance|current_log_distance|near_target_limit|reason
#
# An EVC candidate is allowed precisely when
#
#   |log(target_ev)-log(|after_ev|)| < 0.2*evb
#
# or when it is strictly closer to the target eigenvalue than the current
# graph. The absolute value applies to after_ev exactly as specified.
evc_eigenvalue_gate() {
    local TARGET_EV="$1"
    local CURRENT_EV="$2"
    local AFTER_EV="$3"
    local EVB="$4"

    python3 - "$TARGET_EV" "$CURRENT_EV" "$AFTER_EV" "$EVB" << 'EOF_PY'
import math
import sys

target_ev, current_ev, after_ev, evb = map(float, sys.argv[1:5])

if not all(math.isfinite(value) for value in (target_ev, current_ev, after_ev, evb)):
    raise SystemExit("Eigenvalue gate received a non-finite value")
if target_ev <= 0.0:
    raise SystemExit(f"Target largest eigenvalue must be positive, got {target_ev}")
if current_ev <= 0.0:
    raise SystemExit(f"Current largest eigenvalue must be positive, got {current_ev}")
if evb < 0.0:
    raise SystemExit(f"evb must be nonnegative, got {evb}")

after_abs = abs(after_ev)
current_distance = abs(math.log(target_ev) - math.log(current_ev))
near_target_limit = 0.002 * evb

if after_abs <= 0.0:
    print(
        "DISALLOWED|inf|"
        f"{current_distance:.17g}|{near_target_limit:.17g}|nonpositive_after_ev"
    )
    raise SystemExit(0)

after_distance = abs(math.log(target_ev) - math.log(after_abs))

if after_distance < near_target_limit:
    status = "ALLOWED"
    reason = "within_0.2_evb"
elif after_distance < current_distance:
    status = "ALLOWED"
    reason = "closer_than_current"
else:
    status = "DISALLOWED"
    reason = "too_far_and_not_closer"

print(
    f"{status}|{after_distance:.17g}|{current_distance:.17g}|"
    f"{near_target_limit:.17g}|{reason}"
)
EOF_PY
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
            # Mode 60 uses five nodes rather than a single four-node graphlet.
            # It cannot be filtered using one four-node graphlet column.
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
    if ! run_expfy "$TARGET" "$INPUT_FILE" "$REQUESTED_SAMPLES" 0 \
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
    if ! run_expfy "$TARGET" "$INPUT_FILE" \
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


        local INPUT_ABS OUTPUT_ABS TARGET_ABS VEC_TAR_ABS F0_ABS
        local EVS_ABS BXK4F_ABS BXK4ONE_ABS
        INPUT_ABS=$(realpath "$INPUT_SYNTH")
        OUTPUT_ABS=$(realpath -m "$OUTPUT_SYNTH")
        TARGET_ABS=$(realpath "$TARGET")
        VEC_TAR_ABS=$(realpath vecTar.txt)
        F0_ABS=$(realpath f0.txt)
        EVS_ABS=$(realpath ./evs)
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

        run_bnb "$GRAPHLET_K" "$INPUT_ABS" > "$VEC_INS"

        displ "$VEC_INS" "$VEC_TAR_ABS" "$INPL"
        cat "$VEC_INS" "$VEC_TAR_ABS" > "$INP"

        local had_score
        had_score=$(objective_score "$F0_ABS" "$INPL")

        echo "New $ROUND_OBJECTIVE round with input $INPUT_SYNTH"
        echo "Evaluating round candidates sequentially"

        # ---------------------------------------------------------
        # 1. Generate and evaluate six EVC candidates sequentially: small and
        #    large runs of evc1, evc2, and the union/intersection evc3.
        # ---------------------------------------------------------
        local EVC_DIR="$ROUND_DIR/evc"
        mkdir -p "$EVC_DIR"

        # Compute the fixed target/current reference once for this round. Each
        # EVC worker computes only its own after-change eigenvalue.
        local TARGET_EV_FILE="$EVC_DIR/target.ev"
        local CURRENT_EV_FILE="$EVC_DIR/current.ev"
        local TARGET_EV CURRENT_EV EVC_REFERENCE
        local EVC_REFERENCE_STATUS EVC_REFERENCE_AFTER_DISTANCE
        local CURRENT_EV_DISTANCE EVC_DISTANCE_LIMIT EVC_REFERENCE_REASON

        "$EVS_ABS" "$TARGET_ABS" > "$TARGET_EV_FILE"
        "$EVS_ABS" "$INPUT_ABS" > "$CURRENT_EV_FILE"
        TARGET_EV=$(read_single_eigenvalue "$TARGET_EV_FILE" "target")
        CURRENT_EV=$(read_single_eigenvalue "$CURRENT_EV_FILE" "current")

        # Reuse the candidate-gate calculation to obtain the round's current
        # log-distance and its absolute 0.2*evb threshold.
        EVC_REFERENCE=$(evc_eigenvalue_gate \
            "$TARGET_EV" "$CURRENT_EV" "$CURRENT_EV" "$evb")
        IFS='|' read -r EVC_REFERENCE_STATUS EVC_REFERENCE_AFTER_DISTANCE \
            CURRENT_EV_DISTANCE EVC_DISTANCE_LIMIT EVC_REFERENCE_REASON \
            <<< "$EVC_REFERENCE"

        echo "EVC eigenvalue gate: target_ev=$TARGET_EV, current_ev=$CURRENT_EV, current_log_distance=$CURRENT_EV_DISTANCE, near_target_limit=$EVC_DISTANCE_LIMIT"

        local -a EVC_LABELS=(
            outpev11 outpev12
            outpev21 outpev22
            outpev31 outpev32
        )
        local -a EVC_TOOLS=(
            evc1 evc1
            evc2 evc2
            evc3 evc3
        )
        local -a EVC_AMOUNTS=(
            "$OUTPEV_SMALL_SAMPLES" "$OUTPEV_LARGE_SAMPLES"
            "$OUTPEV_SMALL_SAMPLES" "$OUTPEV_LARGE_SAMPLES"
            "$OUTPEV_SMALL_SAMPLES" "$OUTPEV_LARGE_SAMPLES"
        )
        local IDX LABEL TOOL AMOUNT

        for IDX in "${!EVC_LABELS[@]}"; do
            LABEL="${EVC_LABELS[$IDX]}"
            TOOL="${EVC_TOOLS[$IDX]}"
            AMOUNT="${EVC_AMOUNTS[$IDX]}"

            (
                set -e
                local WORK_DIR="$EVC_DIR/work_$LABEL"
                prepare_serial_worker "$WORK_DIR" "$ROOT_DIR"
                trap 'rm -rf -- "$WORK_DIR"' EXIT
                cd "$WORK_DIR"

                "./$TOOL" "$INPUT_ABS" "$AMOUNT" > candidate.graph

                local AFTER_EV GATE_RESULT GATE_STATUS
                local AFTER_EV_DISTANCE WORKER_CURRENT_DISTANCE
                local WORKER_DISTANCE_LIMIT GATE_REASON

                "$EVS_ABS" candidate.graph > candidate.ev
                AFTER_EV=$(read_single_eigenvalue candidate.ev "$LABEL after-change")
                GATE_RESULT=$(evc_eigenvalue_gate \
                    "$TARGET_EV" "$CURRENT_EV" "$AFTER_EV" "$evb")
                IFS='|' read -r GATE_STATUS AFTER_EV_DISTANCE \
                    WORKER_CURRENT_DISTANCE WORKER_DISTANCE_LIMIT GATE_REASON \
                    <<< "$GATE_RESULT"

                if [ "$GATE_STATUS" = "DISALLOWED" ]; then
                    # A rejected EVC graph is not a failed worker. Publish the
                    # rejection and skip the much more expensive bnb call.
                    printf 'DISALLOWED|NA|NA|%s|%s|%s\n' \
                        "$AFTER_EV" "$AFTER_EV_DISTANCE" "$GATE_REASON" \
                        > "$EVC_DIR/$LABEL.result.tmp"
                    mv "$EVC_DIR/$LABEL.result.tmp" "$EVC_DIR/$LABEL.result"
                    exit 0
                fi

                if [ "$GATE_STATUS" != "ALLOWED" ]; then
                    echo "Unexpected EVC eigenvalue-gate status '$GATE_STATUS' for $LABEL." >&2
                    exit 1
                fi

                run_bnb "$GRAPHLET_K" candidate.graph > candidate.vec
                displ "$VEC_INS" candidate.vec candidate.displ

                local SCORE
                SCORE=$(objective_score candidate.displ "$INPL")

                mv candidate.graph "$EVC_DIR/$LABEL.graph"
                printf 'ALLOWED|%s|%s|%s|%s|%s\n' \
                    "$SCORE" "$EVC_DIR/$LABEL.graph" "$AFTER_EV" \
                    "$AFTER_EV_DISTANCE" "$GATE_REASON" \
                    > "$EVC_DIR/$LABEL.result.tmp"
                mv "$EVC_DIR/$LABEL.result.tmp" "$EVC_DIR/$LABEL.result"
            )
        done

        local BEST_EVC_GRAPH="$INPUT_ABS"
        local BEST_EVC_SCORE="$had_score"
        local EVC_STATUS CUR_SCORE CUR_GRAPH CUR_EV CUR_EV_DISTANCE CUR_REASON

        # Read in the original deterministic order so equal scores are handled
        # exactly as they were in the serial version.
        for LABEL in "${EVC_LABELS[@]}"; do
            IFS='|' read -r EVC_STATUS CUR_SCORE CUR_GRAPH CUR_EV \
                CUR_EV_DISTANCE CUR_REASON < "$EVC_DIR/$LABEL.result"

            case "$EVC_STATUS" in
                ALLOWED)
                    echo "$LABEL eigenvalue allowed ($CUR_REASON): after_ev=$CUR_EV, log_distance=$CUR_EV_DISTANCE; $ROUND_OBJECTIVE=$CUR_SCORE"
                    if float_less "$CUR_SCORE" "$BEST_EVC_SCORE"; then
                        BEST_EVC_SCORE="$CUR_SCORE"
                        BEST_EVC_GRAPH="$CUR_GRAPH"
                    fi
                    ;;
                DISALLOWED)
                    echo "$LABEL eigenvalue disallowed: after_ev=$CUR_EV, log_distance=$CUR_EV_DISTANCE; requires distance < $EVC_DISTANCE_LIMIT or distance < $CURRENT_EV_DISTANCE"
                    ;;
                *)
                    echo "Invalid EVC result status '$EVC_STATUS' for $LABEL." >&2
                    exit 1
                    ;;
            esac
        done

        if [ "$BEST_EVC_GRAPH" != "$INPUT_ABS" ]; then
            echo "Replacing $INPUT_SYNTH with the best EVC candidate ($ROUND_OBJECTIVE $BEST_EVC_SCORE)"
            cp "$BEST_EVC_GRAPH" "$INPUT_ABS"
            had_score="$BEST_EVC_SCORE"
        else
            echo "No EVC candidate improved the input graph."
        fi

        # Recompute the current graph vector after the possible EVC replacement.
        run_bnb "$GRAPHLET_K" "$INPUT_ABS" > "$VEC_INS"
        displ "$VEC_INS" "$VEC_TAR_ABS" "$INPL"
        cat "$VEC_INS" "$VEC_TAR_ABS" > "$INP"

        # ---------------------------------------------------------
        # 2. Read all transformation probes once and score every active pair
        #    in one bxk4f invocation. No shell worker or process per pair.
        # ---------------------------------------------------------
        local RESULTS_FILE="$ROUND_DIR/round_results.tmp"
        local BATCH_INPUT="$ROUND_DIR/bxk4f_batch_input.tmp"
        local -a BATCH_ENABLED=()
        local BATCH_ID BATCH_PROBE BATCH_PROBE_WORDS
        local CAND_COUNT=${#ACTIVE_CANDS[@]}

        for BATCH_ID in "${ACTIVE_CANDS[@]}"; do
            BATCH_ENABLED[$BATCH_ID]=1
        done

        # Every ID through MAX_TRANSFORMATION has edge-change metadata. Disabled
        # entries (unused IDs 25/26 or rarity-filtered modes) have no curve.
        # Current/target are raw frequencies; A/B are the measured log
        # displacements already stored by generate_outinp.
        {
            printf 'BXK4F_BATCH_V1\n%d\n' "$MAX_TRANSFORMATION"
            cat "$INP"
            printf '\n'
            for ((BATCH_ID=1; BATCH_ID<=MAX_TRANSFORMATION; BATCH_ID++)); do
                printf '%d %d %d\n' "$BATCH_ID" "${ec[$BATCH_ID]}" \
                    "${BATCH_ENABLED[$BATCH_ID]:-0}"
                if (( ${BATCH_ENABLED[$BATCH_ID]:-0} )); then
                    BATCH_PROBE="$ROOT_DIR/outinp_val_${BATCH_ID}.txt"
                    if [ ! -s "$BATCH_PROBE" ]; then
                        echo "Missing transformation data: $BATCH_PROBE" >&2
                        exit 1
                    fi
                    BATCH_PROBE_WORDS=$(wc -w < "$BATCH_PROBE")
                    if (( BATCH_PROBE_WORDS != OUTINP_WORDS_EXPECTED )); then
                        echo "Malformed $BATCH_PROBE: expected $OUTINP_WORDS_EXPECTED values for k=$GRAPHLET_K, got $BATCH_PROBE_WORDS." >&2
                        exit 1
                    fi
                    cat "$BATCH_PROBE"
                    printf '\n'
                fi
            done
        } > "$BATCH_INPUT"

        echo "Scoring $((CAND_COUNT * (CAND_COUNT + 1) / 2)) pairs from $CAND_COUNT active transformations in one bxk4f call."
        if [ "$ROUND_OBJECTIVE" = "RSEAS" ]; then
            "$BXK4F_ABS" "$GRAPHLET_K" "$RSEAS_E" "${RSEAS_WEIGHTS[@]}" \
                < "$BATCH_INPUT" > "$RESULTS_FILE.pending"
        else
            "$BXK4F_ABS" "$GRAPHLET_K" \
                < "$BATCH_INPUT" > "$RESULTS_FILE.pending"
        fi
        mv "$RESULTS_FILE.pending" "$RESULTS_FILE"

        had_score=$(objective_score "$F0_ABS" "$INPL")
        echo "$ROUND_OBJECTIVE before modification is $had_score"

        local BAD_ZERO_SCALE_COUNT
        BAD_ZERO_SCALE_COUNT=$(awk -F'|' '
            $1 == "BAD_ZERO_SCALE" { count++ }
            END { print count + 0 }
        ' "$RESULTS_FILE")
        if (( BAD_ZERO_SCALE_COUNT > 0 )); then
            echo "Excluded $BAD_ZERO_SCALE_COUNT mixed pairs because scale=0 implies a redundant candidate; self-pairs were retained."
        fi

        # ---------------------------------------------------------
        # 3. Evaluate the ten most promising nonredundant pairs sequentially.
        # ---------------------------------------------------------
        local TOP_DIR="$ROUND_DIR/top"
        mkdir -p "$TOP_DIR"

        local -a TOP_LINES=()
        mapfile -t TOP_LINES < <(
            awk -F'|' '$1 != "BAD_ZERO_SCALE"' "$RESULTS_FILE" |
                sort -n |
                head -n 10
        )

        if (( ${#TOP_LINES[@]} == 0 )); then
            echo "No pair-scoring results were produced." >&2
            exit 1
        fi

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
                prepare_serial_worker "$WORK_DIR" "$ROOT_DIR"
                trap 'rm -rf -- "$WORK_DIR"' EXIT
                cd "$WORK_DIR"

                # apply_exp2 refers to TARGET dynamically and invokes run_expfy.
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

                run_bnb "$GRAPHLET_K" "$STAGE1" > stage1.vec
                displ "$VEC_INS" stage1.vec stage1.displ
                cat "$INP" stage1.displ > scale_input.tmp
                if [ "$ROUND_OBJECTIVE" = "RSEAS" ]; then
                    "$BXK4ONE_ABS" "$GRAPHLET_K" "$RSEAS_E" "${RSEAS_WEIGHTS[@]}" \
                        < scale_input.tmp > scale_output.tmp
                else
                    "$BXK4ONE_ABS" "$GRAPHLET_K" < scale_input.tmp > scale_output.tmp
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

                run_bnb "$GRAPHLET_K" "$FINAL_CANDIDATE" > final.vec
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
            )
        done

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
		echo "$BEST_VAL $BEST_VAL2" >> wcic.txt

        if float_greater "$BEST_SCORE" "$had_score"; then
            echo "The best evaluated pair was worse than the input ($ROUND_OBJECTIVE $had_score); keeping the input graph."
            cp "$INPUT_ABS" "$OUTPUT_ABS"
        else
            echo "Using the saved best candidate, then refining the same pair $PAIR_REFINEMENT_PASSES times."

            local REFINE_DIR="$ROUND_DIR/refinement"
            mkdir -p "$REFINE_DIR"

            local REFINE_GRAPH="$REFINE_DIR/current_best.graph"
            local REFINE_SCORE="$BEST_SCORE"
            local REFINE_SCALE="$BEST_MAG"
            local REFINE_SCALE2="$BEST_MAG2"
            cp "$BEST_GRAPH" "$REFINE_GRAPH"

            local REFINE_ITER
            local CURRENT_VEC PROBE_VEC CANDIDATE_VEC
            local PROBE_DISPL CANDIDATE_TARGET_DISPL SCALE_INPUT
            local PROBE_GRAPH REFINED_CANDIDATE
            local PROBE_TESTED_SCALE PROBE_TESTED_SCALE2
            local REQUESTED_NEW_SCALE REQUESTED_NEW_SCALE2
            local CANDIDATE_TESTED_SCALE CANDIDATE_TESTED_SCALE2
            local MAG_FACTOR CANDIDATE_SCORE SCALE_INPUT_WORDS

            for ((REFINE_ITER=1; REFINE_ITER<=PAIR_REFINEMENT_PASSES; REFINE_ITER++)); do
                CURRENT_VEC="$REFINE_DIR/current_${REFINE_ITER}.vec"
                PROBE_VEC="$REFINE_DIR/probe_${REFINE_ITER}.vec"
                CANDIDATE_VEC="$REFINE_DIR/candidate_${REFINE_ITER}.vec"
                PROBE_DISPL="$REFINE_DIR/probe_${REFINE_ITER}.displ"
                CANDIDATE_TARGET_DISPL="$REFINE_DIR/candidate_target_${REFINE_ITER}.displ"
                SCALE_INPUT="$REFINE_DIR/scale_input_${REFINE_ITER}.tmp"
                PROBE_GRAPH="$REFINE_DIR/probe_${REFINE_ITER}.graph"
                REFINED_CANDIDATE="$REFINE_DIR/refined_candidate_${REFINE_ITER}.graph"

                # Refresh the objective vector; node-selection weights remain
                # those of the initial synthetic graph, including in k=5 runs.
                run_bnb "$GRAPHLET_K" "$REFINE_GRAPH" > "$CURRENT_VEC"

                # Probe by applying the same pair once more at its latest
                # accepted effective scales.
                apply_exp2 "$REFINE_GRAPH" "$REFINE_SCALE" "$REFINE_SCALE2" \
                    "$BEST_VAL" "$BEST_VAL2" "$PROBE_GRAPH" \
                    PROBE_TESTED_SCALE PROBE_TESTED_SCALE2

                run_bnb "$GRAPHLET_K" "$PROBE_GRAPH" > "$PROBE_VEC"
                displ "$CURRENT_VEC" "$PROBE_VEC" "$PROBE_DISPL"
                cat "$CURRENT_VEC" "$VEC_TAR_ABS" "$PROBE_DISPL" > "$SCALE_INPUT"

                SCALE_INPUT_WORDS=$(wc -w < "$SCALE_INPUT")
                if (( SCALE_INPUT_WORDS != SCALE_INPUT_WORDS_EXPECTED )); then
                    echo "Malformed refinement bxk4one input: expected $SCALE_INPUT_WORDS_EXPECTED values for k=$GRAPHLET_K, got $SCALE_INPUT_WORDS." >&2
                    exit 1
                fi

                if [ "$ROUND_OBJECTIVE" = "RSEAS" ]; then
                    "$BXK4ONE_ABS" "$GRAPHLET_K" "$RSEAS_E" "${RSEAS_WEIGHTS[@]}" \
                        < "$SCALE_INPUT" > "$REFINE_DIR/scale_output_${REFINE_ITER}.tmp"
                else
                    "$BXK4ONE_ABS" "$GRAPHLET_K" \
                        < "$SCALE_INPUT" > "$REFINE_DIR/scale_output_${REFINE_ITER}.tmp"
                fi
                read -r MAG_FACTOR < "$REFINE_DIR/scale_output_${REFINE_ITER}.tmp"

                # bxk4one returns one multiplier for the combined pair response,
                # so x' and y' are the two actually completed probe scales times
                # that common multiplier.
                REQUESTED_NEW_SCALE=$(awk -v S="$PROBE_TESTED_SCALE" -v M="$MAG_FACTOR" \
                    'BEGIN { print int(S * M) }')
                REQUESTED_NEW_SCALE2=$(awk -v S="$PROBE_TESTED_SCALE2" -v M="$MAG_FACTOR" \
                    'BEGIN { print int(S * M) }')

                apply_exp2 "$REFINE_GRAPH" "$REQUESTED_NEW_SCALE" "$REQUESTED_NEW_SCALE2" \
                    "$BEST_VAL" "$BEST_VAL2" "$REFINED_CANDIDATE" \
                    CANDIDATE_TESTED_SCALE CANDIDATE_TESTED_SCALE2

                run_bnb "$GRAPHLET_K" "$REFINED_CANDIDATE" > "$CANDIDATE_VEC"
                displ "$CANDIDATE_VEC" "$VEC_TAR_ABS" "$CANDIDATE_TARGET_DISPL"
                CANDIDATE_SCORE=$(objective_score "$F0_ABS" "$CANDIDATE_TARGET_DISPL")

                echo "Pair refinement $REFINE_ITER/$PAIR_REFINEMENT_PASSES for $BEST_VAL,$BEST_VAL2: probe scales=$PROBE_TESTED_SCALE,$PROBE_TESTED_SCALE2, multiplier=$MAG_FACTOR, candidate scales=$CANDIDATE_TESTED_SCALE,$CANDIDATE_TESTED_SCALE2, $ROUND_OBJECTIVE=$CANDIDATE_SCORE (current=$REFINE_SCORE)."

                if float_less "$CANDIDATE_SCORE" "$REFINE_SCORE"; then
                    mv "$REFINED_CANDIDATE" "$REFINE_GRAPH"
                    REFINE_SCORE="$CANDIDATE_SCORE"
                    REFINE_SCALE="$CANDIDATE_TESTED_SCALE"
                    REFINE_SCALE2="$CANDIDATE_TESTED_SCALE2"
                    echo "Pair refinement $REFINE_ITER accepted."
                else
                    rm -f "$REFINED_CANDIDATE"
                    echo "Pair refinement $REFINE_ITER rejected; Skipping refinement."
					break
                fi
            done

            cp "$REFINE_GRAPH" "$OUTPUT_ABS"
            echo "Final refined pair score: $ROUND_OBJECTIVE=$REFINE_SCORE with latest accepted scales $REFINE_SCALE,$REFINE_SCALE2."
        fi

        echo "Achieved graphlet frequency:"
        run_bnb "$GRAPHLET_K" "$OUTPUT_ABS"
    )
}

generate_outinp() {
        local InS="$1"
        local ROOT_DIR
        ROOT_DIR=$(pwd -P)

        local INS_ABS TARGET_ABS VEC_INS_ABS
        INS_ABS=$(realpath "$InS")
        TARGET_ABS=$(realpath "$TARGET")
        VEC_INS_ABS=$(realpath vecInS.txt)


        echo
        echo "Re-getting transformation results sequentially"
        echo

        local VAL

        for VAL in "${ACTIVE_CANDS[@]}"; do
                (
                        set -e

                        local WORK_DIR
                        WORK_DIR=$(mktemp -d "$ROOT_DIR/.generate_outinp_${VAL}.XXXXXX")
                        trap 'rm -rf -- "$WORK_DIR"' EXIT

                        prepare_serial_worker "$WORK_DIR" "$ROOT_DIR"

                        cd "$WORK_DIR"

                        # Bash uses dynamic scoping, so apply_exp sees this absolute TARGET.
                        local TARGET="$TARGET_ABS"
                        local BASE_SYNTH1="baseline1_val_${VAL}.txt"
                        local BASE_SYNTH2="baseline2_val_${VAL}.txt"
                        local OUT_INP="outinp_val_${VAL}.txt"
                        local SCALE TESTED_SCALE1 TESTED_SCALE2

                        echo "Getting $VAL result"
#set -x
                        SCALE=30000
                        apply_exp "$INS_ABS" "$SCALE" "$VAL" "$BASE_SYNTH1" \
                            TESTED_SCALE1
                        run_bnb "$GRAPHLET_K" "$BASE_SYNTH1" > "vecBS1_val_${VAL}.txt"
                        displ "$VEC_INS_ABS" "vecBS1_val_${VAL}.txt" "trfBS1_val_${VAL}.txt"

                        SCALE=300000
                        apply_exp "$INS_ABS" "$SCALE" "$VAL" "$BASE_SYNTH2" \
                            TESTED_SCALE2
                        run_bnb "$GRAPHLET_K" "$BASE_SYNTH2" > "vecBS2_val_${VAL}.txt"
#set +x
                        displ "$VEC_INS_ABS" "vecBS2_val_${VAL}.txt" "trfBS2_val_${VAL}.txt"

                        # Each transformation file is self-contained:
                        #   line 1: effective tested scales ka kb
                        #   line 2: displacement measured at ka
                        #   line 3: displacement measured at kb
                        printf '%s %s\n' "$TESTED_SCALE1" "$TESTED_SCALE2" > "$OUT_INP"
                        cat "trfBS1_val_${VAL}.txt" "trfBS2_val_${VAL}.txt" >> "$OUT_INP"

                        # Each transformation probe contains:
                        #   ka kb + D values at ka + D values at kb.
                        # Catch malformed probe data here instead of letting bxk4f
                        # fail later with only "incomplete input".
                        local OUT_INP_WORDS
                        OUT_INP_WORDS=$(wc -w < "$OUT_INP")
                        if (( OUT_INP_WORDS != OUTINP_WORDS_EXPECTED )); then
                                echo "Malformed $OUT_INP for transformation $VAL: expected $OUTINP_WORDS_EXPECTED values for k=$GRAPHLET_K, got $OUT_INP_WORDS." >&2
                                exit 1
                        fi

                        # mv is atomic on the same filesystem, so run_round never sees a half-written file.
                        mv "$OUT_INP" "$ROOT_DIR/outinp_val_${VAL}.txt"
                        echo "Finished $VAL result (tested scales: $TESTED_SCALE1, $TESTED_SCALE2)"
                )
        done
}

# Main loop
CURRENT="$INIT_SYNTH"
run_bnb "$GRAPHLET_K" "$TARGET" > vecTar.txt

# Phase 1: the requested number of dynamic-weight RSEAS rounds. Probe curves
# are regenerated on rounds 1, 1+RBGO, 1+2*RBGO, ... within this phase.
for ((i=1; i<=ROUNDS; i++)); do
    set_round_constraints "$i"
    run_bnb "$GRAPHLET_K" "$CURRENT" > vecInS.txt
    set_round_rseas_state "$i" vecInS.txt vecTar.txt

    if (( (i-1) % RBGO == 0 )); then
        refresh_rare_graphlets data_middle1.txt
        generate_outinp "$CURRENT"
    fi

    NEXT="tmp_round${i}.txt"
    run_round "$CURRENT" "$NEXT" RSEAS
    CURRENT="$NEXT"
done

# Phase 2: RMSE starts a new probe-refresh cadence. Thus its first round always
# regenerates outinp, followed by every RBGO-th round within the RMSE phase.
echo
if (( RMSE_ROUNDS > 0 )); then
    echo "Completed $ROUNDS RSEAS rounds. Starting $RMSE_ROUNDS equal-weight RMSE polishing rounds."
else
    echo "Completed $ROUNDS RSEAS rounds. No RMSE polishing rounds were requested."
fi

evb="$TIGHTEST_EVB"
degb="$TIGHTEST_DEGB"
hdgb="$TIGHTEST_HDGB"

for ((j=1; j<=RMSE_ROUNDS; j++)); do
    set_rmse_round_constraints "$j"

    if (( (j-1) % RBGO == 0 )); then
        run_bnb "$GRAPHLET_K" "$CURRENT" > vecInS.txt
        refresh_rare_graphlets data_middle1.txt
        generate_outinp "$CURRENT"
    fi

    TOTAL_INDEX=$((ROUNDS + j))
    NEXT="tmp_round${TOTAL_INDEX}.txt"
    run_round "$CURRENT" "$NEXT" RMSE
    CURRENT="$NEXT"
done

cp "$CURRENT" "$FINAL"
echo
echo "Final output stored in $FINAL after $ROUNDS RSEAS rounds and $RMSE_ROUNDS RMSE rounds; outinp was regenerated every $RBGO round(s) within each phase."
rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
          inp1l_*.txt synth_base_*.txt tmp_round*.txt PIS*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt bres.tmp evh.txt \
          out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt bextxzx.txt suminpxzx.txt data_middle1.txt idk.txt \
          expfy_samples_done.tmp .expfy_samples_*.tmp
rm -rf baseline .generate_outinp_* .run_round.*
