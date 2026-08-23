#!/bin/bash

rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
	  inp1l_*.txt synth_base_*.txt tmp_round*.txt PIS*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt bres.tmp \
	  out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt bextxzx.txt suminpxzx.txt data_middle1.txt idk.txt
rm -rf baseline .generate_outinp_* .run_round.*

set -e

if [ $# -eq 3 ]; then
	TARGET0="$1"
	FINAL="$2"
	ROUNDS="$3"
	TARGET="target_c.txt"

	./cleanup "$TARGET0" "$TARGET"
	./gen_deg2 "$TARGET" synth_base_1.txt

	./bnb 4 synth_base_1.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_1.txt 50000 50000 1 8 2.0 20 50 > synth_base_2.txt
	./bnb 4 synth_base_2.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_2.txt 10000 10000 1 8 1.7 18 40 > synth_base_3.txt
	./bnb 4 synth_base_3.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_3.txt 10000 10000 1 8 1.5 18 30 > synth_base_4.txt
	./bnb 4 synth_base_4.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_4.txt 10000 10000 1 8 1.3 17 26 > synth_base_5.txt
	./bnb 4 synth_base_5.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_5.txt 10000 1000 1 17 1.2 16 24 > synth_base_6.txt
	./bnb 4 synth_base_6.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_6.txt 10000 1000 1 17 1.1 15 23 > synth_base_7.txt
	./bnb 4 synth_base_7.txt > bres.tmp
	cat data_middle.txt > data_middle1.txt
	./expfy "$TARGET" synth_base_7.txt 10000 1000 1 17 1 14 22 > synth_base_8.txt

	PRE_INIT_SYNTH=synth_base_8.txt
elif [ $# -eq 4 ]; then
	TARGET0="$1"
	FINAL="$2"
	PRE_INIT_SYNTH="$3"
	ROUNDS="$4"
	TARGET="target_c.txt"
	./cleanup "$TARGET0" "$TARGET"
else
    echo "Usage: $0 <Target> <Output> <Rounds> | $0 <Target> <Output> <Base_Synth> <Rounds>"
    exit 1
fi

./evs "$PRE_INIT_SYNTH" > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
if [ "$(echo "$EV1 < $EV2 - 14" | bc -l)" -eq 1 ]; then
	./evc1 "$PRE_INIT_SYNTH" 3000 > PIS1.txt
elif [ "$(echo "$EV1 > $EV2 + 14" | bc -l)" -eq 1 ]; then
	./evc2 "$PRE_INIT_SYNTH" 3000 > PIS1.txt
else
	cp "$PRE_INIT_SYNTH" PIS1.txt
fi

./evs PIS1.txt > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
if [ "$(echo "$EV1 < $EV2 - 7" | bc -l)" -eq 1 ]; then
	./evc1 PIS1.txt 1000 > PIS2.txt
elif [ "$(echo "$EV1 > $EV2 + 7" | bc -l)" -eq 1 ]; then
	./evc2 PIS1.txt 1000 > PIS2.txt
else
    cp PIS1.txt PIS2.txt
fi

./evs PIS2.txt > evh.txt
read EV1 < evh.txt
./evs "$TARGET" > evh.txt
read EV2 < evh.txt
if [ "$(echo "$EV1 < $EV2 - 2" | bc -l)" -eq 1 ]; then
	./evc1 PIS2.txt 400 > PIS3.txt
elif [ "$(echo "$EV1 > $EV2 + 2" | bc -l)" -eq 1 ]; then
	./evc2 PIS2.txt 400 > PIS3.txt
else
	cp PIS2.txt PIS3.txt
fi

INIT_SYNTH=PIS3.txt

evb=0.7
degb=5
hdgb=10
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

# Log displacement
displ() {
    python3 - "$1" "$2" "$3" << 'EOF'
import sys, math

A = list(map(float, open(sys.argv[1]).read().split()))
B = list(map(float, open(sys.argv[2]).read().split()))

diff = []
for a, b in zip(A, B):
    if a <= 0 or b <= 0:
        diff.append(0.0)
    else:
        diff.append(math.log(a) - math.log(b))
open(sys.argv[3], "w").write(" ".join(map(str, diff)) + "\n")
EOF
}

# Linear displacement
disp() {
    python3 - "$1" "$2" "$3" << 'EOF'
import sys, math

A = list(map(float, open(sys.argv[1]).read().split()))
B = list(map(float, open(sys.argv[2]).read().split()))

diff = []
for a, b in zip(A, B):
    if a <= 0 or b <= 0:
        diff.append(0.0)
    else:
        diff.append(a - b)
open(sys.argv[3], "w").write(" ".join(map(str, diff)) + "\n")
EOF
}


CMODE=0

mapfile -t weights_table < weights4.txt

# RMSE
rmse() {
    read -a W <<< "${weights_table[$CMODE]}"

    python3 - "$1" "$2" "${W[@]}" << 'EOF'
import sys, math

fileA, fileB, cmode = sys.argv[1], sys.argv[2], int(sys.argv[3])

A = list(map(float, open(fileA).read().split()))
B = list(map(float, open(fileB).read().split()))

w = list(map(float, sys.argv[3:]))

n = min(len(A), len(B), 6)

num = sum(w[i] * (A[i]-B[i])**2 for i in range(n))
rmse = math.sqrt(num)
print(rmse)
EOF
}

# Create an isolated current directory for a parallel worker.  Several of the
# project executables create fixed-name scratch files such as data_middle.txt,
# so merely changing the shell output filename is not enough to prevent races.
prepare_parallel_worker() {
    local WORK_DIR="$1"
    local ROOT_DIR="$2"
    local SHARED_FILE TOOL

    mkdir -p "$WORK_DIR"

    for SHARED_FILE in data_middle.txt data_middle1.txt weights4.txt f0.txt; do
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

# Apply expfy with one type
apply_exp() {
    local INPUT_FILE="$1"
    local SCALE="$2"
    local VAL="$3"
    local OUTPUT_FILE="$4"

	if [ "${ec[$VAL]}" -eq 0 ]; then
	    SCALE=$((SCALE / 600))
	else
	    SCALE=$((SCALE / 6000))
	fi

	if [ "$SCALE" -lt 1 ]; then
	    SCALE=0
	fi

    if [ "$VAL" -ne 25 ] && [ "$VAL" -ne 26 ]; then
        ./expfy "$TARGET" "$INPUT_FILE" "$SCALE" 0 "$VAL" 1 "$evb" "$degb" "$hdgb" > "$OUTPUT_FILE"
    fi
}


# Apply expfy with two types
apply_exp2() {
    local INPUT_FILE="$1"
    local SCALE="$2"
    local SCALE2="$3"
    local VAL="$4"
    local VAL2="$5"
    local OUTPUT_FILE="$6"

	if [ "${ec[$VAL2]}" -eq 0 ]; then
	    SCALE2=$((SCALE2 / 600))
	else
	    SCALE2=$((SCALE2 / 6000))
	fi

	if [ "$SCALE2" -lt 1 ]; then
	    SCALE2=0
	fi

	if [ "${ec[$VAL]}" -eq 0 ]; then
	    SCALE=$((SCALE / 60))
	else
	    SCALE=$((SCALE / 600))
	fi

	if [ "$SCALE" -lt 1 ]; then
	    SCALE=0
	fi

    if [ "$VAL" -ne 25 ] && [ "$VAL" -ne 26 ]; then
        ./expfy "$TARGET" "$INPUT_FILE" "$SCALE" "$SCALE2" "$VAL" "$VAL2" "$evb" "$degb" "$hdgb" > "$OUTPUT_FILE"
    fi
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
        BXK4F_ABS=$(realpath ./bxk4f)
        BXK4ONE_ABS=$(realpath ./bxk4one)

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

        local had_rmse
        had_rmse=$(rmse "$F0_ABS" "$INPL")

        echo "New round with input $INPUT_SYNTH"
        echo "Running independent round work with up to $MAX_JOBS parallel jobs"

        # ---------------------------------------------------------
        # 1. Generate and evaluate the four EVC candidates in parallel.
        # ---------------------------------------------------------
        local EVC_DIR="$ROUND_DIR/evc"
        mkdir -p "$EVC_DIR"

        local -a EVC_LABELS=(outpev11 outpev12 outpev21 outpev22)
        local -a EVC_TOOLS=(evc1 evc1 evc2 evc2)
        local -a EVC_AMOUNTS=(100 200 100 200)
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
                SCORE=$(rmse candidate.displ "$INPL")

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
        local BEST_EVC_RMSE="$had_rmse"
        local CUR_RMSE CUR_GRAPH

        # Read in the original deterministic order so equal scores are handled
        # exactly as they were in the serial version.
        for LABEL in "${EVC_LABELS[@]}"; do
            IFS='|' read -r CUR_RMSE CUR_GRAPH < "$EVC_DIR/$LABEL.result"
            echo "$LABEL RMSE = $CUR_RMSE"

            if float_less "$CUR_RMSE" "$BEST_EVC_RMSE"; then
                BEST_EVC_RMSE="$CUR_RMSE"
                BEST_EVC_GRAPH="$CUR_GRAPH"
            fi
        done

        if [ "$BEST_EVC_GRAPH" != "$INPUT_ABS" ]; then
            echo "Replacing $INPUT_SYNTH with the best EVC candidate (RMSE $BEST_EVC_RMSE)"
            cp "$BEST_EVC_GRAPH" "$INPUT_ABS"
            had_rmse="$BEST_EVC_RMSE"
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

                        cat "$INP" "$OUT_INP" "$OUT_INP2" > pair_input.tmp
                        "$BXK4F_ABS" "${ec[$VAL]}" "${ec[$VAL2]}" "$CMODE" \
                            < pair_input.tmp > pair_output.tmp

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

        had_rmse=$(rmse "$F0_ABS" "$INPL")
        echo "RMSE before modification is $had_rmse"

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
        local sorted_rmse sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2
        local SCALE SCALE2

        for IDX in "${!TOP_LINES[@]}"; do
            RANK=$((IDX + 1))
            LINE="${TOP_LINES[$IDX]}"
            IFS='|' read -r sorted_rmse sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2 <<< "$LINE"

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
                local MAG_FACTOR GOT_RMSE
                local ADJUSTED_SCALE="$SCALE"
                local ADJUSTED_SCALE2="$SCALE2"

                apply_exp2 "$INPUT_ABS" "$ADJUSTED_SCALE" "$ADJUSTED_SCALE2" \
                    "$sorted_VAL" "$sorted_VAL2" "$STAGE1"

                ./bnb 4 "$STAGE1" > stage1.vec
                displ "$VEC_INS" stage1.vec stage1.displ
                cat "$INP" stage1.displ > scale_input.tmp
                "$BXK4ONE_ABS" "$CMODE" < scale_input.tmp > scale_output.tmp

                read -r MAG_FACTOR < scale_output.tmp
                ADJUSTED_SCALE=$(awk -v S="$ADJUSTED_SCALE" -v M="$MAG_FACTOR" \
                    'BEGIN { print int(S * M) }')
                ADJUSTED_SCALE2=$(awk -v S="$ADJUSTED_SCALE2" -v M="$MAG_FACTOR" \
                    'BEGIN { print int(S * M) }')

                apply_exp2 "$INPUT_ABS" "$ADJUSTED_SCALE" "$ADJUSTED_SCALE2" \
                    "$sorted_VAL" "$sorted_VAL2" "$FINAL_CANDIDATE"

                ./bnb 4 "$FINAL_CANDIDATE" > final.vec
                displ "$VEC_INS" final.vec final.displ
                GOT_RMSE=$(rmse final.displ "$INPL")

                printf '%s|%s|%s|%s|%s|%s|%s\n' \
                    "$GOT_RMSE" "$sorted_VAL" "$sorted_VAL2" \
                    "$ADJUSTED_SCALE" "$ADJUSTED_SCALE2" "$sorted_rmse" "$RANK" \
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
        local BEST_RMSE=""
        local GOT_RMSE EXPECTED_RMSE RESULT_FILE

        # Consume results by original rank, preserving the serial tie-breaking.
        for ((RANK=1; RANK<=${#TOP_LINES[@]}; RANK++)); do
            RESULT_FILE="$TOP_DIR/result_$(printf '%02d' "$RANK").txt"
            if [ ! -f "$RESULT_FILE" ]; then
                continue
            fi

            IFS='|' read -r GOT_RMSE sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2 \
                EXPECTED_RMSE _ < "$RESULT_FILE"

            echo "VAL=$sorted_VAL,$sorted_VAL2 MAG=$sorted_MAG,MAG2=$sorted_MAG2. expected RMSE=$EXPECTED_RMSE, got RMSE=$GOT_RMSE."

            if [ -z "$BEST_RMSE" ] || float_less "$GOT_RMSE" "$BEST_RMSE"; then
                BEST_RMSE="$GOT_RMSE"
                BEST_VAL="$sorted_VAL"
                BEST_VAL2="$sorted_VAL2"
                BEST_MAG="$sorted_MAG"
                BEST_MAG2="$sorted_MAG2"
            fi
        done

        if [ -z "$BEST_RMSE" ]; then
            echo "No top-pair candidate completed successfully." >&2
            exit 1
        fi

        SCALE=$(integer_part "$BEST_MAG")
        SCALE2=$(integer_part "$BEST_MAG2")
        if float_greater "$BEST_RMSE" "$had_rmse"; then
            SCALE=0
            SCALE2=0
        fi

        echo "Best pair was $BEST_VAL and $BEST_VAL2 with rmse $BEST_RMSE"
        apply_exp2 "$INPUT_ABS" "$SCALE" "$SCALE2" \
            "$BEST_VAL" "$BEST_VAL2" "$OUTPUT_ABS"

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
			for SHARED_FILE in data_middle.txt data_middle1.txt weights4.txt; do
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
			local SCALE

			echo "Getting $VAL result"
			: > "$OUT_INP"

			SCALE=30000
			apply_exp "$INS_ABS" "$SCALE" "$VAL" "$BASE_SYNTH1"
			./bnb 4 "$BASE_SYNTH1" > "vecBS1_val_${VAL}.txt"
			displ "$VEC_INS_ABS" "vecBS1_val_${VAL}.txt" "trfBS1_val_${VAL}.txt"
			cat "trfBS1_val_${VAL}.txt" >> "$OUT_INP"

			SCALE=300000
			apply_exp "$INS_ABS" "$SCALE" "$VAL" "$BASE_SYNTH2"
			./bnb 4 "$BASE_SYNTH2" > "vecBS2_val_${VAL}.txt"
			displ "$VEC_INS_ABS" "vecBS2_val_${VAL}.txt" "trfBS2_val_${VAL}.txt"
			cat "trfBS2_val_${VAL}.txt" >> "$OUT_INP"

			# mv is atomic on the same filesystem, so run_round never sees a half-written file.
			mv "$OUT_INP" "$ROOT_DIR/outinp_val_${VAL}.txt"
			echo "Finished $VAL result"
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

for ((i=1; i<=ROUNDS; i++)); do
	if (( i == 5 )); then CMODE=1; fi
	if (( i == 9 )); then CMODE=2; fi
	if (( i == 13 )); then CMODE=3; fi
	if (( i == 17 )); then CMODE=4; fi
	if (( i == 21 )); then CMODE=5; fi
	./bnb 4 "$CURRENT" > vecInS.txt

    if (( i == 3 )); then
        evb=0.6
		degb=3
		hdgb=7
    fi
    if (( i == 7 )); then
        evb=1
		degb=17
		hdgb=85
    fi
    if (( i == 11 )); then
        evb=0.9
		degb=16
		hdgb=80
    fi
    if (( i == 15 )); then
        evb=0.7
		degb=15
		hdgb=75
    fi
    if (( i == 23 )); then
        evb=0.5
		degb=14
		hdgb=70
    fi

    if (( (i-1) % 4 == 0 )); then
		./bnb 4 "$CURRENT" > bres.tmp
		cat data_middle.txt > data_middle1.txt
		refresh_rare_graphlets data_middle1.txt
        generate_outinp "$CURRENT"
    fi

    NEXT="tmp_round${i}.txt"
    run_round "$CURRENT" "$NEXT"
    CURRENT="$NEXT"
done

cp "$CURRENT" "$FINAL"
echo
echo "Final output stored in $FINAL"
rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
	  inp1l_*.txt synth_base_*.txt tmp_round*.txt PIS*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt bres.tmp \
	  out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt bextxzx.txt suminpxzx.txt data_middle1.txt idk.txt
rm -rf baseline .generate_outinp_* .run_round.
