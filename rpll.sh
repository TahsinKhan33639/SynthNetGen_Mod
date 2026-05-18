# Usage: ./rpll.sh target output

rm -f baseline_val_*.txt outinp_val_*.txt synth1_*.txt \
      outinp_l_val_*.txt suminp_val_*.txt bext_val_*.txt \
	  inp1l_*.txt synth_base_*.txt tmp_round*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt \
	  out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt
rm -rf baseline

set -e

if [ $# -eq 2 ]; then
	echo "a"
	TARGET0="$1"
	FINAL="$2"
	TARGET="target_c.txt"

	./cleanup "$TARGET0" "$TARGET"
	./gen_deg2 "$TARGET" synth_base_1.txt
	./expfc "$TARGET" synth_base_1.txt 50000 50000 1 8 2.0 50 > synth_base_2.txt
	./expfc "$TARGET" synth_base_2.txt 10000 10000 1 8 1.8 40 > synth_base_3.txt
	./expfc "$TARGET" synth_base_3.txt 10000 10000 1 8 1.6 30 > synth_base_4.txt
	./expfc "$TARGET" synth_base_4.txt 10000 10000 1 8 1.5 24 > synth_base_5.txt
	./expfc "$TARGET" synth_base_5.txt 10000 10000 1 8 1.4 20 > synth_base_6.txt
	./expfc "$TARGET" synth_base_6.txt 10000 10000 1 8 1.3 15 > synth_base_7.txt
	./expfc "$TARGET" synth_base_7.txt 10000 10000 1 8 1.2 12 > synth_base_8.txt
	./expfc "$TARGET" synth_base_8.txt 10000 10000 1 8 1.1 10 > synth_base_9.txt
	./expfc "$TARGET" synth_base_9.txt 10000 10000 1 8 1.0 9 > synth_base_10.txt
	./expfc "$TARGET" synth_base_10.txt 10000 10000 1 8 0.9 8 > synth_base_11.txt
	./expfc "$TARGET" synth_base_11.txt 10000 10000 1 8 0.8 7 > synth_base_12.txt
	./expfc "$TARGET" synth_base_12.txt 10000 10000 1 8 0.7 6 > synth_base_13.txt
	./expfc "$TARGET" synth_base_13.txt 10000 10000 1 8 0.6 5 > synth_base_14.txt
	./expfc "$TARGET" synth_base_14.txt 10000 10000 1 8 0.5 4 > synth_base_15.txt
	./expfc "$TARGET" synth_base_15.txt 10000 10000 1 8 0.5 3 > synth_base_16.txt
	./expfc "$TARGET" synth_base_16.txt 10000 10000 1 8 0.5 2 > synth_base_17.txt

	INIT_SYNTH=synth_base_17.txt
elif [ $# -eq 3 ]; then
	TARGET0="$1"
	FINAL="$2"
	INIT_SYNTH="$3"
	TARGET="target_c.txt"
	./cleanup "$TARGET0" "$TARGET"
else
    echo "Usage: $0 <Target> <Output> or $0 <Target> <Output> <Base_Synth>"
    exit 1
fi

ROUNDS=24
evb=1
degb=4
CANDS=($(seq 1 56 | grep -v -E '^(25|26)$'))

declare -a ec
for i in {1..56}; do
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

# Apply expfc with one type
apply_exp() {
    local INPUT_FILE="$1"
    local SCALE="$2"
    local VAL="$3"
    local OUTPUT_FILE="$4"

	if [ "${ec[$VAL]}" -eq 0 ]; then
	    SCALE=$((SCALE / 60))
	else
	    SCALE=$((SCALE / 600))
	fi

	if [ "$SCALE" -lt 1 ]; then
	    SCALE=0
	fi

    if [ "$VAL" -ne 25 ] && [ "$VAL" -ne 26 ]; then
        ./expfc "$TARGET" "$INPUT_FILE" "$SCALE" 0 "$VAL" 1 "$evb" "$degb" > "$OUTPUT_FILE"
    fi
}


# Apply expfc with two types
apply_exp2() {
    local INPUT_FILE="$1"
    local SCALE="$2"
    local SCALE2="$3"
    local VAL="$4"
    local VAL2="$5"
    local OUTPUT_FILE="$6"

	if [ "${ec[$VAL2]}" -eq 0 ]; then
	    SCALE2=$((SCALE2 / 60))
	else
	    SCALE2=$((SCALE2 / 600))
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
        ./expfc "$TARGET" "$INPUT_FILE" "$SCALE" "$SCALE2" "$VAL" "$VAL2" "$evb" "$degb" > "$OUTPUT_FILE"
    fi
}


# optimization round
run_round() {
    local INPUT_SYNTH="$1"
    local OUTPUT_SYNTH="$2"

	displ vecInS.txt vecTar.txt inpl.txt
	cat vecInS.txt vecTar.txt > inp.txt

    echo "New round with input $INPUT_SYNTH"

    BEST_VAL2=""
	BEST_MAG2=""
    BEST_VAL=""
    BEST_RMSE=""
	BEST_MAG=""

	RESULTS_FILE="round_results.tmp"
	> "$RESULTS_FILE"

	# --- collect results ---
	for VAL in "${CANDS[@]}"; do
		for VAL2 in "${CANDS[@]}"; do
			if (( VAL2 > VAL )); then
            	continue
        	fi
		    SUMINP="suminp_val_${VAL}_${VAL2}.txt"
		    BEXT="bext_val_${VAL}_${VAL2}.txt"
		    OUT_INP="outinp_val_${VAL}.txt"
		    OUT_INP2="outinp_val_${VAL2}.txt"

		    cat inp.txt "$OUT_INP" "$OUT_INP2" > "$SUMINP"
		    ./bxk4ue "${ec[$VAL]}" "${ec[$VAL2]}" "$CMODE" < "$SUMINP" > "$BEXT"

		    read mag mag2 val < "$BEXT"

		    echo "$val|$VAL|$VAL2|$mag|$mag2" >> "$RESULTS_FILE"
		done
	done

	had_rmse=$(rmse f0.txt inpl.txt)
	echo "RMSE before modification is $had_rmse"
	FIRST_TRY=1

	# Compute top 10 pairs
	FIRST_LINE=$(sort -n "$RESULTS_FILE" | head -n 1)

	while IFS="|" read sorted_rmse sorted_VAL sorted_VAL2 sorted_MAG sorted_MAG2; do

	    SCALE=$(echo "($sorted_MAG * 1)/1" | bc)
	    SCALE2=$(echo "($sorted_MAG2 * 1)/1" | bc)

		if [ "$SCALE2" -le 0 ] && [ "$SCALE" -le 0 ]; then
		    if [ "$FIRST_TRY" -eq 1 ]; then
		        echo "First candidate has non-positive scales. Forcing evaluation."
		    else
		        continue
		    fi
		fi
		FIRST_TRY=0
		apply_exp2 "$INPUT_SYNTH" "$SCALE" "$SCALE2" "$sorted_VAL" "$sorted_VAL2" "$OUTPUT_SYNTH"

		./bnb 4 "$OUTPUT_SYNTH" > vecOutty.txt
	    displ vecInS.txt vecOutty.txt out_inplxzx.txt
	    got_rmse=$(rmse out_inplxzx.txt inpl.txt)

	    echo "VAL=$sorted_VAL,$sorted_VAL2 MAG=$sorted_MAG,MAG2=$sorted_MAG2. expected RMSE=$sorted_rmse, got RMSE=$got_rmse."

		if [ -z "$BEST_RMSE" ] || [ "$(echo "$got_rmse < $BEST_RMSE" | bc -l)" -eq 1 ]; then
		    BEST_RMSE="$got_rmse"
		    BEST_VAL="$sorted_VAL"
		    BEST_VAL2="$sorted_VAL2"
		    BEST_MAG="$sorted_MAG"
		    BEST_MAG2="$sorted_MAG2"
		fi

	done < <(sort -n "$RESULTS_FILE" | head -n 10)

    # Apply expfc to make OUTPUT_SYNTH

    SCALE=$(echo "($BEST_MAG * 1)/1" | bc)
    SCALE2=$(echo "($BEST_MAG2 * 1)/1" | bc)
	if [ "$(echo "$BEST_RMSE > $had_rmse" | bc -l)" -eq 1 ]; then
		SCALE=0
		SCALE2=0
	fi
	echo "Best pair was $BEST_VAL and $BEST_VAL2 with rmse $BEST_RMSE"
	apply_exp2 "$INPUT_SYNTH" "$SCALE" "$SCALE2" "$BEST_VAL" "$BEST_VAL2" "$OUTPUT_SYNTH"
	echo "Achieved graphlet frequency:"
	./bnb 4 "$OUTPUT_SYNTH"
}

generate_outinp() {
	local InS="$1"
	echo
    echo "Re-getting transformation results"
	echo

	for VAL in "${CANDS[@]}"; do
		echo "Getting $VAL result"
		BASE_SYNTH1="baseline1_val_${VAL}.txt"
		BASE_SYNTH2="baseline2_val_${VAL}.txt"
		BASE_SYNTH3="baseline3_val_${VAL}.txt"
		BASE_SYNTH4="baseline4_val_${VAL}.txt"
	    OUT_INP="outinp_val_${VAL}.txt"
		> "$OUT_INP"

		SCALE=1500
		apply_exp "$InS" "$SCALE" "$VAL" "$BASE_SYNTH1"
		./bnb 4 "$BASE_SYNTH1" > vecBS1.txt
		displ vecInS.txt vecBS1.txt trfBS.txt
		cat trfBS.txt >> "$OUT_INP"

		SCALE=8000
		apply_exp "$InS" "$SCALE" "$VAL" "$BASE_SYNTH2"
		./bnb 4 "$BASE_SYNTH2" > vecBS2.txt
		displ vecInS.txt vecBS2.txt trfBS.txt
		cat trfBS.txt >> "$OUT_INP"

		SCALE=50000
		apply_exp "$InS" "$SCALE" "$VAL" "$BASE_SYNTH3"
		./bnb 4 "$BASE_SYNTH3" > vecBS3.txt
		displ vecInS.txt vecBS3.txt trfBS.txt
		cat trfBS.txt >> "$OUT_INP"

		SCALE=300000
		apply_exp "$InS" "$SCALE" "$VAL" "$BASE_SYNTH4"
		./bnb 4 "$BASE_SYNTH4" > vecBS4.txt
		displ vecInS.txt vecBS4.txt trfBS.txt
		cat trfBS.txt >> "$OUT_INP"

		rm "$BASE_SYNTH1" "$BASE_SYNTH2" "$BASE_SYNTH3" "$BASE_SYNTH4"
	done
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
        evb=1
		degb=4
    fi
    if (( i == 7 )); then
        evb=1
		degb=3
    fi
    if (( i == 11 )); then
        evb=1
		degb=3
    fi
    if (( i == 15 )); then
        evb=0.7
		degb=3
    fi
    if (( i == 23 )); then
        evb=0.5
		degb=2
    fi

    if (( (i-1) % 4 == 0 )); then
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
	  inp1l_*.txt synth_base_*.txt tmp_round*.txt
rm -f round_results.tmp inp.txt inpl.txt vecBS*.txt vecOutty.txt trfBS.txt target_c.txt smth.txt \
	  out_inplxzx.txt data1.txt data_middle.txt vecInS.txt vecTar.txt
rm -rf baseline
