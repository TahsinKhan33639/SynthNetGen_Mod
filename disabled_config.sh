#!/usr/bin/env bash

# Shared by both pipelines. EVC exclusions are consulted only by run_round.
declare -A DISABLED_TRANSFORMATIONS=()
declare -A DISABLED_EVC=()

load_disabled_transformations() {
    local CONFIG_PATH="$1"
    local LINE TOKEN VALUE LIMIT
    local -a LINES=() IDS=()

    if [[ ! -f "$CONFIG_PATH" || ! -r "$CONFIG_PATH" ]]; then
        echo "Cannot read disabled transformation file: $CONFIG_PATH" >&2
        return 1
    fi
    mapfile -t LINES < "$CONFIG_PATH"
    if (( ${#LINES[@]} != 2 )); then
        echo "$CONFIG_PATH must contain exactly two lines." >&2
        return 1
    fi

    DISABLED_TRANSFORMATIONS=()
    DISABLED_EVC=()
    for LINE in 0 1; do
        LIMIT=95
        if (( LINE == 1 )); then LIMIT=3; fi
        IFS=$' \t\r\v\f' read -r -a IDS <<< "${LINES[$LINE]}"
        for TOKEN in "${IDS[@]}"; do
            if [[ ! "$TOKEN" =~ ^[0-9]+$ ]]; then
                echo "$CONFIG_PATH line $((LINE + 1)): invalid integer '$TOKEN'." >&2
                return 1
            fi
            # Strip leading zeroes before arithmetic to avoid octal parsing
            # and overflow from unusually long numeric tokens.
            VALUE="${TOKEN#"${TOKEN%%[!0]*}"}"
            VALUE="${VALUE:-0}"
            if (( LINE == 1 && ${#IDS[@]} == 1 )) && [[ "$VALUE" == 0 ]]; then
                continue
            fi
            if (( ${#VALUE} > 2 )) || (( VALUE < 1 || VALUE > LIMIT )); then
                echo "$CONFIG_PATH line $((LINE + 1)): IDs must be in [1, $LIMIT]; use 0 alone on line 2 for no EVC exclusions." >&2
                return 1
            fi
            if (( LINE == 0 )); then
                DISABLED_TRANSFORMATIONS[$VALUE]=1
            else
                DISABLED_EVC[$VALUE]=1
            fi
        done
    done

    # Retain the validated contents so each run can create its own snapshot.
    printf -v DISABLED_CONFIG_CONTENT '%s\n%s\n' "${LINES[0]}" "${LINES[1]}"
}
