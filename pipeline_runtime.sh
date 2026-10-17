#!/usr/bin/env bash

# Resolve caller paths before moving each run into its own workspace.
initialize_pipeline_run() {
    local COUNTER="$1"
    local TOOL TMP_PARENT
    local -a REQUIRED_TOOLS=(cleanup evs evc1 evc2 evc3
        bxk4f_rseas bxk4f_rmse bxk4one_rseas bxk4one_rmse "$COUNTER" orca)

    TARGET0=$(realpath -e -- "$TARGET0") || return 1
    FINAL=$(realpath -m -- "$FINAL") || return 1
    if [[ ! -f "$TARGET0" || ! -r "$TARGET0" ]]; then
        echo "Target graph is not a readable file: $TARGET0" >&2
        return 1
    fi
    if [[ -d "$FINAL" || ! -d "$(dirname -- "$FINAL")" ]]; then
        echo "Output must name a file in an existing directory: $FINAL" >&2
        return 1
    fi
    if (( INITIAL_GRAPH_GIVEN )); then
        PRE_INIT_SYNTH=$(realpath -e -- "$PRE_INIT_SYNTH") || return 1
        if [[ ! -f "$PRE_INIT_SYNTH" || ! -r "$PRE_INIT_SYNTH" ]]; then
            echo "Base graph is not a readable file: $PRE_INIT_SYNTH" >&2
            return 1
        fi
    else
        REQUIRED_TOOLS+=(gen_deg2)
    fi
    if [[ "$COUNTER" == onb ]]; then
        REQUIRED_TOOLS+=(otb2)
    fi

    EXPFY_BIN=$(realpath -e -- "${EXPFY_BIN:-$SCMD_PROJECT_DIR/expfy}") || return 1
    if [[ ! -x "$EXPFY_BIN" ]]; then
        echo "Transformation backend is not executable: $EXPFY_BIN. Run ./build.sh first." >&2
        return 1
    fi
    for TOOL in "${REQUIRED_TOOLS[@]}"; do
        if [[ ! -x "$SCMD_PROJECT_DIR/$TOOL" ]]; then
            echo "Missing executable: $SCMD_PROJECT_DIR/$TOOL. Run ./build.sh first." >&2
            return 1
        fi
    done
    if [[ ! -r "$SCMD_PROJECT_DIR/f0.txt" ]]; then
        echo "Missing zero-reference vector: $SCMD_PROJECT_DIR/f0.txt" >&2
        return 1
    fi

    if [[ "$COUNTER" == bnb ]]; then
        BLANT_RESOURCE_DIR="$SCMD_PROJECT_DIR"
        if [[ ! -x "$BLANT_RESOURCE_DIR/blant" && -x "$SCMD_PROJECT_DIR/Synthcheck/blant" ]]; then
            BLANT_RESOURCE_DIR="$SCMD_PROJECT_DIR/Synthcheck"
        fi
        if [[ ! -x "$BLANT_RESOURCE_DIR/blant" ]]; then
            echo "BLANT is missing from the project resource directory." >&2
            return 1
        fi
    fi

    TMP_PARENT=$(realpath -e -- "${TMPDIR:-/tmp}") || return 1
    RUN_DIR=$(mktemp -d "$TMP_PARENT/scmd-run.XXXXXXXX") || return 1
    SCMD_RUN_OWNER=$BASHPID
    EXPFY_LOG=""
    trap cleanup_pipeline_run EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM

    mkdir -- "$RUN_DIR/expfy_logs" || return 1
    EXPFY_LOG=$(mktemp "${FINAL}.expfy.XXXXXXXX.log") || return 1
    echo "expfy stderr log: $EXPFY_LOG" >&2

    for TOOL in "${REQUIRED_TOOLS[@]}"; do
        ln -s -- "$SCMD_PROJECT_DIR/$TOOL" "$RUN_DIR/$TOOL"
    done
    cp -- "$SCMD_PROJECT_DIR/f0.txt" "$RUN_DIR/f0.txt"
    cp -- "$TARGET0" "$RUN_DIR/input_target.el"
    TARGET0="$RUN_DIR/input_target.el"
    if (( INITIAL_GRAPH_GIVEN )); then
        cp -- "$PRE_INIT_SYNTH" "$RUN_DIR/input_base.el"
        PRE_INIT_SYNTH="$RUN_DIR/input_base.el"
    fi
    cd -- "$RUN_DIR"

    if ! EXPFY_MODE_LIMIT=$(invoke_expfy --max-transformation); then
        echo "Cannot query transformation support from $EXPFY_BIN. Rebuild with ./build.sh." >&2
        return 1
    fi
    if ! [[ "$EXPFY_MODE_LIMIT" =~ ^[1-9][0-9]{0,8}$ ]]; then
        echo "Invalid transformation limit from $EXPFY_BIN: '$EXPFY_MODE_LIMIT'." >&2
        return 1
    fi
}

# Separate logs prevent concurrent workers from interleaving diagnostics.
invoke_expfy() {
    local CALL_LOG STATUS
    CALL_LOG=$(mktemp "$RUN_DIR/expfy_logs/call.XXXXXXXX.log") || return 1
    {
        printf 'Working directory: %s\nCommand:' "$PWD"
        printf ' %q' "$EXPFY_BIN" "$@"
        printf '\n'
    } > "$CALL_LOG"
    if "$EXPFY_BIN" "$@" 2>> "$CALL_LOG"; then
        STATUS=0
    else
        STATUS=$?
    fi
    printf '\nExit status: %s\n\n' "$STATUS" >> "$CALL_LOG"
    return "$STATUS"
}

cleanup_pipeline_run() {
    local STATUS=$?
    if [[ "$BASHPID" != "$SCMD_RUN_OWNER" ]]; then
        return "$STATUS"
    fi
    trap - EXIT
    local -a LOG_FILES=("$RUN_DIR"/expfy_logs/*.log)
    if [[ -n "$EXPFY_LOG" && -f "${LOG_FILES[0]}" ]]; then
        if ! cat -- "${LOG_FILES[@]}" > "$EXPFY_LOG"; then
            echo "Cannot save expfy diagnostics to $EXPFY_LOG; private files retained at $RUN_DIR." >&2
            if (( STATUS == 0 )); then STATUS=1; fi
            return "$STATUS"
        fi
    fi
    rm -rf -- "$RUN_DIR"
    return "$STATUS"
}

