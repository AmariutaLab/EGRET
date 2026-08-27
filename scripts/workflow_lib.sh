#!/bin/bash
# --- workflow_lib.sh ------------------------------------------------------
# Shared helpers sourced by every EGRET shell script: --flag parsing, conda
# activation, SLURM submission with retries, and precondition checks.
#
#   in : nothing -- sourced, never executed
#   out: shell functions only
#   run: sourced by all scripts/*.sh
#
# Every caller opens by locating this file, using a three-step fallback:
#
#   scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
#   [[ -r .../workflow_lib.sh ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
#
# BASH_SOURCE alone is not enough because SLURM copies a batch script into a private spool
# directory before running it, so inside a job BASH_SOURCE points at a directory that holds
# no workflow_lib.sh. The fallbacks are the path the submitting process exported, then the
# directory sbatch was invoked from. Each script re-derives this rather than calling a
# helper, because the helper lives in the file being located.
# --------------------------------------------------------------------------

cli_die() {
    echo "ERROR: $*" >&2
    exit 2
}

cli_usage() {
    printf '%b\n' "$USAGE"
}

# Callers define CLI_FLAGS, CLI_VARS, and CLI_REQUIRED arrays, initialize defaults,
# then call parse_cli "$@". Every option takes exactly one value.
parse_cli() {
    local arg value i found var
    while (($#)); do
        arg=$1
        shift
        if [[ "$arg" == "--help" ]]; then
            cli_usage
            exit 0
        fi
        [[ "$arg" == --* ]] || cli_die "positional arguments are not supported: $arg"
        found=0
        for i in "${!CLI_FLAGS[@]}"; do
            if [[ "$arg" == "${CLI_FLAGS[$i]}" ]]; then
                found=1
                (($#)) || cli_die "missing value for $arg"
                value=$1
                shift
                [[ "$value" != --* ]] || cli_die "missing value for $arg"
                var=${CLI_VARS[$i]}
                printf -v "$var" '%s' "$value"
                break
            fi
        done
        ((found)) || cli_die "unknown option: $arg"
    done

    for var in "${CLI_REQUIRED[@]}"; do
        [[ -n "${!var:-}" ]] || cli_die "missing required option: --${var}"
    done
}

activate_workflow_env() {
    local env_name=$1
    [[ -n "$env_name" ]] || cli_die "missing required option: --conda_env"
    # `conda activate` is a shell FUNCTION, not the binary on PATH, so `command -v conda`
    # is not a sufficient test: it succeeds on the binary while the function is undefined.
    # That is the state of every non-interactive shell -- every SLURM job here, and any
    # `./script.sh` run. Source the hook whenever the function is missing.
    if ! declare -f conda >/dev/null 2>&1; then
        local conda_base
        conda_base=$(conda info --base 2>/dev/null) || conda_base=""
        if [[ -n "$conda_base" && -r "${conda_base}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "${conda_base}/etc/profile.d/conda.sh"
        elif [[ -r "${HOME}/.bashrc" ]]; then
            # shellcheck disable=SC1090
            source "${HOME}/.bashrc"
        fi
    fi
    declare -f conda >/dev/null 2>&1 || {
        echo "ERROR: conda is unavailable; cannot activate $env_name" >&2
        exit 1
    }
    conda activate "$env_name"
}

join_by_colon() {
    local IFS=:
    echo "$*"
}

write_job_ids() {
    local path=$1
    shift
    mkdir -p "$(dirname "$path")"
    join_by_colon "$@" > "$path"
}

read_job_ids() {
    local path=$1
    [[ -s "$path" ]] || {
        echo "ERROR: current-run job ID file is missing or empty: $path" >&2
        return 1
    }
    tr -d '[:space:]' < "$path"
}

require_file() {
    local path=$1 label=${2:-file}
    [[ -s "$path" ]] || {
        echo "ERROR: expected $label is missing or empty: $path" >&2
        return 1
    }
}

require_dir_with_match() {
    local dir=$1 pattern=$2 label=${3:-output}
    [[ -d "$dir" ]] && compgen -G "$dir/$pattern" >/dev/null || {
        echo "ERROR: expected $label not found in $dir (pattern: $pattern)" >&2
        return 1
    }
}


# Count files under DIR matching the given find arguments, printing 0 when DIR does not
# exist yet. Callers run under `set -euo pipefail`, where the obvious spelling
#     n=$(find "$dir" -type f 2>/dev/null | wc -l)
# is a trap: on a missing directory find exits 1, 2>/dev/null hides the message, pipefail
# carries the failure through wc, and set -e kills the script with exit 1 and NO output --
# a job that dies in one second with an empty log is this bug.
count_files() {
    local dir=$1
    shift
    [[ -d "$dir" ]] || { echo 0; return 0; }
    local n
    n=$(find "$dir" "$@" 2>/dev/null | wc -l) || n=0
    echo "$n"
}

# sbatch with retries. A busy scheduler intermittently rejects a submission outright even
# with nodes free ("Required partition not available"); under set -e one such rejection
# ends the driver mid-dispatch and orphans every job downstream of it.
submit_with_retry() {
    local attempt=1 max=3 delay=10 out=""
    while true; do
        if out=$(sbatch "$@" 2>&1); then
            printf '%s\n' "$out"
            return 0
        fi
        if (( attempt >= max )); then
            echo "ERROR: sbatch failed after ${max} attempts: $out" >&2
            return 1
        fi
        echo "WARNING: sbatch attempt ${attempt}/${max} failed (${out}); retrying in ${delay}s" >&2
        sleep "$delay"
        attempt=$((attempt + 1)); delay=$((delay * 3))
    done
}
