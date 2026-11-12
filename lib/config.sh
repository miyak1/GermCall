#!/usr/bin/env bash
# =============================================================================
#  GermCall: shared helpers (config reader).  Usage at the top of the main script:
#
#      SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#      source "${SCRIPT_DIR}/lib/config.sh"
#      init_config "$@"
#
#  After init_config you have:
#      CONFIG   absolute path of the YAML file
#      WORKDIR  project.workdir (absolute); the shell is already cd'd into it
#      PREFIX   project.prefix
#      THREADS  project.threads
#      EXTRA_ARGS  array with command-line arguments that were not --config
#  and the functions cfg, cfg_list, path_of, dir_of, tool, require_tool,
#  require_file, banner.
# =============================================================================

GC_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GC_ROOT="$(cd "${GC_LIB_DIR}/.." && pwd)"

# cfg KEY [DEFAULT]  -> scalar value (lists: one per line)
cfg() { python3 "${GC_LIB_DIR}/cfg.py" "${CONFIG}" "$1" "${2-}"; }

# cfg_list KEY -> one item per line (empty if missing)
cfg_list() { python3 "${GC_LIB_DIR}/cfg.py" "${CONFIG}" "$1" ""; }

# path_of VALUE -> absolute path (VALUE relative to WORKDIR unless absolute / ~)
path_of() {
    local v="$1"
    [[ -z "$v" ]] && { echo ""; return; }
    v="${v/#\~/$HOME}"
    if [[ "$v" = /* ]]; then echo "$v"; else echo "${WORKDIR}/${v}"; fi
}

# dir_of KEY DEFAULT -> absolute path of a dirs.* style setting
dir_of() { path_of "$(cfg "$1" "$2")"; }

# tool NAME DEFAULT -> command configured as tools.NAME
tool() { cfg "tools.$1" "$2"; }

# require_tool COMMAND [HINT]
require_tool() {
    if ! command -v "$1" >/dev/null 2>&1 && [[ ! -x "$1" ]]; then
        echo "ERROR: required tool not found: $1" >&2
        [[ -n "${2-}" ]] && echo "       $2" >&2
        echo "       Set its location under 'tools:' in ${CONFIG}" >&2
        exit 127
    fi
}

# require_file FILE [HINT]
require_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: input file not found: $1" >&2
        [[ -n "${2-}" ]] && echo "       $2" >&2
        exit 1
    fi
}

banner() {
    echo "============================================================"
    echo "  $*"
    echo "============================================================"
}

init_config() {
    local cfg_arg=""
    EXTRA_ARGS=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --config|-c) cfg_arg="$2"; shift 2 ;;
            --config=*)  cfg_arg="${1#--config=}"; shift ;;
            *)           EXTRA_ARGS+=("$1"); shift ;;
        esac
    done

    CONFIG="${cfg_arg:-${GERMCALL_CONFIG:-${GC_ROOT}/config/config.yaml}}"
    if [[ ! -f "${CONFIG}" ]]; then
        echo "ERROR: config file not found: ${CONFIG}" >&2
        echo "       Copy config/config.example.yaml, edit it," >&2
        echo "       then pass it with  --config path/to/file.yaml  (or set GERMCALL_CONFIG)." >&2
        exit 1
    fi
    CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"
    export GERMCALL_CONFIG="${CONFIG}"

    command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found" >&2; exit 127; }

    local wd
    wd="$(cfg project.workdir .)"
    wd="${wd/#\~/$HOME}"
    if [[ "$wd" != /* ]]; then
        # a relative workdir is relative to the directory you launched from
        wd="$(pwd)/${wd}"
    fi
    mkdir -p "${wd}"
    WORKDIR="$(cd "${wd}" && pwd)"
    cd "${WORKDIR}"

    PREFIX="$(cfg project.prefix '')"
    if [[ -z "${PREFIX}" ]]; then
        echo "ERROR: project.prefix is empty in ${CONFIG}" >&2
        exit 1
    fi
    THREADS="$(cfg project.threads 4)"
}
