#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TRACE_TOOL="${SCRIPT_DIR}/b73_mo17_trace.py"
DEFAULT_DATA_ROOT="${REPO_ROOT}/../AnchorWave-wfa2-benchmark-data/installation"
DEFAULT_BOUNDED_WORK="${DEFAULT_DATA_ROOT}/bounded-chr1-4000000/work"
DEFAULT_ANCHORWAVE="${REPO_ROOT}/build-dynamic-selector/anchorwave"

usage() {
    cat <<'USAGE'
Usage:
  benchmark_b73_mo17_trace.sh replay [b73_mo17_trace.py replay options]
  benchmark_b73_mo17_trace.sh compare [b73_mo17_trace.py compare options]
  benchmark_b73_mo17_trace.sh bounded --output DIR [options]

The replay command accepts -t/--threads and -M/--max-memory-gib, but the
public `anchorwave ali` command itself has neither option.  They define a
conservative number of independent ali processes; this is not the production
process-wide scheduler.  `bounded` runs the real 4-Mb genoAli pipeline and
passes -t/-M directly to AnchorWave.

bounded options:
  --output DIR             New immutable result directory (required)
  --work-dir DIR           Prepared 4-Mb work directory
  --anchorwave FILE        AnchorWave executable
  -t, --threads INT        AnchorWave threads (default: 12)
  -M, --max-memory-gib N   Process limit in GiB (default: 80)
  -w, --window-width INT   Window/memory parameter (default: 100000)
  -bt, --bt-minutes N      Exact-DP policy: 0 exact-first; positive balanced
                           time ceiling in minutes (default: 30)
  --prefix NAME            Output prefix inside DIR (default: run)

No command downloads data.  Existing result directories are never overwritten.
USAGE
}

find_gnu_time() {
    if command -v gtime >/dev/null 2>&1; then
        command -v gtime
    elif [[ -x /usr/bin/time ]] && /usr/bin/time --version >/dev/null 2>&1; then
        printf '%s\n' /usr/bin/time
    else
        printf 'GNU time is required (gtime on macOS)\n' >&2
        return 1
    fi
}

run_bounded() {
    local output=""
    local work_dir="${DEFAULT_BOUNDED_WORK}"
    local anchorwave="${DEFAULT_ANCHORWAVE}"
    local threads=12
    local max_memory_gib=80
    local window_width=100000
    local bt_minutes=30
    local prefix=run

    while (( $# > 0 )); do
        case "$1" in
            --output) output="$2"; shift 2 ;;
            --work-dir) work_dir="$2"; shift 2 ;;
            --anchorwave) anchorwave="$2"; shift 2 ;;
            -t|--threads) threads="$2"; shift 2 ;;
            -M|--max-memory-gib) max_memory_gib="$2"; shift 2 ;;
            -w|--window-width) window_width="$2"; shift 2 ;;
            -bt|--bt-minutes) bt_minutes="$2"; shift 2 ;;
            --prefix) prefix="$2"; shift 2 ;;
            -h|--help) usage; return 0 ;;
            *) printf 'unknown bounded option: %s\n' "$1" >&2; return 2 ;;
        esac
    done
    if [[ -z "${output}" ]]; then
        printf 'bounded requires --output DIR\n' >&2
        return 2
    fi
    if [[ -e "${output}" ]]; then
        printf 'refusing to overwrite existing result path: %s\n' "${output}" >&2
        return 2
    fi
    if [[ ! -x "${anchorwave}" ]]; then
        printf 'AnchorWave binary is not executable: %s\n' "${anchorwave}" >&2
        return 2
    fi

    local gff="${work_dir}/chr1.4000000.gff3"
    local reference="${work_dir}/B73.chr1.4000000.fa"
    local query="${work_dir}/Mo17.chr1.4000000.fa"
    local cds="${work_dir}/cds.fa"
    local cds_sam="${work_dir}/cds.sam"
    local ref_sam="${work_dir}/ref.sam"
    local input
    for input in "${gff}" "${reference}" "${query}" "${cds}" "${cds_sam}" "${ref_sam}"; do
        if [[ ! -s "${input}" ]]; then
            printf 'required bounded input is missing: %s\n' "${input}" >&2
            return 2
        fi
    done

    mkdir -p "${output}"
    local absolute_output
    absolute_output="$(cd "${output}" && pwd)"
    local result_prefix="${absolute_output}/${prefix}"
    local tmp_suffix=".tmp.$$"
    local gnu_time
    gnu_time="$(find_gnu_time)"
    local memory_args=()
    if [[ "${max_memory_gib}" != "0" && "${max_memory_gib}" != "0.0" ]]; then
        memory_args=(-M "${max_memory_gib}")
    fi

    python3 - "${absolute_output}/scenario.json" "${anchorwave}" "${threads}" \
        "${max_memory_gib}" "${window_width}" "${bt_minutes}" "${work_dir}" <<'PY'
import hashlib, json, pathlib, sys
out, binary, threads, memory, window, bt, work = sys.argv[1:]
h = hashlib.sha256()
with open(binary, "rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        h.update(chunk)
pathlib.Path(out).write_text(json.dumps({
    "format_version": 1,
    "kind": "b73_mo17_bounded_genoali",
    "anchorwave": str(pathlib.Path(binary).resolve()),
    "anchorwave_sha256": h.hexdigest(),
    "requested_threads": int(threads),
    "requested_max_memory_gib": float(memory),
    "window_width": int(window),
    "bt_minutes": float(bt),
    "bounded_work_dir": str(pathlib.Path(work).resolve()),
    "scope": "B73/Mo17 chromosome 1 bases 1-4000000",
}, indent=2, sort_keys=True) + "\n")
PY

    local command=(
        "${anchorwave}" genoAli
        -i "${gff}" -as "${cds}" -r "${reference}"
        -a "${cds_sam}" -ar "${ref_sam}" -s "${query}"
        -n "${result_prefix}.anchors${tmp_suffix}"
        -o "${result_prefix}.maf${tmp_suffix}"
        -f "${result_prefix}.fragment.maf${tmp_suffix}"
        -b "${result_prefix}.methods.bed${tmp_suffix}"
        --alignment-trace "${result_prefix}.alignment-trace.tsv${tmp_suffix}"
        -t "${threads}" -w "${window_width}" -bt "${bt_minutes}" -IV
        "${memory_args[@]}"
    )
    {
        printf 'command='
        printf '%q ' "${command[@]}"
        printf '\n'
    } > "${absolute_output}/command.txt"

    "${gnu_time}" -v -o "${result_prefix}.time.txt${tmp_suffix}" \
        "${command[@]}" > /dev/null 2> "${result_prefix}.stderr.txt${tmp_suffix}"

    mv "${result_prefix}.anchors${tmp_suffix}" "${result_prefix}.anchors"
    mv "${result_prefix}.maf${tmp_suffix}" "${result_prefix}.maf"
    mv "${result_prefix}.fragment.maf${tmp_suffix}" "${result_prefix}.fragment.maf"
    mv "${result_prefix}.methods.bed${tmp_suffix}" "${result_prefix}.methods.bed"
    mv "${result_prefix}.alignment-trace.tsv${tmp_suffix}" \
        "${result_prefix}.alignment-trace.tsv"
    mv "${result_prefix}.time.txt${tmp_suffix}" "${result_prefix}.time.txt"
    mv "${result_prefix}.stderr.txt${tmp_suffix}" "${result_prefix}.stderr.txt"
    python3 "${TRACE_TOOL}" summarize-pipeline \
        --result-dir "${absolute_output}" --prefix "${prefix}"
    printf 'completed\n' > "${absolute_output}/complete"
}

case "${1:-}" in
    replay)
        shift
        gnu_time="$(find_gnu_time)"
        python3 "${TRACE_TOOL}" replay --gnu-time "${gnu_time}" "$@"
        ;;
    compare)
        shift
        python3 "${TRACE_TOOL}" compare "$@"
        ;;
    bounded)
        shift
        run_bounded "$@"
        ;;
    -h|--help|help|"")
        usage
        ;;
    *)
        printf 'unknown command: %s\n' "$1" >&2
        usage >&2
        exit 2
        ;;
esac
