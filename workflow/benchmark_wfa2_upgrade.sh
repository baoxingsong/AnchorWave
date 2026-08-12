#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BENCHMARK_ROOT="${BENCHMARK_ROOT:-${REPO_ROOT}/../AnchorWave-wfa2-benchmark-data/installation}"
DOWNLOAD_DIR="${BENCHMARK_ROOT}/downloads"
WORK_DIR="${BENCHMARK_ROOT}/work"
RESULT_DIR="${BENCHMARK_ROOT}/results"

OLD_BIN="${OLD_BIN:-/private/tmp/anchorwave-wfa-old-build/anchorwave}"
NEW_BIN="${NEW_BIN:-/private/tmp/anchorwave-wfa-new-build/anchorwave}"
MINIMAP2_BIN="${MINIMAP2_BIN:-$(command -v minimap2 || true)}"
ANCHORWAVE_THREADS="${ANCHORWAVE_THREADS:-1}"
WINDOW_WIDTH="${WINDOW_WIDTH:-100000}"
MAX_MEMORY_GIB="${MAX_MEMORY_GIB:-}"

if command -v nproc >/dev/null 2>&1; then
    DEFAULT_MINIMAP2_THREADS="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then
    DEFAULT_MINIMAP2_THREADS="$(sysctl -n hw.logicalcpu 2>/dev/null || true)"
    DEFAULT_MINIMAP2_THREADS="${DEFAULT_MINIMAP2_THREADS:-1}"
else
    DEFAULT_MINIMAP2_THREADS=1
fi
MINIMAP2_THREADS="${MINIMAP2_THREADS:-${DEFAULT_MINIMAP2_THREADS}}"

GFF_GZ="${DOWNLOAD_DIR}/Zea_mays.AGPv4.34.gff3.gz"
REFERENCE_GZ="${DOWNLOAD_DIR}/Zea_mays.AGPv4.dna.toplevel.fa.gz"
QUERY_GZ="${DOWNLOAD_DIR}/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz"
GFF="${WORK_DIR}/Zea_mays.AGPv4.34.gff3"
REFERENCE="${WORK_DIR}/Zea_mays.AGPv4.dna.toplevel.fa"
QUERY="${WORK_DIR}/Zm-Mo17-REFERENCE-CAU-1.0.fa"
CDS="${WORK_DIR}/cds.fa"
CDS_SAM="${WORK_DIR}/cds.sam"
REF_SAM="${WORK_DIR}/ref.sam"

mkdir -p "${DOWNLOAD_DIR}" "${WORK_DIR}" "${RESULT_DIR}"

# Long single-thread maize alignments can run for many hours. On macOS, keep
# the host awake only for this script's lifetime; resource timing still wraps
# AnchorWave directly rather than the caffeinate helper.
if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -dimsu -w "$$" &
fi

usage() {
    printf '%s\n' \
        "Usage: $0 download|prepare|provenance|run-old|run-new|compare|all" \
        "" \
        "Environment:" \
        "  OLD_BIN, NEW_BIN             old/new AnchorWave executables" \
        "  BENCHMARK_ROOT               data/result directory" \
        "  MINIMAP2_BIN                 minimap2 executable" \
        "  MINIMAP2_THREADS             mapping threads (default: all logical CPUs)" \
        "  ANCHORWAVE_THREADS           genoAli threads (default: 1)" \
        "  WINDOW_WIDTH                 -w value (default: 100000)"
}

download_one() {
    local url="$1"
    local target="$2"
    local attempt=1
    if [[ -s "${target}" ]] && gzip -t "${target}" 2>/dev/null; then
        return 0
    fi
    while true; do
        if curl --fail --location --continue-at - --silent --show-error \
                --output "${target}" "${url}" && \
                gzip -t "${target}" 2>/dev/null; then
            return 0
        fi
        if (( attempt >= 20 )); then
            printf 'download failed after %d attempts: %s\n' \
                "${attempt}" "${url}" >&2
            return 1
        fi
        attempt=$((attempt + 1))
        sleep 5
    done
}

download() {
    # Retry in an outer loop. Each new curl process re-reads the current file
    # length for -C -, whereas curl's internal retry can reuse an older offset.
    download_one \
        https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-34/gff3/zea_mays/Zea_mays.AGPv4.34.gff3.gz \
        "${GFF_GZ}"
    download_one \
        https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-34/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz \
        "${REFERENCE_GZ}"
    download_one \
        https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz \
        "${QUERY_GZ}"
}

decompress_once() {
    local source_file="$1"
    local target_file="$2"
    if [[ ! -s "${target_file}" ]]; then
        gzip -dc "${source_file}" > "${target_file}.tmp"
        mv "${target_file}.tmp" "${target_file}"
    fi
}

prepare() {
    if [[ -z "${MINIMAP2_BIN}" || ! -x "${MINIMAP2_BIN}" ]]; then
        printf 'minimap2 is required; set MINIMAP2_BIN\n' >&2
        return 1
    fi
    if [[ ! -x "${NEW_BIN}" ]]; then
        printf 'NEW_BIN is not executable: %s\n' "${NEW_BIN}" >&2
        return 1
    fi

    decompress_once "${GFF_GZ}" "${GFF}"
    decompress_once "${REFERENCE_GZ}" "${REFERENCE}"
    decompress_once "${QUERY_GZ}" "${QUERY}"
    perl -pi -e 's/^>chr/>/' "${QUERY}"

    if [[ ! -s "${CDS}" ]]; then
        "${NEW_BIN}" gff2seq -r "${REFERENCE}" -i "${GFF}" -o "${CDS}.tmp"
        mv "${CDS}.tmp" "${CDS}"
    fi
    if [[ ! -s "${CDS_SAM}" ]]; then
        "${MINIMAP2_BIN}" -x splice -t "${MINIMAP2_THREADS}" -k 12 \
            -a -p 0.4 -N 20 "${QUERY}" "${CDS}" > "${CDS_SAM}.tmp"
        mv "${CDS_SAM}.tmp" "${CDS_SAM}"
    fi
    if [[ ! -s "${REF_SAM}" ]]; then
        "${MINIMAP2_BIN}" -x splice -t "${MINIMAP2_THREADS}" -k 12 \
            -a -p 0.4 -N 20 "${REFERENCE}" "${CDS}" > "${REF_SAM}.tmp"
        mv "${REF_SAM}.tmp" "${REF_SAM}"
    fi
}

time_command() {
    local time_output="$1"
    shift
    if command -v gtime >/dev/null 2>&1; then
        gtime -v -o "${time_output}" "$@"
    elif [[ -x /usr/bin/time ]] && /usr/bin/time --version >/dev/null 2>&1; then
        /usr/bin/time -v -o "${time_output}" "$@"
    else
        printf 'GNU time is required (brew install gnu-time on macOS)\n' >&2
        return 1
    fi
}

record_provenance() {
    local provenance="${RESULT_DIR}/benchmark-provenance.txt"
    {
        printf 'recorded_utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        printf 'repository=%s\n' "${REPO_ROOT}"
        printf 'repository_head=%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || printf unknown)"
        printf 'build_type=Release\n'
        printf 'old_wfa2_snapshot=repository vendored VERSION v2.3 before upgrade\n'
        printf 'wfa2_upstream=main@2fc4a1afac0f624e7020ee5aadb7692b38157eaa\n'
        printf 'old_binary=%s\n' "${OLD_BIN}"
        printf 'new_binary=%s\n' "${NEW_BIN}"
        if [[ -x "${OLD_BIN}" ]]; then
            shasum -a 256 "${OLD_BIN}" | sed 's/^/old_binary_sha256=/'
        fi
        if [[ -x "${NEW_BIN}" ]]; then
            shasum -a 256 "${NEW_BIN}" | sed 's/^/new_binary_sha256=/'
        fi
        printf 'anchorwave_threads=%s\n' "${ANCHORWAVE_THREADS}"
        printf 'window_width=%s\n' "${WINDOW_WIDTH}"
        printf 'max_memory_gib=%s\n' "${MAX_MEMORY_GIB:-unlimited-legacy-mode}"
        printf 'approximate_per_thread_budget_bytes=%s\n' \
            "$((WINDOW_WIDTH * WINDOW_WIDTH))"
        printf 'minimap2=%s\n' "${MINIMAP2_BIN}"
        printf 'minimap2_threads=%s\n' "${MINIMAP2_THREADS}"
        printf 'system='; uname -a
    } > "${provenance}"
}

run_version() {
    local label="$1"
    local binary="$2"
    local prefix="${RESULT_DIR}/${label}"
    local complete="${prefix}.complete"
    if [[ -s "${complete}" ]]; then
        printf '%s benchmark already complete\n' "${label}"
        return 0
    fi
    if [[ ! -x "${binary}" ]]; then
        printf '%s binary is not executable: %s\n' "${label}" "${binary}" >&2
        return 1
    fi
    record_provenance

    local memory_args=()
    if [[ -n "${MAX_MEMORY_GIB}" ]]; then
        memory_args=(-M "${MAX_MEMORY_GIB}")
    fi

    time_command "${prefix}.time.txt" \
        "${binary}" genoAli \
        -i "${GFF}" -as "${CDS}" -r "${REFERENCE}" \
        -a "${CDS_SAM}" -ar "${REF_SAM}" -s "${QUERY}" \
        -n "${prefix}.anchors.tmp" \
        -o "${prefix}.anchorwave.maf.tmp" \
        -f "${prefix}.anchorwave.fragment.maf.tmp" \
        -b "${prefix}.methods.bed.tmp" \
        -t "${ANCHORWAVE_THREADS}" -w "${WINDOW_WIDTH}" -IV \
        "${memory_args[@]}" \
        > /dev/null 2> "${prefix}.stderr.log"

    mv "${prefix}.anchors.tmp" "${prefix}.anchors"
    mv "${prefix}.anchorwave.maf.tmp" "${prefix}.anchorwave.maf"
    mv "${prefix}.anchorwave.fragment.maf.tmp" \
        "${prefix}.anchorwave.fragment.maf"
    mv "${prefix}.methods.bed.tmp" "${prefix}.methods.bed"
    printf 'completed\n' > "${complete}"
}

compare() {
    shasum -a 256 \
        "${RESULT_DIR}/old.anchorwave.maf" \
        "${RESULT_DIR}/new.anchorwave.maf" \
        "${RESULT_DIR}/old.anchorwave.fragment.maf" \
        "${RESULT_DIR}/new.anchorwave.fragment.maf" \
        "${RESULT_DIR}/old.anchors" \
        "${RESULT_DIR}/new.anchors" \
        "${RESULT_DIR}/old.methods.bed" \
        "${RESULT_DIR}/new.methods.bed" \
        > "${RESULT_DIR}/sha256.txt"
    {
        if cmp -s <(tail -n +2 "${RESULT_DIR}/old.anchors") \
                  <(tail -n +2 "${RESULT_DIR}/new.anchors"); then
            printf 'data_verdict=EXACT\n'
        else
            printf 'data_verdict=DIFFERENT\n'
        fi
        if cmp -s <(sed -n '1p' "${RESULT_DIR}/old.anchors") \
                  <(sed -n '1p' "${RESULT_DIR}/new.anchors"); then
            printf 'command_header_identical=yes\n'
        else
            printf 'command_header_identical=no (expected old/new output paths)\n'
        fi
        wc -l "${RESULT_DIR}/old.anchors" "${RESULT_DIR}/new.anchors"
    } > "${RESULT_DIR}/anchor-comparison.txt"
    {
        printf '[old]\n'
        awk -F '\t' 'NF >= 4 {count[$4]++} END {for (method in count) print method, count[method]}' \
            "${RESULT_DIR}/old.methods.bed" | sort
        printf '[new]\n'
        awk -F '\t' 'NF >= 4 {count[$4]++} END {for (method in count) print method, count[method]}' \
            "${RESULT_DIR}/new.methods.bed" | sort
    } > "${RESULT_DIR}/method-counts.txt"
    python3 -B "${SCRIPT_DIR}/compare_maf.py" \
        "${RESULT_DIR}/old.anchorwave.maf" \
        "${RESULT_DIR}/new.anchorwave.maf" \
        --json "${RESULT_DIR}/alignment-comparison.json" --report-only \
        --score-penalties 6,8,2,75,1
    python3 -B "${SCRIPT_DIR}/compare_maf.py" \
        "${RESULT_DIR}/old.anchorwave.fragment.maf" \
        "${RESULT_DIR}/new.anchorwave.fragment.maf" \
        --json "${RESULT_DIR}/fragment-comparison.json" --report-only \
        --score-penalties 6,8,2,75,1
}

case "${1:-}" in
    download) download ;;
    prepare) prepare ;;
    provenance) record_provenance ;;
    run-old) run_version old "${OLD_BIN}" ;;
    run-new) run_version new "${NEW_BIN}" ;;
    compare) compare ;;
    all)
        download
        prepare
        record_provenance
        run_version old "${OLD_BIN}"
        run_version new "${NEW_BIN}"
        compare
        ;;
    *) usage; exit 1 ;;
esac
