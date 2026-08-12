#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FULL_ROOT="${FULL_ROOT:-${REPO_ROOT}/../AnchorWave-wfa2-benchmark-data/installation}"
WINDOW_BASES="${WINDOW_BASES:-4000000}"
BENCHMARK_ROOT="${BENCHMARK_ROOT:-${FULL_ROOT}/bounded-chr1-${WINDOW_BASES}}"
WORK_DIR="${BENCHMARK_ROOT}/work"
RESULT_DIR="${BENCHMARK_ROOT}/results"

OLD_BIN="${OLD_BIN:-/private/tmp/anchorwave-wfa-old-build/anchorwave}"
NEW_BIN="${NEW_BIN:-/private/tmp/anchorwave-wfa-new-build/anchorwave}"
MINIMAP2_BIN="${MINIMAP2_BIN:-$(command -v minimap2 || true)}"
MINIMAP2_THREADS="${MINIMAP2_THREADS:-12}"
ANCHORWAVE_THREADS="${ANCHORWAVE_THREADS:-1}"
WINDOW_WIDTH="${WINDOW_WIDTH:-100000}"
MAX_MEMORY_GIB="${MAX_MEMORY_GIB:-}"

FULL_WORK="${FULL_ROOT}/work"
FULL_GFF="${FULL_WORK}/Zea_mays.AGPv4.34.gff3"
FULL_REFERENCE="${FULL_WORK}/Zea_mays.AGPv4.dna.toplevel.fa"
FULL_QUERY="${FULL_WORK}/Zm-Mo17-REFERENCE-CAU-1.0.fa"

GFF="${WORK_DIR}/chr1.${WINDOW_BASES}.gff3"
REFERENCE="${WORK_DIR}/B73.chr1.${WINDOW_BASES}.fa"
QUERY="${WORK_DIR}/Mo17.chr1.${WINDOW_BASES}.fa"
CDS="${WORK_DIR}/cds.fa"
CDS_SAM="${WORK_DIR}/cds.sam"
REF_SAM="${WORK_DIR}/ref.sam"

mkdir -p "${WORK_DIR}" "${RESULT_DIR}"

if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -dimsu -w "$$" &
fi

require_inputs() {
    local input
    for input in "${FULL_GFF}" "${FULL_REFERENCE}" "${FULL_QUERY}"; do
        if [[ ! -s "${input}" ]]; then
            printf 'required full installation input is missing: %s\n' \
                "${input}" >&2
            return 1
        fi
    done
    for input in "${OLD_BIN}" "${NEW_BIN}" "${MINIMAP2_BIN}"; do
        if [[ -z "${input}" || ! -x "${input}" ]]; then
            printf 'required executable is missing: %s\n' "${input}" >&2
            return 1
        fi
    done
}

extract_chr1_prefix() {
    local source_fasta="$1"
    local target_fasta="$2"
    if [[ -s "${target_fasta}" ]]; then
        return 0
    fi
    awk -v limit="${WINDOW_BASES}" '
        BEGIN { print ">1" }
        NR == 1 { next }
        /^>/ { exit }
        {
            remaining = limit - emitted
            if (remaining <= 0) exit
            if (length($0) <= remaining) {
                print $0
                emitted += length($0)
            } else {
                print substr($0, 1, remaining)
                emitted += remaining
                exit
            }
        }
        END {
            if (emitted != limit) {
                print "failed to extract requested chr1 prefix" > "/dev/stderr"
                exit 2
            }
        }
    ' "${source_fasta}" > "${target_fasta}.tmp"
    mv "${target_fasta}.tmp" "${target_fasta}"
}

extract_gff() {
    if [[ -s "${GFF}" ]]; then
        return 0
    fi
    awk -F '\t' -v limit="${WINDOW_BASES}" '
        /^#/ { print; next }
        $1 == "1" {
            seen = 1
            if ($4 >= 1 && $5 <= limit) print
            next
        }
        seen { exit }
    ' "${FULL_GFF}" > "${GFF}.tmp"
    mv "${GFF}.tmp" "${GFF}"
}

prepare() {
    require_inputs
    extract_chr1_prefix "${FULL_REFERENCE}" "${REFERENCE}"
    extract_chr1_prefix "${FULL_QUERY}" "${QUERY}"
    extract_gff
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
        printf 'GNU time is required\n' >&2
        return 1
    fi
}

run_version() {
    local label="$1"
    local binary="$2"
    local prefix="${RESULT_DIR}/${label}"
    if [[ -s "${prefix}.complete" ]]; then
        return 0
    fi
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
    printf 'completed\n' > "${prefix}.complete"
}

compare() {
    shasum -a 256 "${RESULT_DIR}"/{old,new}.anchorwave.maf \
        "${RESULT_DIR}"/{old,new}.anchorwave.fragment.maf \
        "${RESULT_DIR}"/{old,new}.anchors \
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
    {
        printf '[old]\n'
        awk -F '\t' 'NF >= 4 {count[$4]++} END {for (m in count) print m, count[m]}' \
            "${RESULT_DIR}/old.methods.bed" | sort
        printf '[new]\n'
        awk -F '\t' 'NF >= 4 {count[$4]++} END {for (m in count) print m, count[m]}' \
            "${RESULT_DIR}/new.methods.bed" | sort
    } > "${RESULT_DIR}/method-counts.txt"
}

record_provenance() {
    {
        printf 'scope=chr1:1-%s\n' "${WINDOW_BASES}"
        printf 'full_input_root=%s\n' "${FULL_ROOT}"
        printf 'old_binary_sha256='; shasum -a 256 "${OLD_BIN}" | awk '{print $1}'
        printf 'new_binary_sha256='; shasum -a 256 "${NEW_BIN}" | awk '{print $1}'
        printf 'wfa2_upstream=main@2fc4a1afac0f624e7020ee5aadb7692b38157eaa\n'
        printf 'anchorwave_threads=%s\n' "${ANCHORWAVE_THREADS}"
        printf 'window_width=%s\n' "${WINDOW_WIDTH}"
        printf 'max_memory_gib=%s\n' "${MAX_MEMORY_GIB:-unlimited-legacy-mode}"
        printf 'minimap2_threads=%s\n' "${MINIMAP2_THREADS}"
    } > "${RESULT_DIR}/benchmark-provenance.txt"
}

case "${1:-all}" in
    prepare) prepare ;;
    run-old) run_version old "${OLD_BIN}" ;;
    run-new) run_version new "${NEW_BIN}" ;;
    compare) compare ;;
    all)
        prepare
        record_provenance
        run_version old "${OLD_BIN}"
        run_version new "${NEW_BIN}"
        compare
        ;;
    *)
        printf 'Usage: %s prepare|run-old|run-new|compare|all\n' "$0" >&2
        exit 1
        ;;
esac
