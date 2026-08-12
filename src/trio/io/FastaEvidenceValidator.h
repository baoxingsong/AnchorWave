#pragma once

#include "src/trio/model/AlignmentEvidence.h"

#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

/**
 * One deterministic validation finding.
 *
 * position0 is -1 when a finding is not tied to a residue.  line is zero when
 * no FASTA line is applicable (for example, a missing FASTA mapping).  The
 * machine-readable code is stable; message is intended for humans.
 */
struct FastaEvidenceDiagnostic {
    std::string code;
    std::string taxon;
    std::string fastaPath;
    std::string sequence;
    int64_t position0 = -1;
    std::size_t line = 0;
    std::string message;

    bool operator<(const FastaEvidenceDiagnostic &other) const;
    std::string canonicalString() const;
};

/** Aggregate counters from a single-pass inspection. */
struct FastaEvidenceValidationCounts {
    std::uint64_t taxaRequired = 0;
    std::uint64_t fastaFilesScanned = 0;
    std::uint64_t sequencesScanned = 0;
    std::uint64_t basesScanned = 0;
    std::uint64_t sourcesRequired = 0;
    std::uint64_t sourcesFound = 0;
    std::uint64_t sourcesValidated = 0;
    std::uint64_t residueObservationsRequired = 0;
    std::uint64_t residueObservationsValidated = 0;
};

/**
 * Inspection result.  diagnostics is always sorted independently of FASTA
 * record order, which makes logs and tests reproducible.
 */
struct FastaEvidenceValidationReport {
    FastaEvidenceValidationCounts counts;
    std::vector<FastaEvidenceDiagnostic> diagnostics;

    bool valid() const { return diagnostics.empty(); }
};

class FastaEvidenceValidationError : public std::runtime_error {
public:
    explicit FastaEvidenceValidationError(
        const FastaEvidenceValidationReport &report);

    const FastaEvidenceValidationReport &report() const { return report_; }

private:
    FastaEvidenceValidationReport report_;
};

/**
 * Validate normalized pairwise evidence against the forward FASTA spelling.
 *
 * Requirements are derived from both AlignmentEvidenceSet::sourceSizes and
 * AlignmentEvidenceSet::residues.  Only taxa required by the evidence are
 * opened.  Each required taxon's FASTA is streamed exactly once and whole
 * chromosome strings are never retained.  zlib transparent reading supports
 * ordinary FASTA and gzip/BGZF-compatible sequential FASTA input.
 */
class FastaEvidenceValidator {
public:
    using TaxonFastaMap = std::map<std::string, std::string>;

    /** Return all validation findings without throwing for content errors. */
    static FastaEvidenceValidationReport inspect(
        const AlignmentEvidenceSet &evidence,
        const TaxonFastaMap &fastaByTaxon);

    /**
     * Return counters on success; throw FastaEvidenceValidationError carrying
     * the complete report when any evidence, FASTA, or compatibility check
     * fails.
     */
    static FastaEvidenceValidationCounts validate(
        const AlignmentEvidenceSet &evidence,
        const TaxonFastaMap &fastaByTaxon);
};

}  // namespace trio
}  // namespace anchorwave
