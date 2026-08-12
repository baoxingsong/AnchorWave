#pragma once

#include "src/trio/io/MafEvidenceReader.h"

#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

struct ResidueId {
    std::string taxon;
    std::string occurrencePath;
    std::string sequence;
    int64_t forwardPosition0 = 0;

    bool operator<(const ResidueId &other) const;
    bool operator==(const ResidueId &other) const;
    std::string canonicalString() const;
};

struct ResidueObservation {
    ResidueId id;
    char base = 'N';
    int64_t sourceSize = 0;
};

struct CopyAssignment {
    std::string copyGroup;
    bool hard = false;
    double confidence = 0.0;
    std::string provenance;
};

struct EvidenceProvenance {
    std::string evidenceId;
    std::string pairId;
    std::string runId;
    std::string mafPath;
    std::size_t blockIndex = 0;
    std::size_t column = 0;
    double alignmentScore = 0.0;
};

struct HomologyEvidence {
    ResidueId left;
    ResidueId right;
    EvidenceProvenance provenance;
};

// A pairwise gap is attached to the present residue and to the nearest observed
// forward-strand coordinates bracketing the gap on the other path. A value of
// -1 means that a flank is outside the MAF block and is therefore unknown, not
// biologically absent.
struct AlignedAbsenceEvidence {
    ResidueId present;
    std::string absentTaxon;
    std::string absentOccurrencePath;
    std::string absentSequence;
    int64_t absentLeftFlank0 = -1;
    int64_t absentRightFlank0 = -1;
    EvidenceProvenance provenance;
};

struct PairwiseMafInput {
    std::string leftTaxon;
    std::string rightTaxon;
    std::string mafPath;
    std::string runId;
    double sourceWeight = 1.0;
};

struct AlignmentEvidenceSet {
    std::map<ResidueId, ResidueObservation> residues;
    // Declared MAF source sizes are tracked independently of residue overlap.
    // This catches inconsistent pairwise runs even when their aligned blocks
    // cover disjoint coordinates of the same source sequence.
    std::map<std::pair<std::string, std::string>, int64_t> sourceSizes;
    std::vector<HomologyEvidence> homologies;
    std::vector<AlignedAbsenceEvidence> alignedAbsences;
    std::set<std::string> observedPairs;
};

class AlignmentEvidenceError : public std::runtime_error {
public:
    explicit AlignmentEvidenceError(const std::string &message)
        : std::runtime_error(message) {}
};

class PairwiseEvidenceNormalizer {
public:
    static AlignmentEvidenceSet normalize(const std::vector<PairwiseMafInput> &inputs);

    static void appendBlocks(const PairwiseMafInput &input,
                             const std::vector<MafBlock> &blocks,
                             AlignmentEvidenceSet &destination);

    static std::string canonicalPairId(const std::string &taxonA,
                                       const std::string &taxonB);
};

}  // namespace trio
}  // namespace anchorwave
