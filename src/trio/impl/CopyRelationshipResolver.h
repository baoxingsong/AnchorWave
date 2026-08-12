#pragma once

#include "src/trio/model/AlignmentEvidence.h"
#include "src/trio/model/TrioTypes.h"

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

enum class CopyResolutionStatus {
    RESOLVED,
    ALTERNATIVES,
    UNRESOLVED
};

struct CopyResolutionAudit {
    std::string evidenceId;
    bool selected = false;
    std::string classification;
    std::string reason;
};

struct ResidueCopyAlternative {
    ResidueId residue;
    std::vector<std::string> copyGroups;
    std::vector<double> scores;
};

struct CopyResolutionResult {
    CopyResolutionStatus status = CopyResolutionStatus::RESOLVED;
    std::map<ResidueId, CopyAssignment> assignments;
    std::set<std::string> selectedEvidenceIds;
    std::vector<CopyResolutionAudit> audit;
    std::vector<ResidueCopyAlternative> alternatives;
};

struct CopyResolverOptions {
    bool strictExplicitGroups = false;
    double softTieMargin = 0.05;
};

class CopyResolutionError : public std::runtime_error {
public:
    explicit CopyResolutionError(const std::string &message)
        : std::runtime_error(message) {}
};

// Resolves copy compatibility before multi-sequence graph construction. Hard
// interval membership/exclusion always wins. Anchor-member constraints require
// an anchor-to-interval expansion layer; until that layer exists, resolve()
// rejects them explicitly instead of silently ignoring them. Soft constraints
// rank candidates; unambiguous pairwise components may receive an explicitly
// labelled inferred group. Count capacities are validation constraints only and
// never choose a paralog/homeolog identity.
class CopyRelationshipResolver {
public:
    static CopyResolutionResult resolve(const AlignmentEvidenceSet &evidence,
                                        const CopyConstraintSet &constraints,
                                        const CopyResolverOptions &options);
};

}  // namespace trio
}  // namespace anchorwave
