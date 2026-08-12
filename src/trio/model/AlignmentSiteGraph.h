#pragma once

#include "src/trio/model/AlignmentEvidence.h"

#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

enum class ConsistencyClass {
    CONSISTENT,
    LOCAL_POA_REPAIRED,
    SUPPORTED_BY_TWO_EDGES,
    MISSING_EVIDENCE,
    GAP_SHIFT,
    BOUNDARY_DISCORDANCE,
    RESIDUE_TRANSITIVITY,
    LOW_CONFIDENCE_PAIRWISE,
    NONHOMOLOGOUS_INSERTIONS,
    COPY_COLLISION,
    ORDER_ORIENTATION,
    UNANCHORED_OR_TOO_LARGE
};

std::string consistencyClassName(ConsistencyClass value);

struct EvidenceDisposition {
    std::string evidenceId;
    bool accepted = false;
    ConsistencyClass classification = ConsistencyClass::MISSING_EVIDENCE;
    std::string reason;
};

struct AlignmentSite {
    std::string siteId;
    std::string copyGroup;
    bool copyResolved = false;
    std::vector<ResidueObservation> observations;
    std::set<std::string> directPairSupport;
    std::vector<std::string> evidenceIds;
    std::set<std::string> alignedAbsentPaths;
    std::vector<std::string> absenceEvidenceIds;
    ConsistencyClass consistency = ConsistencyClass::MISSING_EVIDENCE;
};

struct PathSegment {
    std::string pathId;
    std::string occurrencePath;
    std::string taxon;
    std::string sequence;
    int64_t start0 = 0;
    int64_t end0 = 0;
    int64_t sourceSize = 0;
    std::vector<std::string> siteIds;
    std::string sequenceText;
};

struct GraphConflict {
    std::string conflictId;
    ConsistencyClass classification = ConsistencyClass::MISSING_EVIDENCE;
    std::vector<std::string> evidenceIds;
    std::vector<ResidueId> residues;
    bool poaEligible = false;
    std::string disposition;
};

class AlignmentSiteGraph {
public:
    const std::map<std::string, AlignmentSite> &sites() const { return sites_; }
    const std::vector<PathSegment> &pathSegments() const { return pathSegments_; }
    const std::vector<EvidenceDisposition> &evidenceAudit() const { return evidenceAudit_; }
    const std::vector<GraphConflict> &conflicts() const { return conflicts_; }
    const std::string &ingroup1() const { return ingroup1_; }
    const std::string &ingroup2() const { return ingroup2_; }
    const std::string &primaryOutgroup() const { return primaryOutgroup_; }

    void validateInvariants() const;
    std::string coreHash() const;

private:
    friend class AlignmentSiteGraphBuilder;
    friend class GraphRepairEngine;
    std::map<std::string, AlignmentSite> sites_;
    std::vector<PathSegment> pathSegments_;
    std::vector<EvidenceDisposition> evidenceAudit_;
    std::vector<GraphConflict> conflicts_;
    std::string ingroup1_;
    std::string ingroup2_;
    std::string primaryOutgroup_;
};

class AlignmentGraphError : public std::runtime_error {
public:
    explicit AlignmentGraphError(const std::string &message)
        : std::runtime_error(message) {}
};

struct AlignmentGraphBuildOptions {
    std::string ingroup1;
    std::string ingroup2;
    std::string primaryOutgroup;
    std::map<ResidueId, CopyAssignment> copyAssignments;
    std::set<std::string> allowedEvidenceIds;
    bool restrictToAllowedEvidence = false;
    bool requireResolvedCopyAssignments = false;
    bool requirePrimaryTriangle = true;
};

class AlignmentSiteGraphBuilder {
public:
    static AlignmentSiteGraph build(const AlignmentEvidenceSet &evidence,
                                    const AlignmentGraphBuildOptions &options);
};

}  // namespace trio
}  // namespace anchorwave
