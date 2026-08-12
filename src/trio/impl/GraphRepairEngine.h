#pragma once

#include "src/trio/impl/LocalizedPoaRepair.h"
#include "src/trio/model/AlignmentSiteGraph.h"

#include <cstddef>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

struct GraphRepairOptions {
    LocalizedPoaOptions poa;
    std::size_t maximumCoreSites = 10000;
};

struct GraphRepairAudit {
    std::string conflictId;
    std::string repairId;
    std::string disposition;
    std::string immutableLeftSite;
    std::string immutableRightSite;
    std::vector<std::string> oldCoreSites;
    std::vector<std::string> newCoreSites;
    std::string outsideHashBefore;
    std::string outsideHashAfter;
    double baselineScore = 0.0;
    double repairedScore = 0.0;
    double scoreDelta = 0.0;
};

struct GraphRepairResult {
    AlignmentSiteGraph graph;
    std::vector<GraphRepairAudit> audit;
};

// Finds the smallest POA-eligible core bracketed by concordant sites shared by
// all implicated occurrence paths, adds any other path that traverses both
// boundaries, aligns only the core, validates exact source spelling, and then
// applies the candidate atomically to a graph copy.
class GraphRepairEngine {
public:
    static GraphRepairResult repairEligibleConflicts(
        const AlignmentSiteGraph &input,
        const GraphRepairOptions &options = GraphRepairOptions());
};

}  // namespace trio
}  // namespace anchorwave

