#pragma once

#include "src/trio/model/AlignmentSiteGraph.h"

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

class GfaWriter {
public:
    static void write(const AlignmentSiteGraph &graph, std::ostream &output);
    static void writeMetadata(const AlignmentSiteGraph &graph, std::ostream &output);
};

struct MafExportOmission {
    std::string componentId;
    std::string reason;
    std::vector<std::string> siteIds;
};

struct MafExportResult {
    std::size_t blocksWritten = 0;
    std::vector<MafExportOmission> omissions;
};

// MAF is a lossy projection. Only components with one deterministic linear
// site order and explicit observed-base/aligned-absence states are emitted.
// Runs split whenever path coverage changes; missing evidence is never '-'.
class MafGraphExporter {
public:
    static MafExportResult write(const AlignmentSiteGraph &graph,
                                 std::ostream &output);
    static void writeOmissions(const MafExportResult &result,
                               const std::string &coreHash,
                               std::ostream &output);
};

}  // namespace trio
}  // namespace anchorwave
