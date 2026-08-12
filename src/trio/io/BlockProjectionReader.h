#pragma once

#include "src/trio/evolution/AncestralAdjacency.h"
#include "src/trio/model/TrioTypes.h"

#include <string>

namespace anchorwave {
namespace trio {

// Versioned external macro-synteny projection input. Coordinates are forward,
// zero-based, half-open extant-genome coordinates. This reader deliberately
// does not consume ancestral sequence calls.
class BlockProjectionReader {
public:
    static BlockProjectionSet read(const std::string &path,
                                   const TaxonManifest &taxa,
                                   const AlignmentSiteGraph &graph);
};

}  // namespace trio
}  // namespace anchorwave
