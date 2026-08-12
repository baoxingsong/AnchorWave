#pragma once

#include "TrioTypes.h"

#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

struct ManifestReaderOptions {
    bool validatePaths;

    ManifestReaderOptions() : validatePaths(false) {}
};

struct PairwiseManifestReaderOptions : public ManifestReaderOptions {
    PairwiseScope scope;

    PairwiseManifestReaderOptions() : scope(PairwiseScope::TRIANGLES) {}
};

struct TrioConfig {
    std::string taxaManifestPath;
    std::string speciesTreePath;
    std::string pairwiseManifestPath;
    std::string copyRelationsPath;
    std::string outputPrefix;
    PairwiseScope pairwiseScope;
    std::vector<std::string> defaultCopyCapacitySpecs;
    bool validateInputPaths;

    TrioConfig()
        : pairwiseScope(PairwiseScope::TRIANGLES), validateInputPaths(false) {}
};

}  // namespace trio
}  // namespace anchorwave
