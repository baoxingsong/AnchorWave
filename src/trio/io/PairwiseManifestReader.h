#pragma once

#include "../model/TrioConfig.h"

#include <string>

namespace anchorwave {
namespace trio {

class PairwiseManifestReader {
public:
    explicit PairwiseManifestReader(
        const PairwiseManifestReaderOptions &options = PairwiseManifestReaderOptions());

    PairwiseManifest read(const std::string &path,
                          const TaxonManifest &taxa) const;

    static void validateRequiredPairs(const PairwiseManifest &manifest,
                                      const TaxonManifest &taxa,
                                      PairwiseScope scope,
                                      const std::string &manifestPath);

private:
    PairwiseManifestReaderOptions options_;
};

}  // namespace trio
}  // namespace anchorwave
