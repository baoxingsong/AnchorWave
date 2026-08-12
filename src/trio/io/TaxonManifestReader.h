#pragma once

#include "../model/TrioConfig.h"

#include <string>

namespace anchorwave {
namespace trio {

class TaxonManifestReader {
public:
    explicit TaxonManifestReader(
        const ManifestReaderOptions &options = ManifestReaderOptions());

    TaxonManifest read(const std::string &path) const;

private:
    ManifestReaderOptions options_;
};

}  // namespace trio
}  // namespace anchorwave
