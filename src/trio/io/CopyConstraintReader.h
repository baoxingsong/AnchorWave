#pragma once

#include "../model/TrioConfig.h"

#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

class CopyConstraintReader {
public:
    explicit CopyConstraintReader(
        const ManifestReaderOptions &options = ManifestReaderOptions());

    CopyConstraintSet read(
        const std::string &path,
        const TaxonManifest &taxa,
        const std::vector<std::string> &defaultCapacitySpecs =
            std::vector<std::string>()) const;

    static CopyCapacity parseDefaultCapacitySpec(
        const std::string &specification,
        const SourceLocation &location = SourceLocation("<copy-number>", 0));

private:
    ManifestReaderOptions options_;
};

}  // namespace trio
}  // namespace anchorwave
