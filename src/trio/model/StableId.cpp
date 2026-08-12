#include "src/trio/model/StableId.h"

#include <iomanip>
#include <sstream>

namespace anchorwave {
namespace trio {

uint64_t fnv1a64(const std::string &value) {
    uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char byte : value) {
        hash ^= static_cast<uint64_t>(byte);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

std::string stableId(const std::string &prefix,
                     const std::vector<std::string> &fields) {
    std::ostringstream framed;
    framed << "TrioAnchorGraph-v1";
    for (const std::string &field : fields) {
        framed << '|' << field.size() << ':' << field;
    }
    std::ostringstream result;
    result << prefix << '_' << std::hex << std::setw(16) << std::setfill('0')
           << fnv1a64(framed.str());
    return result.str();
}

}  // namespace trio
}  // namespace anchorwave

