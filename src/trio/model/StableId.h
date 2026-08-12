#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

uint64_t fnv1a64(const std::string &value);

// Creates a deterministic, insertion-order-independent identifier. Fields are
// length-delimited before hashing so ["ab", "c"] differs from ["a", "bc"].
std::string stableId(const std::string &prefix,
                     const std::vector<std::string> &fields);

}  // namespace trio
}  // namespace anchorwave

