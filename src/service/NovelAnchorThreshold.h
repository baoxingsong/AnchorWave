#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace anchorwave {

inline uint64_t novelAnchorMinimumArea(int64_t minimumLength) {
    if (minimumLength <= 0) {
        throw std::invalid_argument("novel-anchor minimum length must be positive");
    }
    const uint64_t length = static_cast<uint64_t>(minimumLength);
    if (length > std::numeric_limits<uint64_t>::max() / length) {
        throw std::overflow_error("novel-anchor minimum length is too large to square");
    }
    return length * length;
}

// Test refLength * queryLength > minimumArea without multiplying the two
// sequence lengths.  This keeps the comparison defined even for unusually
// large inputs.
inline bool novelAnchorAreaExceeds(uint64_t refLength,
                                   uint64_t queryLength,
                                   uint64_t minimumArea) {
    return queryLength != 0 && refLength > minimumArea / queryLength;
}

}  // namespace anchorwave
