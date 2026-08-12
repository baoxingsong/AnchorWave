#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

namespace anchorwave {
namespace trio {

// Scores are maximized. Gap parameters use g(l) = open + l * extend, exactly
// matching the convention documented by TrioAnchorGraph. All gap values are
// normally non-positive. This bounded-window DP is also the test oracle for a
// future WFA/graph adapter.
struct TwoPieceAffineScoring {
    double match = 2.0;
    double mismatch = -4.0;
    double gapOpen1 = -4.0;
    double gapExtend1 = -2.0;
    double gapOpen2 = -80.0;
    double gapExtend2 = -1.0;

    void validate() const;
    double gapScore(std::size_t length) const;
};

struct PairwiseAlignment {
    std::string first;
    std::string second;
    double score = 0.0;
};

class TwoPieceAffineError : public std::runtime_error {
public:
    explicit TwoPieceAffineError(const std::string &message)
        : std::runtime_error(message) {}
};

class TwoPieceAffineAligner {
public:
    static PairwiseAlignment global(const std::string &first,
                                    const std::string &second,
                                    const TwoPieceAffineScoring &scoring,
                                    std::size_t maximumMatrixCells = 5000000);

    static double scoreAlignedPair(const std::string &firstAligned,
                                   const std::string &secondAligned,
                                   const TwoPieceAffineScoring &scoring);
};

}  // namespace trio
}  // namespace anchorwave
