#pragma once

#include "src/trio/alignment/TwoPieceAffine.h"

#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

struct PoaSequence {
    std::string pathId;
    std::string sequence;
};

struct LocalizedPoaOptions {
    TwoPieceAffineScoring scoring;
    std::size_t maxSequences = 32;
    std::size_t maxSequenceLength = 10000;
    std::size_t maxTotalBases = 100000;
    std::size_t maxPairwiseDpCells = 5000000;
    double minimumScoreDelta = 0.0;
};

struct LocalizedPoaRequest {
    std::string conflictId;
    std::string immutableLeftSite;
    std::string immutableRightSite;
    std::vector<PoaSequence> sequences;

    // Optional imported MSA, keyed by path ID. If supplied, it is used only as
    // the baseline objective; the ungapped rows must equal sequences exactly.
    std::map<std::string, std::string> baselineRows;
};

struct LocalizedPoaPatch {
    std::string repairId;
    std::string conflictId;
    std::string immutableLeftSite;
    std::string immutableRightSite;
    std::vector<std::string> pathIds;
    std::vector<std::string> alignedRows;
    double baselineScore = 0.0;
    double repairedScore = 0.0;
    double scoreDelta = 0.0;
    bool accepted = false;
    std::string disposition;
    std::string alignmentHash;
};

class LocalizedPoaError : public std::runtime_error {
public:
    explicit LocalizedPoaError(const std::string &message)
        : std::runtime_error(message) {}
};

// A deterministic localized column-POA. Each profile column is a POA site and
// distinct bases are allele branches. Multiple guide seeds are evaluated and
// the best sum-of-pairs two-piece-affine result is selected. The method returns
// a candidate patch and never mutates the authoritative graph directly.
class LocalizedPoaRepair {
public:
    static LocalizedPoaPatch repair(const LocalizedPoaRequest &request,
                                    const LocalizedPoaOptions &options);
};

}  // namespace trio
}  // namespace anchorwave
