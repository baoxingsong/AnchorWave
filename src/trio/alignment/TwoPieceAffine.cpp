#include "src/trio/alignment/TwoPieceAffine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <vector>

namespace anchorwave {
namespace trio {
namespace {

enum State : unsigned char { MATCH = 0, GAP_SECOND_1 = 1, GAP_SECOND_2 = 2,
                             GAP_FIRST_1 = 3, GAP_FIRST_2 = 4, NO_STATE = 255 };
const std::size_t STATE_COUNT = 5;

struct Cell {
    std::array<double, STATE_COUNT> score;
    std::array<unsigned char, STATE_COUNT> previous;
    Cell() {
        score.fill(-std::numeric_limits<double>::infinity());
        previous.fill(NO_STATE);
    }
};

bool better(double candidate, double current) {
    return candidate > current + 1e-12;
}

void offer(Cell &cell, State state, double candidate, State previous) {
    if (better(candidate, cell.score[state])) {
        cell.score[state] = candidate;
        cell.previous[state] = previous;
    }
}

double substitution(char left, char right,
                    const TwoPieceAffineScoring &scoring) {
    return left == right ? scoring.match : scoring.mismatch;
}

}  // namespace

void TwoPieceAffineScoring::validate() const {
    const double values[] = {match, mismatch, gapOpen1, gapExtend1,
                             gapOpen2, gapExtend2};
    for (double value : values) {
        if (!std::isfinite(value)) {
            throw TwoPieceAffineError("two-piece-affine scores must be finite");
        }
    }
    if (gapOpen1 > 0.0 || gapExtend1 > 0.0 ||
        gapOpen2 > 0.0 || gapExtend2 > 0.0) {
        throw TwoPieceAffineError("gap open/extension scores must be non-positive");
    }
}

double TwoPieceAffineScoring::gapScore(std::size_t length) const {
    validate();
    if (length == 0) return 0.0;
    const double first = gapOpen1 + static_cast<double>(length) * gapExtend1;
    const double second = gapOpen2 + static_cast<double>(length) * gapExtend2;
    return std::max(first, second);
}

PairwiseAlignment TwoPieceAffineAligner::global(
    const std::string &first, const std::string &second,
    const TwoPieceAffineScoring &scoring, std::size_t maximumMatrixCells) {
    scoring.validate();
    if (maximumMatrixCells == 0) {
        throw TwoPieceAffineError("maximum alignment matrix cell count must be positive");
    }
    const std::size_t columns = second.size() + 1;
    if (first.size() > std::numeric_limits<std::size_t>::max() / columns) {
        throw TwoPieceAffineError("alignment matrix size overflows size_t");
    }
    const std::size_t rows = first.size() + 1;
    if (rows > std::numeric_limits<std::size_t>::max() / columns) {
        throw TwoPieceAffineError("alignment matrix size overflows size_t");
    }
    const std::size_t matrixCells = rows * columns;
    if (matrixCells > maximumMatrixCells) {
        std::ostringstream message;
        message << "alignment matrix requires " << matrixCells
                << " cells, exceeding configured maximum " << maximumMatrixCells;
        throw TwoPieceAffineError(message.str());
    }
    std::vector<Cell> matrix(matrixCells);
    const auto cell = [&](std::size_t i, std::size_t j) -> Cell & {
        return matrix[i * columns + j];
    };
    cell(0, 0).score[MATCH] = 0.0;

    for (std::size_t i = 0; i <= first.size(); ++i) {
        for (std::size_t j = 0; j <= second.size(); ++j) {
            if (i == 0 && j == 0) continue;
            Cell &current = cell(i, j);
            if (i > 0 && j > 0) {
                const Cell &diagonal = cell(i - 1, j - 1);
                for (unsigned char state = 0; state < STATE_COUNT; ++state) {
                    offer(current, MATCH,
                          diagonal.score[state] +
                              substitution(first[i - 1], second[j - 1], scoring),
                          static_cast<State>(state));
                }
            }
            if (i > 0) {
                const Cell &above = cell(i - 1, j);
                offer(current, GAP_SECOND_1,
                      above.score[GAP_SECOND_1] + scoring.gapExtend1,
                      GAP_SECOND_1);
                offer(current, GAP_SECOND_2,
                      above.score[GAP_SECOND_2] + scoring.gapExtend2,
                      GAP_SECOND_2);
                const State openingSources[] = {MATCH, GAP_FIRST_1, GAP_FIRST_2};
                for (State source : openingSources) {
                    offer(current, GAP_SECOND_1,
                          above.score[source] + scoring.gapOpen1 + scoring.gapExtend1,
                          source);
                    offer(current, GAP_SECOND_2,
                          above.score[source] + scoring.gapOpen2 + scoring.gapExtend2,
                          source);
                }
            }
            if (j > 0) {
                const Cell &left = cell(i, j - 1);
                offer(current, GAP_FIRST_1,
                      left.score[GAP_FIRST_1] + scoring.gapExtend1,
                      GAP_FIRST_1);
                offer(current, GAP_FIRST_2,
                      left.score[GAP_FIRST_2] + scoring.gapExtend2,
                      GAP_FIRST_2);
                const State openingSources[] = {MATCH, GAP_SECOND_1, GAP_SECOND_2};
                for (State source : openingSources) {
                    offer(current, GAP_FIRST_1,
                          left.score[source] + scoring.gapOpen1 + scoring.gapExtend1,
                          source);
                    offer(current, GAP_FIRST_2,
                          left.score[source] + scoring.gapOpen2 + scoring.gapExtend2,
                          source);
                }
            }
        }
    }

    const Cell &last = cell(first.size(), second.size());
    State state = MATCH;
    for (unsigned char candidate = 1; candidate < STATE_COUNT; ++candidate) {
        if (better(last.score[candidate], last.score[state])) {
            state = static_cast<State>(candidate);
        }
    }
    if (!std::isfinite(last.score[state])) {
        throw TwoPieceAffineError("no finite global alignment was found");
    }
    PairwiseAlignment result;
    result.score = last.score[state];
    std::size_t i = first.size();
    std::size_t j = second.size();
    while (i > 0 || j > 0) {
        const Cell &current = cell(i, j);
        const State previous = static_cast<State>(current.previous[state]);
        if (state == MATCH) {
            if (i == 0 || j == 0) throw TwoPieceAffineError("invalid MATCH backtrace");
            result.first.push_back(first[--i]);
            result.second.push_back(second[--j]);
        } else if (state == GAP_SECOND_1 || state == GAP_SECOND_2) {
            if (i == 0) throw TwoPieceAffineError("invalid first-sequence gap backtrace");
            result.first.push_back(first[--i]);
            result.second.push_back('-');
        } else if (state == GAP_FIRST_1 || state == GAP_FIRST_2) {
            if (j == 0) throw TwoPieceAffineError("invalid second-sequence gap backtrace");
            result.first.push_back('-');
            result.second.push_back(second[--j]);
        } else {
            throw TwoPieceAffineError("invalid two-piece-affine backtrace state");
        }
        state = previous;
    }
    std::reverse(result.first.begin(), result.first.end());
    std::reverse(result.second.begin(), result.second.end());
    return result;
}

double TwoPieceAffineAligner::scoreAlignedPair(
    const std::string &firstAligned, const std::string &secondAligned,
    const TwoPieceAffineScoring &scoring) {
    scoring.validate();
    if (firstAligned.size() != secondAligned.size()) {
        throw TwoPieceAffineError("aligned strings have different widths");
    }
    double score = 0.0;
    std::size_t column = 0;
    while (column < firstAligned.size()) {
        if (firstAligned[column] == '-' && secondAligned[column] == '-') {
            ++column;
            continue;
        }
        if (firstAligned[column] != '-' && secondAligned[column] != '-') {
            score += substitution(firstAligned[column], secondAligned[column], scoring);
            ++column;
            continue;
        }
        const bool gapInFirst = firstAligned[column] == '-';
        std::size_t length = 0;
        while (column < firstAligned.size()) {
            if (firstAligned[column] == '-' && secondAligned[column] == '-') {
                ++column;
                continue;
            }
            const bool sameGap = gapInFirst ? firstAligned[column] == '-'
                                            : secondAligned[column] == '-';
            if (!sameGap || (gapInFirst ? secondAligned[column] == '-'
                                        : firstAligned[column] == '-')) {
                break;
            }
            ++length;
            ++column;
        }
        score += scoring.gapScore(length);
    }
    return score;
}

}  // namespace trio
}  // namespace anchorwave
