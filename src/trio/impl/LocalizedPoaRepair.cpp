#include "src/trio/impl/LocalizedPoaRepair.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

std::string ungapped(const std::string &sequence) {
    std::string result;
    for (char base : sequence) {
        if (base != '-') result.push_back(base);
    }
    return result;
}

char consensusAt(const std::vector<std::string> &rows, std::size_t column) {
    std::map<char, std::size_t> counts;
    for (const std::string &row : rows) {
        if (row[column] != '-') ++counts[row[column]];
    }
    if (counts.empty()) throw LocalizedPoaError("profile contains an all-gap column");
    char bestBase = counts.begin()->first;
    std::size_t bestCount = counts.begin()->second;
    for (const auto &entry : counts) {
        if (entry.second > bestCount ||
            (entry.second == bestCount && entry.first < bestBase)) {
            bestBase = entry.first;
            bestCount = entry.second;
        }
    }
    return bestBase;
}

std::string consensus(const std::vector<std::string> &rows) {
    if (rows.empty()) return std::string();
    std::string result;
    for (std::size_t column = 0; column < rows.front().size(); ++column) {
        result.push_back(consensusAt(rows, column));
    }
    return result;
}

void addSequenceToProfile(std::vector<std::string> &rows,
                          const std::string &sequence,
                          const TwoPieceAffineScoring &scoring,
                          std::size_t maximumMatrixCells) {
    if (rows.empty()) {
        rows.push_back(sequence);
        return;
    }
    const std::string profileConsensus = consensus(rows);
    PairwiseAlignment alignment;
    try {
        alignment = TwoPieceAffineAligner::global(
            profileConsensus, sequence, scoring, maximumMatrixCells);
    } catch (const TwoPieceAffineError &error) {
        throw LocalizedPoaError(std::string("localized POA resource limit: ") +
                                error.what());
    }
    std::vector<std::string> expanded(rows.size());
    std::string added;
    std::size_t oldColumn = 0;
    for (std::size_t column = 0; column < alignment.first.size(); ++column) {
        if (alignment.first[column] == '-') {
            for (std::string &row : expanded) row.push_back('-');
        } else {
            if (oldColumn >= rows.front().size()) {
                throw LocalizedPoaError("profile expansion exceeded old alignment width");
            }
            for (std::size_t row = 0; row < rows.size(); ++row) {
                expanded[row].push_back(rows[row][oldColumn]);
            }
            ++oldColumn;
        }
        added.push_back(alignment.second[column]);
    }
    if (oldColumn != rows.front().size()) {
        throw LocalizedPoaError("profile expansion did not consume all old columns");
    }
    rows.swap(expanded);
    rows.push_back(added);
}

double sumOfPairs(const std::vector<std::string> &rows,
                  const TwoPieceAffineScoring &scoring) {
    double result = 0.0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (std::size_t j = i + 1; j < rows.size(); ++j) {
            result += TwoPieceAffineAligner::scoreAlignedPair(rows[i], rows[j], scoring);
        }
    }
    return result;
}

void removeAllGapColumns(std::vector<std::string> &rows) {
    if (rows.empty()) return;
    std::vector<bool> keep(rows.front().size(), false);
    for (std::size_t column = 0; column < keep.size(); ++column) {
        for (const std::string &row : rows) {
            if (row[column] != '-') {
                keep[column] = true;
                break;
            }
        }
    }
    for (std::string &row : rows) {
        std::string compact;
        for (std::size_t column = 0; column < keep.size(); ++column) {
            if (keep[column]) compact.push_back(row[column]);
        }
        row.swap(compact);
    }
}

std::string serialize(const std::vector<std::string> &pathIds,
                      const std::vector<std::string> &rows) {
    std::ostringstream result;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        result << pathIds[i].size() << ':' << pathIds[i] << '='
               << rows[i].size() << ':' << rows[i] << ';';
    }
    return result.str();
}

std::vector<std::string> rowsInSortedPathOrder(
    const std::vector<std::string> &guideNames,
    const std::vector<std::string> &guideRows,
    const std::vector<std::string> &sortedNames) {
    std::vector<std::string> result;
    for (const std::string &name : sortedNames) {
        const auto found = std::find(guideNames.begin(), guideNames.end(), name);
        if (found == guideNames.end()) throw LocalizedPoaError("guide row name was lost");
        result.push_back(guideRows[static_cast<std::size_t>(found - guideNames.begin())]);
    }
    return result;
}

double baselineScore(const LocalizedPoaRequest &request,
                     const std::vector<std::string> &sortedNames,
                     const TwoPieceAffineScoring &scoring,
                     bool &hasBaseline) {
    hasBaseline = !request.baselineRows.empty();
    if (!hasBaseline) return -std::numeric_limits<double>::infinity();
    if (request.baselineRows.size() != request.sequences.size()) {
        throw LocalizedPoaError("baseline MSA does not contain every repair path");
    }
    std::vector<std::string> rows;
    std::size_t width = std::numeric_limits<std::size_t>::max();
    for (const std::string &name : sortedNames) {
        const auto row = request.baselineRows.find(name);
        if (row == request.baselineRows.end()) {
            throw LocalizedPoaError("baseline MSA is missing path '" + name + "'");
        }
        if (width == std::numeric_limits<std::size_t>::max()) width = row->second.size();
        if (row->second.size() != width) {
            throw LocalizedPoaError("baseline MSA rows have different widths");
        }
        rows.push_back(row->second);
    }
    for (std::size_t i = 0; i < sortedNames.size(); ++i) {
        const auto sequence = std::find_if(
            request.sequences.begin(), request.sequences.end(),
            [&](const PoaSequence &value) { return value.pathId == sortedNames[i]; });
        if (sequence == request.sequences.end() || ungapped(rows[i]) != sequence->sequence) {
            throw LocalizedPoaError("baseline MSA does not spell source path '" +
                                    sortedNames[i] + "'");
        }
    }
    return sumOfPairs(rows, scoring);
}

}  // namespace

LocalizedPoaPatch LocalizedPoaRepair::repair(
    const LocalizedPoaRequest &request, const LocalizedPoaOptions &options) {
    options.scoring.validate();
    if (request.conflictId.empty() || request.immutableLeftSite.empty() ||
        request.immutableRightSite.empty() ||
        request.immutableLeftSite == request.immutableRightSite) {
        throw LocalizedPoaError("repair requires a conflict ID and two distinct pinned sites");
    }
    if (request.sequences.size() < 2 ||
        request.sequences.size() > options.maxSequences) {
        throw LocalizedPoaError("repair sequence count is outside configured bounds");
    }
    if (options.maxPairwiseDpCells == 0) {
        throw LocalizedPoaError("maximum pairwise DP cell count must be positive");
    }
    std::vector<PoaSequence> sequences = request.sequences;
    std::sort(sequences.begin(), sequences.end(),
              [](const PoaSequence &a, const PoaSequence &b) {
                  return a.pathId < b.pathId;
              });
    std::size_t totalBases = 0;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        if (sequences[i].pathId.empty()) throw LocalizedPoaError("empty POA path ID");
        if (i && sequences[i - 1].pathId == sequences[i].pathId) {
            throw LocalizedPoaError("duplicate POA path ID '" + sequences[i].pathId + "'");
        }
        if (sequences[i].sequence.size() > options.maxSequenceLength) {
            throw LocalizedPoaError("POA sequence exceeds maximum local-window length");
        }
        if (sequences[i].sequence.find('-') != std::string::npos) {
            throw LocalizedPoaError("POA input sequence must be ungapped");
        }
        if (sequences[i].sequence.size() >
            std::numeric_limits<std::size_t>::max() - totalBases) {
            throw LocalizedPoaError("POA total base count overflows size_t");
        }
        totalBases += sequences[i].sequence.size();
        if (totalBases > options.maxTotalBases) {
            throw LocalizedPoaError("POA repair exceeds maximum total bases");
        }
    }
    std::vector<std::string> sortedNames;
    for (const PoaSequence &sequence : sequences) sortedNames.push_back(sequence.pathId);

    std::vector<std::string> bestRows;
    double bestScore = -std::numeric_limits<double>::infinity();
    std::string bestSerialization;
    for (std::size_t seed = 0; seed < sequences.size(); ++seed) {
        std::vector<std::string> guideNames(1, sequences[seed].pathId);
        std::vector<std::string> guideRows(1, sequences[seed].sequence);
        for (std::size_t i = 0; i < sequences.size(); ++i) {
            if (i == seed) continue;
            addSequenceToProfile(guideRows, sequences[i].sequence, options.scoring,
                                 options.maxPairwiseDpCells);
            guideNames.push_back(sequences[i].pathId);
        }
        std::vector<std::string> candidate =
            rowsInSortedPathOrder(guideNames, guideRows, sortedNames);
        removeAllGapColumns(candidate);
        const double score = sumOfPairs(candidate, options.scoring);
        const std::string canonical = serialize(sortedNames, candidate);
        if (score > bestScore + 1e-12 ||
            (std::fabs(score - bestScore) <= 1e-12 &&
             (bestSerialization.empty() || canonical < bestSerialization))) {
            bestScore = score;
            bestRows = candidate;
            bestSerialization = canonical;
        }
    }

    bool hasBaseline = false;
    const double oldScore = baselineScore(request, sortedNames, options.scoring, hasBaseline);
    LocalizedPoaPatch patch;
    patch.conflictId = request.conflictId;
    patch.immutableLeftSite = request.immutableLeftSite;
    patch.immutableRightSite = request.immutableRightSite;
    patch.pathIds = sortedNames;
    patch.alignedRows = bestRows;
    patch.baselineScore = oldScore;
    patch.repairedScore = bestScore;
    patch.scoreDelta = hasBaseline ? bestScore - oldScore
                                   : std::numeric_limits<double>::infinity();
    patch.accepted = !hasBaseline || patch.scoreDelta + 1e-12 >= options.minimumScoreDelta;
    patch.disposition = patch.accepted ? "candidate_patch_accepted"
                                       : "candidate_patch_below_minimum_delta";
    patch.alignmentHash = stableId("poa", {bestSerialization});
    patch.repairId = stableId("repair", {request.conflictId,
                                          request.immutableLeftSite,
                                          request.immutableRightSite,
                                          patch.alignmentHash});
    return patch;
}

}  // namespace trio
}  // namespace anchorwave
