#include "src/trio/model/AlignmentEvidence.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <tuple>

namespace anchorwave {
namespace trio {
namespace {

char normalizedBase(char base) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
}

char complementBase(char base) {
    switch (normalizedBase(base)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': case 'U': return 'A';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'S': return 'S';
        case 'W': return 'W';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'V': return 'B';
        case 'N': case 'X': case '.': case '*': return normalizedBase(base);
        default: return normalizedBase(base);
    }
}

char forwardBase(const MafRow &row, char alignmentBase) {
    // MAF text is written in the orientation declared by the row. ResidueId
    // coordinates are normalized to the forward source orientation, so its
    // spelling must be normalized too. Otherwise the same residue seen on '+'
    // in one pair and '-' in another appears to have conflicting bases.
    return row.strand == '-' ? complementBase(alignmentBase)
                             : normalizedBase(alignmentBase);
}

double blockScore(const MafBlock &block) {
    const auto found = block.attributes.find("score");
    if (found == block.attributes.end() || found->second.empty()) {
        return 0.0;
    }
    char *end = nullptr;
    const double value = std::strtod(found->second.c_str(), &end);
    if (end == found->second.c_str() || *end != '\0') {
        return 0.0;
    }
    return value;
}

std::string occurrencePath(const MafRow &row) {
    // This is deliberately a raw occurrence path, not an inferred copy ID.
    // Explicit copy relations can later split or group these source paths.
    return row.taxon + "|" + row.source;
}

ResidueId residueAt(const MafRow &row, int64_t position0) {
    ResidueId result;
    result.taxon = row.taxon;
    result.occurrencePath = occurrencePath(row);
    result.sequence = row.source;
    result.forwardPosition0 = position0;
    return result;
}

void addResidue(const ResidueId &id, char base, int64_t sourceSize,
                AlignmentEvidenceSet &destination, const std::string &where) {
    const char normalized = normalizedBase(base);
    const auto found = destination.residues.find(id);
    if (found != destination.residues.end() && found->second.base != normalized &&
        found->second.base != 'N' && normalized != 'N') {
        std::ostringstream message;
        message << where << ": source residue " << id.canonicalString()
                << " has conflicting bases '" << found->second.base
                << "' and '" << normalized << "' across MAF evidence";
        throw AlignmentEvidenceError(message.str());
    }
    if (found != destination.residues.end() && found->second.sourceSize != sourceSize) {
        throw AlignmentEvidenceError(where + ": source size changes for residue " +
                                     id.canonicalString());
    }
    ResidueObservation observation;
    observation.id = id;
    observation.base = found == destination.residues.end() || found->second.base == 'N'
                           ? normalized
                           : found->second.base;
    observation.sourceSize = sourceSize;
    destination.residues[id] = observation;
}

std::vector<int64_t> positionsByColumn(const MafRow &row) {
    std::vector<int64_t> positions(row.text.size(), -1);
    int64_t offset = 0;
    for (std::size_t column = 0; column < row.text.size(); ++column) {
        if (row.text[column] == '-') {
            continue;
        }
        positions[column] = row.strand == '+'
                                ? row.start0 + offset
                                : row.sourceSize - row.start0 - 1 - offset;
        ++offset;
    }
    return positions;
}

void registerSourceSize(const MafRow &row, AlignmentEvidenceSet &destination,
                        const std::string &where) {
    const std::pair<std::string, std::string> key(row.taxon, row.source);
    const auto found = destination.sourceSizes.find(key);
    if (found != destination.sourceSizes.end() && found->second != row.sourceSize) {
        std::ostringstream message;
        message << where << ": MAF source size for " << row.taxon << '|'
                << row.source << " changes from " << found->second << " to "
                << row.sourceSize;
        throw AlignmentEvidenceError(message.str());
    }
    destination.sourceSizes[key] = row.sourceSize;
}

int64_t nearestLeft(const std::vector<int64_t> &positions, std::size_t column) {
    for (std::size_t i = column; i > 0; --i) {
        if (positions[i - 1] >= 0) {
            return positions[i - 1];
        }
    }
    return -1;
}

int64_t nearestRight(const std::vector<int64_t> &positions, std::size_t column) {
    for (std::size_t i = column + 1; i < positions.size(); ++i) {
        if (positions[i] >= 0) {
            return positions[i];
        }
    }
    return -1;
}

EvidenceProvenance provenanceFor(const PairwiseMafInput &input,
                                 const std::string &pairId,
                                 const MafBlock &block, std::size_t column,
                                 const ResidueId &left, const std::string &right) {
    EvidenceProvenance provenance;
    provenance.pairId = pairId;
    provenance.runId = input.runId.empty() ? input.mafPath : input.runId;
    provenance.mafPath = input.mafPath;
    provenance.blockIndex = block.blockIndex;
    provenance.column = column;
    provenance.alignmentScore = blockScore(block) * input.sourceWeight;
    provenance.evidenceId = stableId(
        "ev", {provenance.runId, pairId, std::to_string(block.blockIndex),
               std::to_string(column), left.canonicalString(), right});
    return provenance;
}

}  // namespace

bool ResidueId::operator<(const ResidueId &other) const {
    return std::tie(taxon, occurrencePath, sequence, forwardPosition0) <
           std::tie(other.taxon, other.occurrencePath, other.sequence,
                    other.forwardPosition0);
}

bool ResidueId::operator==(const ResidueId &other) const {
    return taxon == other.taxon && occurrencePath == other.occurrencePath &&
           sequence == other.sequence && forwardPosition0 == other.forwardPosition0;
}

std::string ResidueId::canonicalString() const {
    std::ostringstream value;
    value << taxon << '|' << occurrencePath << '|' << sequence << ':' << forwardPosition0;
    return value.str();
}

std::string PairwiseEvidenceNormalizer::canonicalPairId(
    const std::string &taxonA, const std::string &taxonB) {
    if (taxonA.empty() || taxonB.empty() || taxonA == taxonB) {
        throw AlignmentEvidenceError("pairwise evidence requires two distinct non-empty taxa");
    }
    return taxonA < taxonB ? taxonA + "--" + taxonB : taxonB + "--" + taxonA;
}

AlignmentEvidenceSet PairwiseEvidenceNormalizer::normalize(
    const std::vector<PairwiseMafInput> &inputs) {
    AlignmentEvidenceSet result;
    std::set<std::string> runIds;
    for (const PairwiseMafInput &input : inputs) {
        const std::string pairId = canonicalPairId(input.leftTaxon, input.rightTaxon);
        const std::string runId = input.runId.empty() ? input.mafPath : input.runId;
        if (runId.empty()) {
            throw AlignmentEvidenceError(pairId + ": empty MAF path and run ID");
        }
        if (!runIds.insert(runId).second) {
            throw AlignmentEvidenceError("duplicate pairwise run ID '" + runId + "'");
        }
        const std::vector<MafBlock> blocks = MafEvidenceReader::readFile(
            input.mafPath, input.leftTaxon, input.rightTaxon);
        appendBlocks(input, blocks, result);
    }
    return result;
}

void PairwiseEvidenceNormalizer::appendBlocks(
    const PairwiseMafInput &input, const std::vector<MafBlock> &blocks,
    AlignmentEvidenceSet &destination) {
    const std::string pairId = canonicalPairId(input.leftTaxon, input.rightTaxon);
    destination.observedPairs.insert(pairId);
    for (const MafBlock &block : blocks) {
        if (block.rows.size() != 2) {
            throw AlignmentEvidenceError(input.mafPath + ": normalized pairwise block is not two-row");
        }
        const MafRow &leftRow = block.rows[0];
        const MafRow &rightRow = block.rows[1];
        registerSourceSize(leftRow, destination, input.mafPath);
        registerSourceSize(rightRow, destination, input.mafPath);
        const std::vector<int64_t> leftPositions = positionsByColumn(leftRow);
        const std::vector<int64_t> rightPositions = positionsByColumn(rightRow);
        for (std::size_t column = 0; column < leftRow.text.size(); ++column) {
            const bool hasLeft = leftPositions[column] >= 0;
            const bool hasRight = rightPositions[column] >= 0;
            if (!hasLeft && !hasRight) {
                throw AlignmentEvidenceError(input.mafPath + ": double-gap MAF column is invalid");
            }
            if (hasLeft) {
                addResidue(residueAt(leftRow, leftPositions[column]),
                           forwardBase(leftRow, leftRow.text[column]),
                           leftRow.sourceSize, destination, input.mafPath);
            }
            if (hasRight) {
                addResidue(residueAt(rightRow, rightPositions[column]),
                           forwardBase(rightRow, rightRow.text[column]),
                           rightRow.sourceSize, destination, input.mafPath);
            }
            if (hasLeft && hasRight) {
                HomologyEvidence evidence;
                evidence.left = residueAt(leftRow, leftPositions[column]);
                evidence.right = residueAt(rightRow, rightPositions[column]);
                const std::string rightKey = evidence.right.canonicalString();
                evidence.provenance = provenanceFor(input, pairId, block, column,
                                                     evidence.left, rightKey);
                destination.homologies.push_back(evidence);
            } else {
                const MafRow &presentRow = hasLeft ? leftRow : rightRow;
                const MafRow &absentRow = hasLeft ? rightRow : leftRow;
                const std::vector<int64_t> &presentPositions = hasLeft ? leftPositions : rightPositions;
                const std::vector<int64_t> &absentPositions = hasLeft ? rightPositions : leftPositions;
                AlignedAbsenceEvidence evidence;
                evidence.present = residueAt(presentRow, presentPositions[column]);
                evidence.absentTaxon = absentRow.taxon;
                evidence.absentOccurrencePath = occurrencePath(absentRow);
                evidence.absentSequence = absentRow.source;
                evidence.absentLeftFlank0 = nearestLeft(absentPositions, column);
                evidence.absentRightFlank0 = nearestRight(absentPositions, column);
                const std::string rightKey = evidence.absentTaxon + "|gap|" +
                                             std::to_string(evidence.absentLeftFlank0) + "|" +
                                             std::to_string(evidence.absentRightFlank0);
                evidence.provenance = provenanceFor(input, pairId, block, column,
                                                     evidence.present, rightKey);
                destination.alignedAbsences.push_back(evidence);
            }
        }
    }
    std::sort(destination.homologies.begin(), destination.homologies.end(),
              [](const HomologyEvidence &a, const HomologyEvidence &b) {
                  return a.provenance.evidenceId < b.provenance.evidenceId;
              });
    std::sort(destination.alignedAbsences.begin(), destination.alignedAbsences.end(),
              [](const AlignedAbsenceEvidence &a, const AlignedAbsenceEvidence &b) {
                  return a.provenance.evidenceId < b.provenance.evidenceId;
              });
}

}  // namespace trio
}  // namespace anchorwave
