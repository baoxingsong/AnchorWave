#include "src/trio/impl/CopyRelationshipResolver.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <tuple>

namespace anchorwave {
namespace trio {
namespace {

struct Candidate {
    std::string group;
    ConstraintStrength strength;
    CopyRelation relation;
    double confidence;
    std::string provenance;
};

struct ResidueConstraintState {
    std::vector<Candidate> candidates;
    std::set<std::string> excludedGroups;
    std::string selectedGroup;
    bool selectedHard = false;
    double selectedConfidence = 0.0;
    std::string provenance;
};

std::string groupId(const CopyConstraint &constraint) {
    return constraint.familyId + ":" + constraint.ancestralCopyId;
}

std::string constraintProvenance(const CopyConstraint &constraint) {
    std::ostringstream result;
    result << constraint.sourceLocation.path << ':' << constraint.sourceLocation.line
           << ':' << constraint.source;
    return result.str();
}

std::string residueKey(const ResidueId &residue) {
    return residue.taxon + "\t" + residue.sequence;
}

std::string memberKey(const CopyConstraint &constraint) {
    std::ostringstream result;
    result << constraint.familyId << '\t' << constraint.taxonId << '\t'
           << static_cast<int>(constraint.memberType) << '\t';
    if (constraint.memberType == CopyMemberType::ANCHOR) {
        result << constraint.memberId;
    } else {
        result << constraint.interval.sequence << ':' << constraint.interval.start0
               << '-' << constraint.interval.end0 << ':' << constraint.interval.strand;
    }
    return result.str();
}

void rejectUnexpandedAnchorConstraints(const CopyConstraintSet &constraints) {
    for (const CopyConstraint &constraint : constraints.records) {
        if (constraint.recordType != CopyRecordType::MEMBER ||
            constraint.memberType != CopyMemberType::ANCHOR) {
            continue;
        }
        std::ostringstream message;
        message << "anchor MEMBER constraint";
        if (!constraint.memberId.empty()) {
            message << " '" << constraint.memberId << "'";
        }
        if (!constraint.taxonId.empty()) {
            message << " for taxon '" << constraint.taxonId << "'";
        }
        if (!constraint.sourceLocation.path.empty()) {
            message << " at " << constraint.sourceLocation.path;
            if (constraint.sourceLocation.line != 0) {
                message << ':' << constraint.sourceLocation.line;
            }
        }
        message << " cannot be resolved because anchor-to-interval expansion is not "
                   "implemented; provide an equivalent member_type=interval constraint";
        throw CopyResolutionError(message.str());
    }
}

void validateCapacities(const CopyConstraintSet &constraints) {
    std::map<std::pair<std::string, std::string>, std::set<std::string>> hardMembers;
    std::map<std::pair<std::string, std::string>, int64_t> hardCounts;
    for (const CopyConstraint &constraint : constraints.records) {
        const std::pair<std::string, std::string> key(constraint.familyId,
                                                      constraint.taxonId);
        if (constraint.recordType == CopyRecordType::MEMBER &&
            constraint.strength == ConstraintStrength::HARD &&
            constraint.relation != CopyRelation::EXCLUDED) {
            hardMembers[key].insert(memberKey(constraint));
        }
        if (constraint.recordType == CopyRecordType::COUNT &&
            constraint.strength == ConstraintStrength::HARD) {
            hardCounts[key] = constraint.expectedCount;
        }
    }
    for (const auto &entry : hardMembers) {
        auto count = hardCounts.find(entry.first);
        if (count != hardCounts.end() &&
            static_cast<int64_t>(entry.second.size()) > count->second) {
            std::ostringstream message;
            message << "hard MEMBER occurrences exceed hard COUNT for family '"
                    << entry.first.first << "', taxon '" << entry.first.second << "'";
            throw CopyResolutionError(message.str());
        }
        int64_t fallback = 0;
        if (constraints.findDefaultCapacity(entry.first.first, entry.first.second,
                                            fallback) &&
            static_cast<int64_t>(entry.second.size()) > fallback) {
            std::ostringstream message;
            message << "hard MEMBER occurrences exceed default capacity for family '"
                    << entry.first.first << "', taxon '" << entry.first.second << "'";
            throw CopyResolutionError(message.str());
        }
    }
}

std::map<ResidueId, ResidueConstraintState> mapIntervalConstraints(
    const AlignmentEvidenceSet &evidence, const CopyConstraintSet &constraints) {
    std::map<ResidueId, ResidueConstraintState> result;
    for (const auto &entry : evidence.residues) result[entry.first];

    std::map<std::string, std::vector<ResidueId>> residuesBySequence;
    for (const auto &entry : evidence.residues) {
        residuesBySequence[residueKey(entry.first)].push_back(entry.first);
    }
    for (auto &entry : residuesBySequence) {
        std::sort(entry.second.begin(), entry.second.end(),
                  [](const ResidueId &a, const ResidueId &b) {
                      return a.forwardPosition0 < b.forwardPosition0;
                  });
    }

    for (const CopyConstraint &constraint : constraints.records) {
        if (constraint.recordType != CopyRecordType::MEMBER ||
            constraint.memberType != CopyMemberType::INTERVAL) {
            continue;
        }
        ResidueId keyResidue;
        keyResidue.taxon = constraint.taxonId;
        keyResidue.sequence = constraint.interval.sequence;
        const auto sequence = residuesBySequence.find(residueKey(keyResidue));
        if (sequence == residuesBySequence.end()) continue;
        const auto begin = std::lower_bound(
            sequence->second.begin(), sequence->second.end(), constraint.interval.start0,
            [](const ResidueId &residue, int64_t coordinate) {
                return residue.forwardPosition0 < coordinate;
            });
        for (auto residue = begin; residue != sequence->second.end() &&
                                    residue->forwardPosition0 < constraint.interval.end0;
             ++residue) {
            Candidate candidate;
            candidate.group = groupId(constraint);
            candidate.strength = constraint.strength;
            candidate.relation = constraint.relation;
            candidate.confidence = constraint.confidence;
            candidate.provenance = constraintProvenance(constraint);
            if (constraint.relation == CopyRelation::EXCLUDED) {
                if (constraint.strength == ConstraintStrength::HARD) {
                    result[*residue].excludedGroups.insert(candidate.group);
                }
            } else {
                result[*residue].candidates.push_back(candidate);
            }
        }
    }
    return result;
}

void selectExplicitCandidates(std::map<ResidueId, ResidueConstraintState> &states,
                              double softTieMargin,
                              std::vector<ResidueCopyAlternative> &alternatives) {
    for (auto &entry : states) {
        ResidueConstraintState &state = entry.second;
        std::map<std::string, double> hardScores;
        std::map<std::string, double> softScores;
        std::map<std::string, std::string> provenance;
        for (const Candidate &candidate : state.candidates) {
            std::map<std::string, double> &scores =
                candidate.strength == ConstraintStrength::HARD ? hardScores : softScores;
            scores[candidate.group] += candidate.confidence;
            provenance[candidate.group] = candidate.provenance;
        }
        if (hardScores.size() > 1) {
            throw CopyResolutionError("residue " + entry.first.canonicalString() +
                                      " belongs to incompatible overlapping hard copy groups");
        }
        if (!hardScores.empty()) {
            state.selectedGroup = hardScores.begin()->first;
            state.selectedHard = true;
            state.selectedConfidence = hardScores.begin()->second;
            state.provenance = provenance[state.selectedGroup];
            if (state.excludedGroups.count(state.selectedGroup)) {
                throw CopyResolutionError("residue " + entry.first.canonicalString() +
                                          " is both included in and excluded from a hard group");
            }
            continue;
        }
        if (softScores.empty()) continue;
        std::vector<std::pair<std::string, double>> ranked(softScores.begin(),
                                                           softScores.end());
        std::sort(ranked.begin(), ranked.end(),
                  [](const std::pair<std::string, double> &a,
                     const std::pair<std::string, double> &b) {
                      if (a.second != b.second) return a.second > b.second;
                      return a.first < b.first;
                  });
        if (ranked.size() > 1 &&
            ranked[0].second - ranked[1].second <= softTieMargin + 1e-12) {
            ResidueCopyAlternative alternative;
            alternative.residue = entry.first;
            for (const auto &candidate : ranked) {
                alternative.copyGroups.push_back(candidate.first);
                alternative.scores.push_back(candidate.second);
            }
            alternatives.push_back(alternative);
            continue;
        }
        state.selectedGroup = ranked.front().first;
        state.selectedConfidence = ranked.front().second;
        state.provenance = provenance[state.selectedGroup];
    }
}

class ResolutionSet {
public:
    ResolutionSet(const std::vector<ResidueId> &residues,
                  const std::map<ResidueId, ResidueConstraintState> &states)
        : parent_(residues.size()), rank_(residues.size(), 0),
          occurrences_(residues.size()), groups_(residues.size()),
          exclusions_(residues.size()) {
        std::iota(parent_.begin(), parent_.end(), 0);
        for (std::size_t i = 0; i < residues.size(); ++i) {
            occurrences_[i][residues[i].occurrencePath] = residues[i];
            const auto state = states.find(residues[i]);
            if (state != states.end()) {
                if (!state->second.selectedGroup.empty()) {
                    groups_[i].insert(state->second.selectedGroup);
                }
                exclusions_[i] = state->second.excludedGroups;
            }
        }
    }
    std::size_t find(std::size_t value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }
    bool compatible(std::size_t left, std::size_t right, std::string &reason) {
        left = find(left); right = find(right);
        if (left == right) return true;
        for (const auto &entry : occurrences_[left]) {
            const auto other = occurrences_[right].find(entry.first);
            if (other != occurrences_[right].end() && !(entry.second == other->second)) {
                reason = "one-to-many collision on occurrence path '" + entry.first + "'";
                return false;
            }
        }
        if (!groups_[left].empty() && !groups_[right].empty() &&
            groups_[left] != groups_[right]) {
            reason = "candidate edge joins incompatible explicit copy groups";
            return false;
        }
        for (const std::string &group : groups_[left]) {
            if (exclusions_[right].count(group)) {
                reason = "candidate edge contradicts a hard copy exclusion";
                return false;
            }
        }
        for (const std::string &group : groups_[right]) {
            if (exclusions_[left].count(group)) {
                reason = "candidate edge contradicts a hard copy exclusion";
                return false;
            }
        }
        return true;
    }
    void join(std::size_t left, std::size_t right) {
        left = find(left); right = find(right);
        if (left == right) return;
        if (rank_[left] < rank_[right]) std::swap(left, right);
        parent_[right] = left;
        occurrences_[left].insert(occurrences_[right].begin(), occurrences_[right].end());
        groups_[left].insert(groups_[right].begin(), groups_[right].end());
        exclusions_[left].insert(exclusions_[right].begin(), exclusions_[right].end());
        occurrences_[right].clear(); groups_[right].clear(); exclusions_[right].clear();
        if (rank_[left] == rank_[right]) ++rank_[left];
    }
    const std::set<std::string> &groups(std::size_t value) {
        return groups_[find(value)];
    }
private:
    std::vector<std::size_t> parent_;
    std::vector<unsigned char> rank_;
    std::vector<std::map<std::string, ResidueId>> occurrences_;
    std::vector<std::set<std::string>> groups_;
    std::vector<std::set<std::string>> exclusions_;
};

}  // namespace

CopyResolutionResult CopyRelationshipResolver::resolve(
    const AlignmentEvidenceSet &evidence, const CopyConstraintSet &constraints,
    const CopyResolverOptions &options) {
    if (options.softTieMargin < 0.0 || !std::isfinite(options.softTieMargin)) {
        throw CopyResolutionError("soft copy tie margin must be finite and non-negative");
    }
    rejectUnexpandedAnchorConstraints(constraints);
    validateCapacities(constraints);
    std::map<ResidueId, ResidueConstraintState> states =
        mapIntervalConstraints(evidence, constraints);
    CopyResolutionResult result;
    selectExplicitCandidates(states, options.softTieMargin, result.alternatives);
    if (!result.alternatives.empty()) result.status = CopyResolutionStatus::ALTERNATIVES;

    std::vector<ResidueId> residues;
    std::map<ResidueId, std::size_t> index;
    for (const auto &entry : evidence.residues) {
        index[entry.first] = residues.size();
        residues.push_back(entry.first);
    }
    if (options.strictExplicitGroups) {
        for (const ResidueId &residue : residues) {
            if (states[residue].selectedGroup.empty() || !states[residue].selectedHard) {
                throw CopyResolutionError("strict copy mode requires a hard explicit group for " +
                                          residue.canonicalString());
            }
        }
    }

    ResolutionSet components(residues, states);
    std::vector<HomologyEvidence> edges = evidence.homologies;
    std::sort(edges.begin(), edges.end(),
              [&](const HomologyEvidence &a, const HomologyEvidence &b) {
                  const ResidueConstraintState &al = states[a.left];
                  const ResidueConstraintState &ar = states[a.right];
                  const ResidueConstraintState &bl = states[b.left];
                  const ResidueConstraintState &br = states[b.right];
                  const int explicitA = !al.selectedGroup.empty() + !ar.selectedGroup.empty();
                  const int explicitB = !bl.selectedGroup.empty() + !br.selectedGroup.empty();
                  if (explicitA != explicitB) return explicitA > explicitB;
                  if (a.provenance.alignmentScore != b.provenance.alignmentScore) {
                      return a.provenance.alignmentScore > b.provenance.alignmentScore;
                  }
                  return a.provenance.evidenceId < b.provenance.evidenceId;
              });
    for (const HomologyEvidence &edge : edges) {
        CopyResolutionAudit audit;
        audit.evidenceId = edge.provenance.evidenceId;
        std::string reason;
        if (!components.compatible(index.at(edge.left), index.at(edge.right), reason)) {
            audit.classification = "COPY_COLLISION";
            audit.reason = reason;
            result.audit.push_back(audit);
            result.status = CopyResolutionStatus::ALTERNATIVES;
            continue;
        }
        components.join(index.at(edge.left), index.at(edge.right));
        audit.selected = true;
        audit.classification = "SELECTED";
        audit.reason = "copy-compatible pairwise relationship";
        result.audit.push_back(audit);
        result.selectedEvidenceIds.insert(edge.provenance.evidenceId);
    }

    std::map<std::size_t, std::vector<ResidueId>> componentResidues;
    for (std::size_t i = 0; i < residues.size(); ++i) {
        componentResidues[components.find(i)].push_back(residues[i]);
    }
    for (auto &entry : componentResidues) {
        std::sort(entry.second.begin(), entry.second.end());
        const std::set<std::string> explicitGroups = components.groups(entry.first);
        std::string selectedGroup;
        bool hard = false;
        double confidence = 0.0;
        std::string provenance;
        if (!explicitGroups.empty()) {
            selectedGroup = *explicitGroups.begin();
            for (const ResidueId &residue : entry.second) {
                const ResidueConstraintState &state = states[residue];
                if (state.selectedGroup == selectedGroup) {
                    hard |= state.selectedHard;
                    confidence = std::max(confidence, state.selectedConfidence);
                    if (!state.provenance.empty()) provenance = state.provenance;
                }
            }
        } else {
            std::vector<std::string> keys;
            for (const ResidueId &residue : entry.second) keys.push_back(residue.canonicalString());
            selectedGroup = stableId("inferred_copy", keys);
            provenance = "unambiguous_copy-compatible_pairwise_component";
            confidence = entry.second.size() >= 3 ? 1.0 : 0.5;
        }
        for (const ResidueId &residue : entry.second) {
            CopyAssignment assignment;
            assignment.copyGroup = selectedGroup;
            assignment.hard = hard;
            assignment.confidence = confidence;
            assignment.provenance = provenance;
            result.assignments[residue] = assignment;
        }
    }
    std::sort(result.audit.begin(), result.audit.end(),
              [](const CopyResolutionAudit &a, const CopyResolutionAudit &b) {
                  return a.evidenceId < b.evidenceId;
              });
    return result;
}

}  // namespace trio
}  // namespace anchorwave
