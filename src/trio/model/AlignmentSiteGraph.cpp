#include "src/trio/model/AlignmentSiteGraph.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <tuple>

namespace anchorwave {
namespace trio {
namespace {

class DisjointSet {
public:
    explicit DisjointSet(const std::vector<ResidueId> &residues)
        : parent_(residues.size()), rank_(residues.size(), 0),
          occurrenceMembers_(residues.size()) {
        std::iota(parent_.begin(), parent_.end(), 0);
        for (std::size_t i = 0; i < residues.size(); ++i) {
            occurrenceMembers_[i][residues[i].occurrencePath] = residues[i];
        }
    }
    std::size_t find(std::size_t value) {
        if (parent_[value] != value) {
            parent_[value] = find(parent_[value]);
        }
        return parent_[value];
    }
    void join(std::size_t left, std::size_t right) {
        left = find(left);
        right = find(right);
        if (left == right) {
            return;
        }
        if (rank_[left] < rank_[right]) {
            std::swap(left, right);
        }
        parent_[right] = left;
        occurrenceMembers_[left].insert(occurrenceMembers_[right].begin(),
                                        occurrenceMembers_[right].end());
        occurrenceMembers_[right].clear();
        if (rank_[left] == rank_[right]) {
            ++rank_[left];
        }
    }
    bool collision(std::size_t left, std::size_t right,
                   std::string &occurrencePath) {
        left = find(left);
        right = find(right);
        if (left == right) {
            return false;
        }
        const std::map<std::string, ResidueId> *smaller = &occurrenceMembers_[left];
        const std::map<std::string, ResidueId> *larger = &occurrenceMembers_[right];
        if (smaller->size() > larger->size()) {
            std::swap(smaller, larger);
        }
        for (const auto &entry : *smaller) {
            const auto found = larger->find(entry.first);
            if (found != larger->end() && !(entry.second == found->second)) {
                occurrencePath = entry.first;
                return true;
            }
        }
        return false;
    }
private:
    std::vector<std::size_t> parent_;
    std::vector<unsigned char> rank_;
    std::vector<std::map<std::string, ResidueId>> occurrenceMembers_;
};

bool assignmentConflict(const ResidueId &left, const ResidueId &right,
                        const std::map<ResidueId, CopyAssignment> &assignments,
                        std::string &reason) {
    const auto a = assignments.find(left);
    const auto b = assignments.find(right);
    if (a == assignments.end() || b == assignments.end() ||
        a->second.copyGroup.empty() || b->second.copyGroup.empty() ||
        a->second.copyGroup == b->second.copyGroup) {
        return false;
    }
    if (a->second.hard || b->second.hard) {
        reason = "hard copy groups disagree: " + a->second.copyGroup + " versus " +
                 b->second.copyGroup;
        return true;
    }
    return false;
}

std::vector<std::string> requiredPairs(const AlignmentGraphBuildOptions &options) {
    return {
        PairwiseEvidenceNormalizer::canonicalPairId(options.ingroup1, options.ingroup2),
        PairwiseEvidenceNormalizer::canonicalPairId(options.ingroup1, options.primaryOutgroup),
        PairwiseEvidenceNormalizer::canonicalPairId(options.ingroup2, options.primaryOutgroup)};
}

ConsistencyClass classifySite(const AlignmentSite &site,
                              const AlignmentGraphBuildOptions &options) {
    std::set<std::string> taxa;
    for (const ResidueObservation &observation : site.observations) {
        taxa.insert(observation.id.taxon);
    }
    const bool allTaxa = taxa.count(options.ingroup1) && taxa.count(options.ingroup2) &&
                         taxa.count(options.primaryOutgroup);
    if (!allTaxa) {
        return ConsistencyClass::MISSING_EVIDENCE;
    }
    std::size_t supported = 0;
    for (const std::string &pair : requiredPairs(options)) {
        supported += site.directPairSupport.count(pair);
    }
    if (supported == 3) {
        return ConsistencyClass::CONSISTENT;
    }
    if (supported == 2) {
        return ConsistencyClass::SUPPORTED_BY_TWO_EDGES;
    }
    return ConsistencyClass::MISSING_EVIDENCE;
}

std::string residuesKey(std::vector<ResidueObservation> observations) {
    std::sort(observations.begin(), observations.end(),
              [](const ResidueObservation &a, const ResidueObservation &b) {
                  return a.id < b.id;
              });
    std::ostringstream result;
    for (const ResidueObservation &observation : observations) {
        result << observation.id.canonicalString() << '=' << observation.base << ';';
    }
    return result.str();
}

GraphConflict makeConflict(ConsistencyClass classification,
                           const std::vector<std::string> &evidenceIds,
                           const std::vector<ResidueId> &residues,
                           bool poaEligible,
                           const std::string &disposition) {
    std::vector<std::string> keys;
    keys.push_back(consistencyClassName(classification));
    keys.insert(keys.end(), evidenceIds.begin(), evidenceIds.end());
    std::vector<std::string> residueKeys;
    for (const ResidueId &residue : residues) {
        residueKeys.push_back(residue.canonicalString());
    }
    std::sort(residueKeys.begin(), residueKeys.end());
    keys.insert(keys.end(), residueKeys.begin(), residueKeys.end());
    GraphConflict conflict;
    conflict.conflictId = stableId("conflict", keys);
    conflict.classification = classification;
    conflict.evidenceIds = evidenceIds;
    conflict.residues = residues;
    conflict.poaEligible = poaEligible;
    conflict.disposition = disposition;
    return conflict;
}

}  // namespace

std::string consistencyClassName(ConsistencyClass value) {
    switch (value) {
        case ConsistencyClass::CONSISTENT: return "CONSISTENT";
        case ConsistencyClass::LOCAL_POA_REPAIRED: return "LOCAL_POA_REPAIRED";
        case ConsistencyClass::SUPPORTED_BY_TWO_EDGES: return "SUPPORTED_BY_TWO_EDGES";
        case ConsistencyClass::MISSING_EVIDENCE: return "MISSING_EVIDENCE";
        case ConsistencyClass::GAP_SHIFT: return "GAP_SHIFT";
        case ConsistencyClass::BOUNDARY_DISCORDANCE: return "BOUNDARY_DISCORDANCE";
        case ConsistencyClass::RESIDUE_TRANSITIVITY: return "RESIDUE_TRANSITIVITY";
        case ConsistencyClass::LOW_CONFIDENCE_PAIRWISE: return "LOW_CONFIDENCE_PAIRWISE";
        case ConsistencyClass::NONHOMOLOGOUS_INSERTIONS: return "NONHOMOLOGOUS_INSERTIONS";
        case ConsistencyClass::COPY_COLLISION: return "COPY_COLLISION";
        case ConsistencyClass::ORDER_ORIENTATION: return "ORDER_ORIENTATION";
        case ConsistencyClass::UNANCHORED_OR_TOO_LARGE: return "UNANCHORED_OR_TOO_LARGE";
    }
    return "UNKNOWN";
}

AlignmentSiteGraph AlignmentSiteGraphBuilder::build(
    const AlignmentEvidenceSet &evidence,
    const AlignmentGraphBuildOptions &options) {
    if (options.ingroup1.empty() || options.ingroup2.empty() ||
        options.primaryOutgroup.empty() || options.ingroup1 == options.ingroup2 ||
        options.ingroup1 == options.primaryOutgroup ||
        options.ingroup2 == options.primaryOutgroup) {
        throw AlignmentGraphError("two distinct ingroups and a distinct primary outgroup are required");
    }
    if (options.requirePrimaryTriangle) {
        for (const std::string &pair : requiredPairs(options)) {
            if (!evidence.observedPairs.count(pair)) {
                throw AlignmentGraphError("missing required primary-trio pairwise evidence: " + pair);
            }
        }
    }

    AlignmentSiteGraph graph;
    graph.ingroup1_ = options.ingroup1;
    graph.ingroup2_ = options.ingroup2;
    graph.primaryOutgroup_ = options.primaryOutgroup;

    std::vector<ResidueId> residues;
    residues.reserve(evidence.residues.size());
    std::map<ResidueId, std::size_t> index;
    for (const auto &entry : evidence.residues) {
        index[entry.first] = residues.size();
        residues.push_back(entry.first);
    }
    DisjointSet set(residues);
    std::vector<HomologyEvidence> edges = evidence.homologies;
    std::sort(edges.begin(), edges.end(),
              [](const HomologyEvidence &a, const HomologyEvidence &b) {
                  return a.provenance.evidenceId < b.provenance.evidenceId;
              });

    for (const HomologyEvidence &edge : edges) {
        EvidenceDisposition disposition;
        disposition.evidenceId = edge.provenance.evidenceId;
        const auto leftIndex = index.find(edge.left);
        const auto rightIndex = index.find(edge.right);
        if (leftIndex == index.end() || rightIndex == index.end()) {
            throw AlignmentGraphError("homology evidence references an unknown residue");
        }
        if (options.restrictToAllowedEvidence &&
            !options.allowedEvidenceIds.count(edge.provenance.evidenceId)) {
            disposition.classification = ConsistencyClass::COPY_COLLISION;
            disposition.reason = "edge was not selected by the copy resolver";
            graph.evidenceAudit_.push_back(disposition);
            continue;
        }
        if (options.requireResolvedCopyAssignments &&
            (!options.copyAssignments.count(edge.left) ||
             !options.copyAssignments.count(edge.right))) {
            disposition.classification = ConsistencyClass::COPY_COLLISION;
            disposition.reason = "copy relationship is unresolved before graph construction";
            graph.evidenceAudit_.push_back(disposition);
            continue;
        }
        std::string copyReason;
        if (assignmentConflict(edge.left, edge.right, options.copyAssignments, copyReason)) {
            disposition.classification = ConsistencyClass::COPY_COLLISION;
            disposition.reason = copyReason;
            graph.evidenceAudit_.push_back(disposition);
            graph.conflicts_.push_back(makeConflict(
                ConsistencyClass::COPY_COLLISION, {edge.provenance.evidenceId},
                {edge.left, edge.right}, false, "return_to_copy_resolver"));
            continue;
        }

        const std::size_t leftRoot = set.find(leftIndex->second);
        const std::size_t rightRoot = set.find(rightIndex->second);
        if (leftRoot != rightRoot) {
            std::string collisionPath;
            if (set.collision(leftRoot, rightRoot, collisionPath)) {
                disposition.classification = ConsistencyClass::RESIDUE_TRANSITIVITY;
                disposition.reason = "union would place two residues from occurrence path '" +
                                     collisionPath + "' in one site";
                graph.evidenceAudit_.push_back(disposition);
                graph.conflicts_.push_back(makeConflict(
                    ConsistencyClass::RESIDUE_TRANSITIVITY,
                    {edge.provenance.evidenceId}, {edge.left, edge.right}, true,
                    "localized_poa_candidate"));
                continue;
            }
            set.join(leftRoot, rightRoot);
        }
        disposition.accepted = true;
        disposition.classification = ConsistencyClass::CONSISTENT;
        disposition.reason = "copy-compatible homology evidence accepted";
        graph.evidenceAudit_.push_back(disposition);
    }

    std::map<std::size_t, std::vector<ResidueObservation>> components;
    for (std::size_t i = 0; i < residues.size(); ++i) {
        components[set.find(i)].push_back(evidence.residues.at(residues[i]));
    }
    std::map<ResidueId, std::string> residueToSite;
    for (auto &component : components) {
        AlignmentSite site;
        site.observations = component.second;
        std::sort(site.observations.begin(), site.observations.end(),
                  [](const ResidueObservation &a, const ResidueObservation &b) {
                      return a.id < b.id;
                  });
        site.siteId = stableId("site", {residuesKey(site.observations)});
        std::set<std::string> copyGroups;
        for (const ResidueObservation &observation : site.observations) {
            const auto assignment = options.copyAssignments.find(observation.id);
            if (assignment != options.copyAssignments.end() &&
                !assignment->second.copyGroup.empty()) {
                copyGroups.insert(assignment->second.copyGroup);
            }
        }
        if (copyGroups.size() > 1) {
            throw AlignmentGraphError(
                "one alignment site contains incompatible resolved copy groups");
        }
        if (!copyGroups.empty()) {
            site.copyGroup = *copyGroups.begin();
            site.copyResolved = true;
        }
        for (const ResidueObservation &observation : site.observations) {
            residueToSite[observation.id] = site.siteId;
        }
        graph.sites_[site.siteId] = site;
    }

    std::set<std::string> acceptedIds;
    for (const EvidenceDisposition &disposition : graph.evidenceAudit_) {
        if (disposition.accepted) {
            acceptedIds.insert(disposition.evidenceId);
        }
    }
    for (const HomologyEvidence &edge : edges) {
        if (!acceptedIds.count(edge.provenance.evidenceId)) {
            continue;
        }
        const std::string &siteId = residueToSite.at(edge.left);
        AlignmentSite &site = graph.sites_.at(siteId);
        site.directPairSupport.insert(edge.provenance.pairId);
        site.evidenceIds.push_back(edge.provenance.evidenceId);
    }
    for (auto &entry : graph.sites_) {
        AlignmentSite &site = entry.second;
        std::sort(site.evidenceIds.begin(), site.evidenceIds.end());
        site.evidenceIds.erase(std::unique(site.evidenceIds.begin(), site.evidenceIds.end()),
                               site.evidenceIds.end());
        site.consistency = classifySite(site, options);
    }

    // A direct aligned-gap statement contradicted by an accepted transitive
    // residue from the absent taxon is a sequence-level inconsistency. It is
    // retained in the audit and may be repaired only inside a bounded window.
    for (const AlignedAbsenceEvidence &absence : evidence.alignedAbsences) {
        EvidenceDisposition disposition;
        disposition.evidenceId = absence.provenance.evidenceId;
        const auto siteFound = residueToSite.find(absence.present);
        if (siteFound == residueToSite.end()) {
            disposition.classification = ConsistencyClass::MISSING_EVIDENCE;
            disposition.reason = "gap evidence references a filtered present residue";
            graph.evidenceAudit_.push_back(disposition);
            continue;
        }
        // A terminal gap has no pair of aligned flanks and cannot distinguish
        // biological absence from unaligned/missing sequence. Keep it missing;
        // only an internally bracketed gap is an explicit absence observation.
        if (absence.absentLeftFlank0 < 0 || absence.absentRightFlank0 < 0) {
            disposition.classification = ConsistencyClass::MISSING_EVIDENCE;
            disposition.reason = "aligned gap lacks two observed flanks";
            graph.evidenceAudit_.push_back(disposition);
            continue;
        }
        const AlignmentSite &site = graph.sites_.at(siteFound->second);
        std::vector<ResidueId> contradictory;
        contradictory.push_back(absence.present);
        for (const ResidueObservation &observation : site.observations) {
            if (observation.id.taxon == absence.absentTaxon) {
                contradictory.push_back(observation.id);
            }
        }
        if (contradictory.size() > 1) {
            disposition.classification = ConsistencyClass::RESIDUE_TRANSITIVITY;
            disposition.reason = "aligned absence contradicts a transitive residue";
            graph.evidenceAudit_.push_back(disposition);
            graph.conflicts_.push_back(makeConflict(
                ConsistencyClass::RESIDUE_TRANSITIVITY,
                {absence.provenance.evidenceId}, contradictory, true,
                    "localized_poa_candidate"));
        } else {
            AlignmentSite &mutableSite = graph.sites_.at(siteFound->second);
            mutableSite.alignedAbsentPaths.insert(absence.absentOccurrencePath);
            mutableSite.absenceEvidenceIds.push_back(absence.provenance.evidenceId);
            disposition.accepted = true;
            disposition.classification = ConsistencyClass::CONSISTENT;
            disposition.reason = "internally flanked aligned absence accepted";
            graph.evidenceAudit_.push_back(disposition);
        }
    }
    std::sort(graph.evidenceAudit_.begin(), graph.evidenceAudit_.end(),
              [](const EvidenceDisposition &a, const EvidenceDisposition &b) {
                  return a.evidenceId < b.evidenceId;
              });
    for (auto &entry : graph.sites_) {
        std::sort(entry.second.absenceEvidenceIds.begin(),
                  entry.second.absenceEvidenceIds.end());
        entry.second.absenceEvidenceIds.erase(
            std::unique(entry.second.absenceEvidenceIds.begin(),
                        entry.second.absenceEvidenceIds.end()),
            entry.second.absenceEvidenceIds.end());
    }

    std::map<std::string, std::vector<ResidueObservation>> paths;
    for (const auto &entry : evidence.residues) {
        paths[entry.first.occurrencePath].push_back(entry.second);
    }
    for (auto &entry : paths) {
        std::vector<ResidueObservation> &path = entry.second;
        std::sort(path.begin(), path.end(),
                  [](const ResidueObservation &a, const ResidueObservation &b) {
                      return a.id.forwardPosition0 < b.id.forwardPosition0;
                  });
        std::size_t begin = 0;
        while (begin < path.size()) {
            std::size_t end = begin + 1;
            while (end < path.size() && path[end].id.forwardPosition0 ==
                                            path[end - 1].id.forwardPosition0 + 1) {
                ++end;
            }
            PathSegment segment;
            segment.taxon = path[begin].id.taxon;
            segment.occurrencePath = path[begin].id.occurrencePath;
            segment.sequence = path[begin].id.sequence;
            segment.start0 = path[begin].id.forwardPosition0;
            segment.end0 = path[end - 1].id.forwardPosition0 + 1;
            segment.sourceSize = path[begin].sourceSize;
            segment.pathId = stableId("path", {entry.first,
                                                std::to_string(segment.start0),
                                                std::to_string(segment.end0)});
            for (std::size_t i = begin; i < end; ++i) {
                segment.siteIds.push_back(residueToSite.at(path[i].id));
                segment.sequenceText.push_back(path[i].base);
            }
            graph.pathSegments_.push_back(segment);
            begin = end;
        }
    }
    std::sort(graph.pathSegments_.begin(), graph.pathSegments_.end(),
              [](const PathSegment &a, const PathSegment &b) {
                  return a.pathId < b.pathId;
              });
    std::sort(graph.conflicts_.begin(), graph.conflicts_.end(),
              [](const GraphConflict &a, const GraphConflict &b) {
                  return a.conflictId < b.conflictId;
              });
    graph.validateInvariants();
    return graph;
}

void AlignmentSiteGraph::validateInvariants() const {
    std::set<ResidueId> seen;
    for (const auto &entry : sites_) {
        if (entry.first != entry.second.siteId) {
            throw AlignmentGraphError("site map key differs from stable site ID");
        }
        std::set<std::string> occurrencePaths;
        for (const ResidueObservation &observation : entry.second.observations) {
            if (observation.sourceSize < 0 || observation.id.forwardPosition0 < 0 ||
                observation.id.forwardPosition0 >= observation.sourceSize) {
                throw AlignmentGraphError("alignment-site residue is outside its source sequence");
            }
            if (!seen.insert(observation.id).second) {
                throw AlignmentGraphError("one source residue occurs in more than one alignment site");
            }
            if (!occurrencePaths.insert(observation.id.occurrencePath).second) {
                throw AlignmentGraphError("one alignment site contains two residues from an occurrence path");
            }
        }
    }
    for (const PathSegment &segment : pathSegments_) {
        if (segment.start0 < 0 || segment.end0 < segment.start0 ||
            segment.end0 > segment.sourceSize ||
            static_cast<int64_t>(segment.siteIds.size()) != segment.end0 - segment.start0 ||
            segment.sequenceText.size() != segment.siteIds.size()) {
            throw AlignmentGraphError("path segment coordinate or spelling invariant failed");
        }
        for (const std::string &siteId : segment.siteIds) {
            if (!sites_.count(siteId)) {
                throw AlignmentGraphError("path segment references an unknown alignment site");
            }
        }
    }
}

std::string AlignmentSiteGraph::coreHash() const {
    std::vector<std::string> fields;
    for (const auto &entry : sites_) {
        std::ostringstream site;
        site << entry.first << ':' << residuesKey(entry.second.observations);
        site << ":copy=" << (entry.second.copyResolved ? entry.second.copyGroup : ".");
        for (const std::string &path : entry.second.alignedAbsentPaths) {
            site << ":gap=" << path;
        }
        fields.push_back(site.str());
    }
    for (const PathSegment &segment : pathSegments_) {
        std::ostringstream path;
        path << segment.pathId << ':' << segment.occurrencePath << ':'
             << segment.taxon << ':' << segment.sequence << ':'
             << segment.start0 << ':' << segment.end0 << ':' << segment.sourceSize
             << ':' << segment.sequenceText;
        for (const std::string &site : segment.siteIds) {
            path << ':' << site;
        }
        fields.push_back(path.str());
    }
    return stableId("core", fields);
}

}  // namespace trio
}  // namespace anchorwave
