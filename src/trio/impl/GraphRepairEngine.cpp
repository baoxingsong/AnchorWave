#include "src/trio/impl/GraphRepairEngine.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

struct WindowPath {
    std::size_t segmentIndex = 0;
    std::size_t leftIndex = 0;
    std::size_t rightIndex = 0;
    std::vector<ResidueObservation> coreObservations;
};

struct ConflictWindow {
    std::string leftSite;
    std::string rightSite;
    std::set<std::string> oldCoreSites;
    std::map<std::string, WindowPath> paths;
    LocalizedPoaRequest request;
};

const ResidueObservation *observationForPath(const AlignmentSite &site,
                                             const std::string &path) {
    for (const ResidueObservation &observation : site.observations) {
        if (observation.id.occurrencePath == path) return &observation;
    }
    return nullptr;
}

bool isBoundary(const AlignmentSiteGraph &graph, const std::string &siteId) {
    const AlignmentSite &site = graph.sites().at(siteId);
    return site.consistency == ConsistencyClass::CONSISTENT && site.copyResolved;
}

std::map<ResidueId, std::string> residueSites(const AlignmentSiteGraph &graph) {
    std::map<ResidueId, std::string> result;
    for (const auto &site : graph.sites()) {
        for (const ResidueObservation &observation : site.second.observations)
            result[observation.id] = site.first;
    }
    return result;
}

bool findIndex(const PathSegment &path, const std::string &site,
               std::size_t &index) {
    const auto found = std::find(path.siteIds.begin(), path.siteIds.end(), site);
    if (found == path.siteIds.end()) return false;
    index = static_cast<std::size_t>(found - path.siteIds.begin());
    return true;
}

std::string selectSharedBoundary(
    const AlignmentSiteGraph &graph,
    const std::map<std::string, std::pair<const PathSegment *, std::size_t>> &seeds,
    bool left) {
    std::map<std::string, std::size_t> count;
    std::map<std::string, std::size_t> worstDistance;
    for (const auto &entry : seeds) {
        const PathSegment &path = *entry.second.first;
        const std::size_t seed = entry.second.second;
        if (left) {
            for (std::size_t i = 0; i < seed; ++i) {
                if (!isBoundary(graph, path.siteIds[i])) continue;
                ++count[path.siteIds[i]];
                worstDistance[path.siteIds[i]] = std::max(
                    worstDistance[path.siteIds[i]], seed - i);
            }
        } else {
            for (std::size_t i = seed + 1; i < path.siteIds.size(); ++i) {
                if (!isBoundary(graph, path.siteIds[i])) continue;
                ++count[path.siteIds[i]];
                worstDistance[path.siteIds[i]] = std::max(
                    worstDistance[path.siteIds[i]], i - seed);
            }
        }
    }
    std::string best;
    std::size_t bestDistance = std::numeric_limits<std::size_t>::max();
    for (const auto &entry : count) {
        if (entry.second != seeds.size()) continue;
        const std::size_t distance = worstDistance[entry.first];
        if (distance < bestDistance || (distance == bestDistance && entry.first < best)) {
            best = entry.first; bestDistance = distance;
        }
    }
    return best;
}

std::string outsideHash(const AlignmentSiteGraph &graph,
                        const std::set<std::string> &coreSites) {
    std::vector<std::string> fields;
    for (const auto &entry : graph.sites()) {
        if (coreSites.count(entry.first)) continue;
        std::ostringstream site;
        site << entry.first << ':' << entry.second.copyGroup << ':'
             << consistencyClassName(entry.second.consistency);
        for (const ResidueObservation &observation : entry.second.observations)
            site << ':' << observation.id.canonicalString() << '=' << observation.base;
        fields.push_back(site.str());
    }
    for (const PathSegment &path : graph.pathSegments()) {
        std::ostringstream value;
        value << path.pathId << ':' << path.sequenceText;
        for (const std::string &site : path.siteIds)
            if (!coreSites.count(site)) value << ':' << site;
        fields.push_back(value.str());
    }
    return stableId("outside", fields);
}

bool buildWindow(const AlignmentSiteGraph &graph, const GraphConflict &conflict,
                 std::size_t maximumCoreSites, ConflictWindow &window,
                 std::string &failure) {
    const std::map<ResidueId, std::string> residueToSite = residueSites(graph);
    std::map<std::string, std::set<std::string>> seedSitesByPath;
    for (const ResidueId &residue : conflict.residues) {
        const auto site = residueToSite.find(residue);
        if (site == residueToSite.end()) {
            failure = "conflict_residue_missing_from_graph"; return false;
        }
        seedSitesByPath[residue.occurrencePath].insert(site->second);
    }
    if (seedSitesByPath.size() < 2) {
        failure = "fewer_than_two_conflict_paths"; return false;
    }
    std::map<std::string, std::pair<const PathSegment *, std::size_t>> seedPositions;
    for (const auto &entry : seedSitesByPath) {
        const PathSegment *selected = nullptr;
        std::size_t selectedIndex = 0;
        for (const PathSegment &path : graph.pathSegments()) {
            if (path.occurrencePath != entry.first) continue;
            for (const std::string &seedSite : entry.second) {
                std::size_t index = 0;
                if (findIndex(path, seedSite, index)) {
                    if (selected != nullptr && selected != &path) {
                        failure = "conflict_spans_multiple_path_segments"; return false;
                    }
                    selected = &path;
                    selectedIndex = index;
                }
            }
        }
        if (selected == nullptr) {
            failure = "conflict_path_segment_not_found"; return false;
        }
        seedPositions[entry.first] = std::make_pair(selected, selectedIndex);
    }
    window.leftSite = selectSharedBoundary(graph, seedPositions, true);
    window.rightSite = selectSharedBoundary(graph, seedPositions, false);
    if (window.leftSite.empty() || window.rightSite.empty()) {
        failure = "no_shared_concordant_flanks"; return false;
    }

    for (std::size_t segmentIndex = 0; segmentIndex < graph.pathSegments().size();
         ++segmentIndex) {
        const PathSegment &path = graph.pathSegments()[segmentIndex];
        std::size_t left = 0, right = 0;
        if (!findIndex(path, window.leftSite, left) ||
            !findIndex(path, window.rightSite, right) || left >= right) continue;
        WindowPath windowPath;
        windowPath.segmentIndex = segmentIndex;
        windowPath.leftIndex = left; windowPath.rightIndex = right;
        std::string sequence;
        for (std::size_t i = left + 1; i < right; ++i) {
            const std::string &siteId = path.siteIds[i];
            window.oldCoreSites.insert(siteId);
            const ResidueObservation *observation =
                observationForPath(graph.sites().at(siteId), path.occurrencePath);
            if (observation == nullptr) {
                failure = "path_site_lacks_source_observation"; return false;
            }
            windowPath.coreObservations.push_back(*observation);
            sequence.push_back(observation->base);
        }
        window.paths[path.occurrencePath] = windowPath;
        window.request.sequences.push_back({path.occurrencePath, sequence});
    }
    for (const auto &entry : seedSitesByPath) {
        if (!window.paths.count(entry.first)) {
            failure = "conflict_path_does_not_traverse_both_flanks"; return false;
        }
    }
    if (window.oldCoreSites.empty() || window.oldCoreSites.size() > maximumCoreSites) {
        failure = window.oldCoreSites.empty() ? "empty_mutable_core" : "core_too_large";
        return false;
    }
    // Every observation in the mutable core must belong to a path included in
    // this repair. Otherwise a structural/missing boundary was crossed.
    for (const std::string &siteId : window.oldCoreSites) {
        for (const ResidueObservation &observation : graph.sites().at(siteId).observations) {
            if (!window.paths.count(observation.id.occurrencePath)) {
                failure = "core_contains_path_without_shared_flanks"; return false;
            }
        }
    }

    // Use the current graph projection as a baseline only when it has a unique
    // complete order. The union of extant path orders is topologically sorted.
    std::map<std::string, std::set<std::string>> edges;
    std::map<std::string, std::size_t> indegree;
    for (const std::string &site : window.oldCoreSites) indegree[site] = 0;
    for (const auto &entry : window.paths) {
        const PathSegment &path = graph.pathSegments()[entry.second.segmentIndex];
        for (std::size_t i = entry.second.leftIndex + 2; i < entry.second.rightIndex; ++i) {
            const std::string &a = path.siteIds[i - 1];
            const std::string &b = path.siteIds[i];
            if (a != b && edges[a].insert(b).second) ++indegree[b];
        }
    }
    std::set<std::string> ready;
    for (const auto &entry : indegree) if (entry.second == 0) ready.insert(entry.first);
    std::vector<std::string> order;
    bool unique = true;
    while (!ready.empty()) {
        if (ready.size() != 1) unique = false;
        const std::string current = *ready.begin(); ready.erase(ready.begin());
        order.push_back(current);
        for (const std::string &next : edges[current])
            if (--indegree[next] == 0) ready.insert(next);
    }
    if (unique && order.size() == window.oldCoreSites.size()) {
        bool complete = true;
        for (const auto &entry : window.paths) {
            std::string row;
            for (const std::string &siteId : order) {
                const AlignmentSite &site = graph.sites().at(siteId);
                const ResidueObservation *observation =
                    observationForPath(site, entry.first);
                if (observation != nullptr) row.push_back(observation->base);
                else if (site.alignedAbsentPaths.count(entry.first)) row.push_back('-');
                else { complete = false; break; }
            }
            if (!complete) break;
            window.request.baselineRows[entry.first] = row;
        }
        if (!complete) window.request.baselineRows.clear();
    }
    window.request.conflictId = conflict.conflictId;
    window.request.immutableLeftSite = window.leftSite;
    window.request.immutableRightSite = window.rightSite;
    return true;
}

void applyPatch(std::map<std::string, AlignmentSite> &sites,
                std::vector<PathSegment> &pathSegments,
                const ConflictWindow &window,
                const LocalizedPoaPatch &patch,
                std::set<std::string> &newCoreSites) {
    std::map<std::string, std::size_t> rowIndex;
    for (std::size_t i = 0; i < patch.pathIds.size(); ++i) rowIndex[patch.pathIds[i]] = i;
    std::map<std::string, std::size_t> consumed;
    std::map<std::string, std::vector<std::string>> replacementPathSites;
    std::set<std::string> copyGroups;
    for (const std::string &siteId : window.oldCoreSites) {
        const AlignmentSite &site = sites.at(siteId);
        if (site.copyResolved) copyGroups.insert(site.copyGroup);
    }
    if (copyGroups.size() > 1) {
        throw AlignmentGraphError("localized POA core crosses resolved copy groups");
    }
    const std::size_t width = patch.alignedRows.empty() ? 0 : patch.alignedRows[0].size();
    for (std::size_t column = 0; column < width; ++column) {
        AlignmentSite site;
        site.consistency = ConsistencyClass::LOCAL_POA_REPAIRED;
        site.copyResolved = !copyGroups.empty();
        if (site.copyResolved) site.copyGroup = *copyGroups.begin();
        site.evidenceIds.push_back(patch.repairId);
        std::vector<std::string> residueKeys;
        for (const std::string &path : patch.pathIds) {
            const char base = patch.alignedRows[rowIndex[path]][column];
            if (base == '-') {
                site.alignedAbsentPaths.insert(path);
                site.absenceEvidenceIds.push_back(patch.repairId);
                continue;
            }
            const WindowPath &windowPath = window.paths.at(path);
            const std::size_t offset = consumed[path]++;
            if (offset >= windowPath.coreObservations.size())
                throw AlignmentGraphError("POA patch consumes too many source residues");
            const ResidueObservation &observation = windowPath.coreObservations[offset];
            if (observation.base != base)
                throw AlignmentGraphError("POA patch changes an extant source base");
            site.observations.push_back(observation);
            residueKeys.push_back(observation.id.canonicalString());
        }
        if (site.observations.empty())
            throw AlignmentGraphError("POA patch created an all-gap alignment site");
        std::sort(site.observations.begin(), site.observations.end(),
                  [](const ResidueObservation &a, const ResidueObservation &b) {
                      return a.id < b.id;
                  });
        std::sort(residueKeys.begin(), residueKeys.end());
        std::vector<std::string> idFields = {patch.repairId, std::to_string(column)};
        idFields.insert(idFields.end(), residueKeys.begin(), residueKeys.end());
        site.siteId = stableId("repair_site", idFields);
        newCoreSites.insert(site.siteId);
        sites[site.siteId] = site;
        for (const ResidueObservation &observation : site.observations)
            replacementPathSites[observation.id.occurrencePath].push_back(site.siteId);
    }
    for (const auto &entry : window.paths) {
        if (consumed[entry.first] != entry.second.coreObservations.size())
            throw AlignmentGraphError("POA patch does not consume every source residue");
    }
    for (const auto &entry : window.paths) {
        PathSegment &path = pathSegments[entry.second.segmentIndex];
        std::vector<std::string> updated;
        updated.insert(updated.end(), path.siteIds.begin(),
                       path.siteIds.begin() + static_cast<std::ptrdiff_t>(entry.second.leftIndex + 1));
        const std::vector<std::string> &replacement = replacementPathSites[entry.first];
        updated.insert(updated.end(), replacement.begin(), replacement.end());
        updated.insert(updated.end(),
                       path.siteIds.begin() + static_cast<std::ptrdiff_t>(entry.second.rightIndex),
                       path.siteIds.end());
        path.siteIds.swap(updated);
    }
    for (const std::string &siteId : window.oldCoreSites) sites.erase(siteId);
}

}  // namespace

GraphRepairResult GraphRepairEngine::repairEligibleConflicts(
    const AlignmentSiteGraph &input, const GraphRepairOptions &options) {
    GraphRepairResult result;
    result.graph = input;
    std::set<std::string> alreadyRepairedResidues;
    for (GraphConflict &conflict : result.graph.conflicts_) {
        if (!conflict.poaEligible) continue;
        bool overlaps = false;
        for (const ResidueId &residue : conflict.residues)
            overlaps |= alreadyRepairedResidues.count(residue.canonicalString());
        GraphRepairAudit audit;
        audit.conflictId = conflict.conflictId;
        if (overlaps) {
            audit.disposition = "overlaps_previously_repaired_core";
            conflict.disposition = audit.disposition;
            result.audit.push_back(audit);
            continue;
        }
        ConflictWindow window;
        std::string failure;
        if (!buildWindow(result.graph, conflict, options.maximumCoreSites,
                         window, failure)) {
            audit.disposition = failure;
            conflict.disposition = failure;
            result.audit.push_back(audit);
            continue;
        }
        audit.immutableLeftSite = window.leftSite;
        audit.immutableRightSite = window.rightSite;
        audit.oldCoreSites.assign(window.oldCoreSites.begin(), window.oldCoreSites.end());
        audit.outsideHashBefore = outsideHash(result.graph, window.oldCoreSites);
        const LocalizedPoaPatch patch = LocalizedPoaRepair::repair(window.request, options.poa);
        audit.repairId = patch.repairId;
        audit.baselineScore = patch.baselineScore;
        audit.repairedScore = patch.repairedScore;
        audit.scoreDelta = patch.scoreDelta;
        if (!patch.accepted) {
            audit.disposition = patch.disposition;
            conflict.disposition = patch.disposition;
            result.audit.push_back(audit);
            continue;
        }
        std::set<std::string> newCore;
        applyPatch(result.graph.sites_, result.graph.pathSegments_, window, patch, newCore);
        result.graph.validateInvariants();
        audit.newCoreSites.assign(newCore.begin(), newCore.end());
        audit.outsideHashAfter = outsideHash(result.graph, newCore);
        if (audit.outsideHashBefore != audit.outsideHashAfter) {
            throw AlignmentGraphError("localized POA changed graph state outside mutable core");
        }
        audit.disposition = "localized_poa_applied";
        conflict.disposition = audit.disposition;
        for (const ResidueId &residue : conflict.residues)
            alreadyRepairedResidues.insert(residue.canonicalString());
        result.audit.push_back(audit);
    }
    return result;
}

}  // namespace trio
}  // namespace anchorwave
