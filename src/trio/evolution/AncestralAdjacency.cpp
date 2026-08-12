#include "src/trio/evolution/AncestralAdjacency.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <set>
#include <sstream>
#include <utility>

namespace anchorwave {
namespace trio {
namespace {

const char *const kCandidateOnlyScope = "CANDIDATE_BLOCK_ADJACENCIES_ONLY";

struct CandidateKey {
    OrientedBlockEndpoint first;
    OrientedBlockEndpoint second;

    bool operator<(const CandidateKey &other) const {
        if (first < other.first) return true;
        if (other.first < first) return false;
        return second < other.second;
    }
};

struct TaxonBlockKey {
    std::string taxon;
    std::string blockId;

    bool operator<(const TaxonBlockKey &other) const {
        if (taxon != other.taxon) return taxon < other.taxon;
        return blockId < other.blockId;
    }
};

struct TaxonPathKey {
    std::string taxon;
    std::string path;

    bool operator<(const TaxonPathKey &other) const {
        if (taxon != other.taxon) return taxon < other.taxon;
        return path < other.path;
    }
};

CandidateKey canonicalCandidate(OrientedBlockEndpoint left,
                                OrientedBlockEndpoint right) {
    if (right < left) std::swap(left, right);
    CandidateKey result;
    result.first = left;
    result.second = right;
    return result;
}

OrientedBlockEndpoint genomicLeftEnd(const ExtantBlockProjection &projection) {
    OrientedBlockEndpoint result;
    result.blockId = projection.blockId;
    result.end = projection.orientation == BlockOrientation::FORWARD
                     ? BlockEnd::HEAD
                     : BlockEnd::TAIL;
    return result;
}

OrientedBlockEndpoint genomicRightEnd(const ExtantBlockProjection &projection) {
    OrientedBlockEndpoint result;
    result.blockId = projection.blockId;
    result.end = projection.orientation == BlockOrientation::FORWARD
                     ? BlockEnd::TAIL
                     : BlockEnd::HEAD;
    return result;
}

std::string adjacencyId(const CandidateKey &candidate) {
    return stableId("adjacency",
                    {candidate.first.blockId, blockEndName(candidate.first.end),
                     candidate.second.blockId, blockEndName(candidate.second.end)});
}

std::string join(const std::vector<std::string> &values, const char delimiter) {
    std::ostringstream output;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) output << delimiter;
        output << values[i];
    }
    return output.str();
}

std::string cleanTsv(std::string value) {
    for (char &character : value) {
        if (character == '\t' || character == '\n' || character == '\r') {
            character = ' ';
        }
    }
    return value;
}

BlockProjectionIssue makeIssue(BlockProjectionIssueKind kind,
                               const std::string &taxon,
                               const std::string &path,
                               std::vector<std::string> blockIds,
                               std::vector<std::string> projectionIds,
                               const std::string &detail) {
    std::sort(blockIds.begin(), blockIds.end());
    blockIds.erase(std::unique(blockIds.begin(), blockIds.end()), blockIds.end());
    std::sort(projectionIds.begin(), projectionIds.end());
    projectionIds.erase(std::unique(projectionIds.begin(), projectionIds.end()),
                        projectionIds.end());
    std::vector<std::string> fields = {blockProjectionIssueKindName(kind), taxon, path};
    fields.insert(fields.end(), blockIds.begin(), blockIds.end());
    fields.insert(fields.end(), projectionIds.begin(), projectionIds.end());
    BlockProjectionIssue issue;
    issue.issueId = stableId("adjacency_issue", fields);
    issue.kind = kind;
    issue.taxon = taxon;
    issue.occurrencePath = path;
    issue.blockIds = blockIds;
    issue.projectionIds = projectionIds;
    issue.detail = detail;
    return issue;
}

bool issueLess(const BlockProjectionIssue &left,
               const BlockProjectionIssue &right) {
    return left.issueId < right.issueId;
}

void sortAndDeduplicateIssues(std::vector<BlockProjectionIssue> &issues) {
    std::sort(issues.begin(), issues.end(), issueLess);
    issues.erase(std::unique(issues.begin(), issues.end(),
                             [](const BlockProjectionIssue &left,
                                const BlockProjectionIssue &right) {
                                 return left.issueId == right.issueId;
                             }),
                 issues.end());
}

bool projectionLess(const ExtantBlockProjection &left,
                    const ExtantBlockProjection &right) {
    if (left.taxon != right.taxon) return left.taxon < right.taxon;
    if (left.occurrencePath != right.occurrencePath) {
        return left.occurrencePath < right.occurrencePath;
    }
    if (left.start0 != right.start0) return left.start0 < right.start0;
    if (left.end0 != right.end0) return left.end0 < right.end0;
    if (left.blockId != right.blockId) return left.blockId < right.blockId;
    return left.projectionId < right.projectionId;
}

bool matchesAt(const std::vector<std::string> &path,
               const std::vector<std::string> &block,
               std::size_t start,
               bool reverse) {
    if (start + block.size() > path.size()) return false;
    for (std::size_t offset = 0; offset < block.size(); ++offset) {
        const std::size_t blockIndex = reverse ? block.size() - offset - 1 : offset;
        if (path[start + offset] != block[blockIndex]) return false;
    }
    return true;
}

std::pair<std::string, std::string> unorderedBlockPair(
    const AncestralAdjacencyCall &call) {
    if (call.endpoint2.blockId < call.endpoint1.blockId) {
        return std::make_pair(call.endpoint2.blockId, call.endpoint1.blockId);
    }
    return std::make_pair(call.endpoint1.blockId, call.endpoint2.blockId);
}

std::vector<std::string> leavesWithObservation(
    const AncestralAdjacencyCall &call, PresenceObservation requested) {
    std::vector<std::string> result;
    for (const auto &entry : call.leafObservations) {
        if (entry.second == requested) result.push_back(entry.first);
    }
    return result;
}

}  // namespace

std::string blockOrientationName(BlockOrientation orientation) {
    switch (orientation) {
        case BlockOrientation::FORWARD: return "+";
        case BlockOrientation::REVERSE: return "-";
        case BlockOrientation::UNKNOWN: return "?";
    }
    return "?";
}

std::string blockEndName(BlockEnd end) {
    switch (end) {
        case BlockEnd::HEAD: return "HEAD";
        case BlockEnd::TAIL: return "TAIL";
    }
    return "UNKNOWN";
}

bool OrientedBlockEndpoint::operator<(const OrientedBlockEndpoint &other) const {
    if (blockId != other.blockId) return blockId < other.blockId;
    return static_cast<int>(end) < static_cast<int>(other.end);
}

bool OrientedBlockEndpoint::operator==(const OrientedBlockEndpoint &other) const {
    return blockId == other.blockId && end == other.end;
}

std::string blockProjectionIssueKindName(BlockProjectionIssueKind kind) {
    switch (kind) {
        case BlockProjectionIssueKind::UNPROJECTED_BLOCK:
            return "UNPROJECTED_BLOCK";
        case BlockProjectionIssueKind::UNKNOWN_ORIENTATION:
            return "UNKNOWN_ORIENTATION";
        case BlockProjectionIssueKind::DUPLICATE_BLOCK_COPY:
            return "DUPLICATE_BLOCK_COPY";
        case BlockProjectionIssueKind::OVERLAPPING_PROJECTIONS:
            return "OVERLAPPING_PROJECTIONS";
        case BlockProjectionIssueKind::SELF_ADJACENCY:
            return "SELF_ADJACENCY";
    }
    return "UNKNOWN";
}

void BlockProjectionSet::validateAgainst(const AlignmentSiteGraph &graph) const {
    if (immutableCoreHash != graph.coreHash()) {
        throw AncestralAdjacencyError(
            "block projections refer to a different alignment graph core");
    }
    std::map<std::pair<std::string, std::string>, int64_t> graphPaths;
    for (const PathSegment &path : graph.pathSegments()) {
        const std::pair<std::string, std::string> key(path.taxon,
                                                      path.occurrencePath);
        const auto existing = graphPaths.find(key);
        if (existing != graphPaths.end() && existing->second != path.sourceSize) {
            throw AncestralAdjacencyError(
                "graph path segments disagree on source size: " +
                path.occurrencePath);
        }
        graphPaths[key] = path.sourceSize;
    }
    std::set<std::string> projectionIds;
    for (const ExtantBlockProjection &projection : projections) {
        if (projection.projectionId.empty() || projection.taxon.empty() ||
            projection.occurrencePath.empty() || projection.blockId.empty()) {
            throw AncestralAdjacencyError(
                "a block projection has an empty required identifier");
        }
        if (projection.start0 < 0 || projection.end0 <= projection.start0) {
            throw AncestralAdjacencyError(
                "a block projection has an invalid half-open interval");
        }
        if (!projectionIds.insert(projection.projectionId).second) {
            throw AncestralAdjacencyError("duplicate block projection ID: " +
                                          projection.projectionId);
        }
        const auto graphPath = graphPaths.find(std::make_pair(
            projection.taxon, projection.occurrencePath));
        if (graphPath == graphPaths.end()) {
            throw AncestralAdjacencyError(
                "block projection names a taxon/path absent from the graph: " +
                projection.taxon + "/" + projection.occurrencePath);
        }
        if (projection.end0 > graphPath->second) {
            throw AncestralAdjacencyError(
                "block projection exceeds graph source size: " +
                projection.projectionId);
        }
        const std::string prefix = projection.taxon + "|";
        if (projection.occurrencePath.compare(0, prefix.size(), prefix) != 0) {
            throw AncestralAdjacencyError(
                "block projection taxon and occurrence path disagree: " +
                projection.projectionId);
        }
    }
}

BlockProjectionSet AncestralAdjacencyInference::projectBlocks(
    const AlignmentSiteGraph &graph,
    const std::vector<StructuralBlockDefinition> &definitions,
    const BlockProjectionOptions &options) {
    const std::string frozenHash = graph.coreHash();
    std::vector<StructuralBlockDefinition> orderedDefinitions = definitions;
    std::sort(orderedDefinitions.begin(), orderedDefinitions.end(),
              [](const StructuralBlockDefinition &left,
                 const StructuralBlockDefinition &right) {
                  return left.blockId < right.blockId;
              });
    std::set<std::string> blockIds;
    std::map<std::string, std::string> siteOwner;
    for (const StructuralBlockDefinition &definition : orderedDefinitions) {
        if (definition.blockId.empty() || definition.siteIds.empty()) {
            throw AncestralAdjacencyError(
                "every structural block needs an ID and at least one site");
        }
        if (!blockIds.insert(definition.blockId).second) {
            throw AncestralAdjacencyError("duplicate structural block ID: " +
                                          definition.blockId);
        }
        std::set<std::string> copyGroups;
        std::set<std::string> localSites;
        for (const std::string &siteId : definition.siteIds) {
            const auto site = graph.sites().find(siteId);
            if (site == graph.sites().end()) {
                throw AncestralAdjacencyError(
                    "structural block references an unknown graph site: " + siteId);
            }
            if (!localSites.insert(siteId).second) {
                throw AncestralAdjacencyError(
                    "structural block repeats graph site: " + siteId);
            }
            const auto owner = siteOwner.find(siteId);
            if (owner != siteOwner.end()) {
                throw AncestralAdjacencyError(
                    "graph site belongs to more than one structural block: " + siteId);
            }
            siteOwner[siteId] = definition.blockId;
            if (options.requireCopyResolved &&
                (!site->second.copyResolved || site->second.copyGroup.empty())) {
                throw AncestralAdjacencyError(
                    "structural block contains a copy-unresolved graph site: " + siteId);
            }
            if (site->second.copyResolved && !site->second.copyGroup.empty()) {
                copyGroups.insert(site->second.copyGroup);
            }
        }
        if (options.requireSingleCopyGroup && copyGroups.size() > 1) {
            throw AncestralAdjacencyError(
                "structural block spans more than one resolved copy group: " +
                definition.blockId);
        }
    }

    std::vector<PathSegment> paths = graph.pathSegments();
    std::sort(paths.begin(), paths.end(),
              [](const PathSegment &left, const PathSegment &right) {
                  return left.pathId < right.pathId;
              });
    BlockProjectionSet result;
    result.immutableCoreHash = frozenHash;
    std::map<std::string, std::size_t> projectionCounts;
    for (const StructuralBlockDefinition &definition : orderedDefinitions) {
        for (const PathSegment &path : paths) {
            if (definition.siteIds.size() > path.siteIds.size()) continue;
            for (std::size_t start = 0;
                 start + definition.siteIds.size() <= path.siteIds.size(); ++start) {
                const bool forward = matchesAt(path.siteIds, definition.siteIds,
                                               start, false);
                const bool reverse = matchesAt(path.siteIds, definition.siteIds,
                                               start, true);
                if (!forward && !reverse) continue;
                ExtantBlockProjection projection;
                projection.taxon = path.taxon;
                projection.occurrencePath = path.occurrencePath;
                projection.blockId = definition.blockId;
                projection.start0 = path.start0 + static_cast<int64_t>(start);
                projection.end0 = projection.start0 +
                                  static_cast<int64_t>(definition.siteIds.size());
                projection.orientation = forward == reverse
                                             ? BlockOrientation::UNKNOWN
                                             : (forward ? BlockOrientation::FORWARD
                                                        : BlockOrientation::REVERSE);
                projection.projectionId = stableId(
                    "block_projection",
                    {frozenHash, path.pathId, definition.blockId,
                     std::to_string(projection.start0),
                     std::to_string(projection.end0),
                     blockOrientationName(projection.orientation)});
                result.projections.push_back(projection);
                ++projectionCounts[definition.blockId];
                if (projection.orientation == BlockOrientation::UNKNOWN) {
                    result.issues.push_back(makeIssue(
                        BlockProjectionIssueKind::UNKNOWN_ORIENTATION,
                        projection.taxon, projection.occurrencePath,
                        {projection.blockId}, {projection.projectionId},
                        "block direction cannot be resolved from this projection"));
                }
            }
        }
        if (projectionCounts[definition.blockId] == 0) {
            result.issues.push_back(makeIssue(
                BlockProjectionIssueKind::UNPROJECTED_BLOCK, "", "",
                {definition.blockId}, {},
                "structural block has no exact contiguous graph-path projection"));
        }
    }
    std::sort(result.projections.begin(), result.projections.end(), projectionLess);
    sortAndDeduplicateIssues(result.issues);
    result.validateAgainst(graph);
    if (graph.coreHash() != frozenHash) {
        throw AncestralAdjacencyError(
            "alignment graph changed during read-only block projection");
    }
    return result;
}

AncestralAdjacencyReport AncestralAdjacencyInference::infer(
    const AlignmentSiteGraph &graph,
    const Phylogeny &tree,
    const BlockProjectionSet &projectionSet,
    const AncestralAdjacencyConfig &config) {
    if (config.targetNode.empty()) {
        throw AncestralAdjacencyError("target ancestor node must be named");
    }
    tree.internalNodeId(config.targetNode);
    projectionSet.validateAgainst(graph);
    const std::string frozenHash = graph.coreHash();
    const std::vector<std::string> leaves = tree.leafNames();
    const std::set<std::string> leafSet(leaves.begin(), leaves.end());

    AncestralAdjacencyReport report;
    report.immutableCoreHash = frozenHash;
    report.targetNode = config.targetNode;
    report.inferenceScope = kCandidateOnlyScope;
    report.issues = projectionSet.issues;

    std::vector<ExtantBlockProjection> projections = projectionSet.projections;
    std::sort(projections.begin(), projections.end(), projectionLess);
    std::map<TaxonBlockKey, std::vector<std::size_t>> byTaxonBlock;
    std::map<TaxonPathKey, std::vector<std::size_t>> byTaxonPath;
    for (std::size_t i = 0; i < projections.size(); ++i) {
        const ExtantBlockProjection &projection = projections[i];
        if (!leafSet.count(projection.taxon)) {
            throw AncestralAdjacencyError(
                "block projection taxon is absent from the species tree: " +
                projection.taxon);
        }
        byTaxonBlock[{projection.taxon, projection.blockId}].push_back(i);
        byTaxonPath[{projection.taxon, projection.occurrencePath}].push_back(i);
    }

    std::set<std::size_t> unusable;
    for (const auto &entry : byTaxonBlock) {
        if (entry.second.size() <= 1) continue;
        std::vector<std::string> ids;
        for (std::size_t index : entry.second) {
            unusable.insert(index);
            ids.push_back(projections[index].projectionId);
        }
        report.issues.push_back(makeIssue(
            BlockProjectionIssueKind::DUPLICATE_BLOCK_COPY,
            entry.first.taxon, "", {entry.first.blockId}, ids,
            "more than one projection of a copy-resolved block makes adjacency state missing"));
    }
    for (std::size_t i = 0; i < projections.size(); ++i) {
        if (projections[i].orientation == BlockOrientation::UNKNOWN) {
            unusable.insert(i);
            report.issues.push_back(makeIssue(
                BlockProjectionIssueKind::UNKNOWN_ORIENTATION,
                projections[i].taxon, projections[i].occurrencePath,
                {projections[i].blockId}, {projections[i].projectionId},
                "unoriented projection cannot provide a HEAD/TAIL adjacency"));
        }
    }
    for (auto &entry : byTaxonPath) {
        std::vector<std::size_t> &indices = entry.second;
        std::sort(indices.begin(), indices.end(),
                  [&](std::size_t left, std::size_t right) {
                      return projectionLess(projections[left], projections[right]);
                  });
        for (std::size_t i = 0; i < indices.size(); ++i) {
            for (std::size_t j = i + 1; j < indices.size(); ++j) {
                const ExtantBlockProjection &left = projections[indices[i]];
                const ExtantBlockProjection &right = projections[indices[j]];
                if (right.start0 >= left.end0) break;
                unusable.insert(indices[i]);
                unusable.insert(indices[j]);
                report.issues.push_back(makeIssue(
                    BlockProjectionIssueKind::OVERLAPPING_PROJECTIONS,
                    entry.first.taxon, entry.first.path,
                    {left.blockId, right.blockId},
                    {left.projectionId, right.projectionId},
                    "overlapping structural-block projections cannot define a path adjacency"));
            }
        }
    }

    std::map<CandidateKey, std::set<std::string>> observedInLeaves;
    for (const auto &entry : byTaxonPath) {
        const std::vector<std::size_t> &indices = entry.second;
        for (std::size_t i = 1; i < indices.size(); ++i) {
            const std::size_t leftIndex = indices[i - 1];
            const std::size_t rightIndex = indices[i];
            if (unusable.count(leftIndex) || unusable.count(rightIndex)) continue;
            const ExtantBlockProjection &left = projections[leftIndex];
            const ExtantBlockProjection &right = projections[rightIndex];
            if (left.blockId == right.blockId) {
                report.issues.push_back(makeIssue(
                    BlockProjectionIssueKind::SELF_ADJACENCY,
                    entry.first.taxon, entry.first.path, {left.blockId},
                    {left.projectionId, right.projectionId},
                    "a repeated block identity cannot define a copy-resolved adjacency"));
                continue;
            }
            observedInLeaves[canonicalCandidate(genomicRightEnd(left),
                                                genomicLeftEnd(right))]
                .insert(entry.first.taxon);
        }
    }

    for (const auto &candidateEntry : observedInLeaves) {
        const CandidateKey &candidate = candidateEntry.first;
        AncestralAdjacencyCall call;
        call.adjacencyId = adjacencyId(candidate);
        call.endpoint1 = candidate.first;
        call.endpoint2 = candidate.second;
        for (const std::string &leaf : leaves) {
            const TaxonBlockKey firstKey = {leaf, candidate.first.blockId};
            const TaxonBlockKey secondKey = {leaf, candidate.second.blockId};
            const auto first = byTaxonBlock.find(firstKey);
            const auto second = byTaxonBlock.find(secondKey);
            bool callable = first != byTaxonBlock.end() && first->second.size() == 1 &&
                            second != byTaxonBlock.end() && second->second.size() == 1;
            if (callable) {
                callable = !unusable.count(first->second.front()) &&
                           !unusable.count(second->second.front());
            }
            if (!callable) {
                call.leafObservations[leaf] = PresenceObservation::MISSING;
            } else if (candidateEntry.second.count(leaf)) {
                call.leafObservations[leaf] = PresenceObservation::PRESENT;
            } else {
                call.leafObservations[leaf] = PresenceObservation::ABSENT;
            }
        }
        call.presence = tree.inferPresence(config.targetNode,
                                           call.leafObservations,
                                           config.adjacencyModel);
        report.calls.push_back(call);
    }
    std::sort(report.calls.begin(), report.calls.end(),
              [](const AncestralAdjacencyCall &left,
                 const AncestralAdjacencyCall &right) {
                  return left.adjacencyId < right.adjacencyId;
              });

    std::map<std::pair<std::string, std::string>, std::vector<std::size_t>>
        callsByBlockPair;
    for (std::size_t i = 0; i < report.calls.size(); ++i) {
        callsByBlockPair[unorderedBlockPair(report.calls[i])].push_back(i);
    }
    for (const auto &entry : callsByBlockPair) {
        if (entry.second.size() <= 1) continue;
        for (std::size_t index : entry.second) {
            report.calls[index].orientationConflict = true;
        }
    }

    std::map<OrientedBlockEndpoint, std::vector<std::size_t>> activeByEndpoint;
    for (std::size_t i = 0; i < report.calls.size(); ++i) {
        const AncestralAdjacencyCall &call = report.calls[i];
        const bool active = (!call.presence.isAmbiguous() &&
                             call.presence.selectedState == 1) ||
                            (call.presence.isAmbiguous() &&
                             config.ambiguousCallsCanConflict);
        if (!active) continue;
        activeByEndpoint[call.endpoint1].push_back(i);
        activeByEndpoint[call.endpoint2].push_back(i);
    }
    for (const auto &entry : activeByEndpoint) {
        if (entry.second.size() <= 1) continue;
        for (std::size_t index : entry.second) {
            report.calls[index].endpointConflict = true;
        }
    }

    for (AncestralAdjacencyCall &call : report.calls) {
        call.supportedWithoutLocalConflict =
            !call.presence.isAmbiguous() && call.presence.selectedState == 1 &&
            !call.orientationConflict && !call.endpointConflict;
        if (call.presence.isAmbiguous()) {
            call.decision = call.endpointConflict
                                ? "PRESENCE_AMBIGUOUS_ENDPOINT_CONFLICT"
                                : "PRESENCE_AMBIGUOUS";
        } else if (call.presence.selectedState == 0) {
            call.decision = "ANCESTRAL_ADJACENCY_ABSENT";
        } else if (call.endpointConflict) {
            call.decision = "SUPPORTED_ENDPOINT_CONFLICT";
        } else if (call.orientationConflict) {
            call.decision = "SUPPORTED_ORIENTATION_CONFLICT";
        } else {
            call.decision = "SUPPORTED_CANDIDATE_ADJACENCY";
        }
    }

    sortAndDeduplicateIssues(report.issues);
    report.validateAgainst(graph);
    if (graph.coreHash() != frozenHash) {
        throw AncestralAdjacencyError(
            "alignment graph changed during read-only adjacency inference");
    }
    return report;
}

void AncestralAdjacencyReport::validateAgainst(
    const AlignmentSiteGraph &graph) const {
    if (immutableCoreHash != graph.coreHash()) {
        throw AncestralAdjacencyError(
            "ancestral adjacency report refers to a different alignment graph core");
    }
    if (targetNode.empty() || inferenceScope != kCandidateOnlyScope) {
        throw AncestralAdjacencyError("invalid ancestral adjacency report metadata");
    }
    std::set<std::string> ids;
    for (const AncestralAdjacencyCall &call : calls) {
        if (call.endpoint2 < call.endpoint1 ||
            call.endpoint1.blockId.empty() || call.endpoint2.blockId.empty() ||
            call.endpoint1.blockId == call.endpoint2.blockId) {
            throw AncestralAdjacencyError(
                "ancestral adjacency has invalid or non-canonical endpoints");
        }
        const CandidateKey candidate = {call.endpoint1, call.endpoint2};
        if (call.adjacencyId != adjacencyId(candidate) ||
            !ids.insert(call.adjacencyId).second) {
            throw AncestralAdjacencyError(
                "ancestral adjacency has a duplicate or unstable ID");
        }
        const bool expectedSupported =
            !call.presence.isAmbiguous() && call.presence.selectedState == 1 &&
            !call.orientationConflict && !call.endpointConflict;
        if (call.supportedWithoutLocalConflict != expectedSupported) {
            throw AncestralAdjacencyError(
                "ancestral adjacency support flag is inconsistent with its call");
        }
    }
}

void AncestralAdjacencyWriters::writeProjections(
    const BlockProjectionSet &projectionSet, std::ostream &output) {
    output << "projection_id\ttaxon\toccurrence_path\tblock_id\tstart0\tend0"
              "\torientation\tcore_hash\n";
    std::vector<ExtantBlockProjection> projections = projectionSet.projections;
    std::sort(projections.begin(), projections.end(), projectionLess);
    for (const ExtantBlockProjection &projection : projections) {
        output << projection.projectionId << '\t' << projection.taxon << '\t'
               << projection.occurrencePath << '\t' << projection.blockId << '\t'
               << projection.start0 << '\t' << projection.end0 << '\t'
               << blockOrientationName(projection.orientation) << '\t'
               << projectionSet.immutableCoreHash << '\n';
    }
}

void AncestralAdjacencyWriters::writeCalls(
    const AncestralAdjacencyReport &report, std::ostream &output) {
    output << "adjacency_id\ttarget_node\tendpoint1_block\tendpoint1_end"
              "\tendpoint2_block\tendpoint2_end\toptimal_states\tselected_state"
              "\tbest_cost\tconfidence_margin\tpresent_leaves\tabsent_leaves"
              "\tmissing_leaves\tpresence_ambiguous\torientation_conflict\tendpoint_conflict"
              "\tsupported_without_local_conflict\tdecision\tinference_scope"
              "\tchromosome_reconstruction\tcore_hash\n";
    for (const AncestralAdjacencyCall &call : report.calls) {
        std::vector<std::string> states;
        for (int state : call.presence.optimalStates) {
            states.push_back(std::to_string(state));
        }
        output << call.adjacencyId << '\t' << report.targetNode << '\t'
               << call.endpoint1.blockId << '\t' << blockEndName(call.endpoint1.end)
               << '\t' << call.endpoint2.blockId << '\t'
               << blockEndName(call.endpoint2.end) << '\t' << join(states, ',')
               << '\t'
               << (call.presence.isAmbiguous()
                       ? "."
                       : std::to_string(call.presence.selectedState))
               << '\t'
               << call.presence.bestCost << '\t' << call.presence.confidenceMargin
               << '\t' << join(leavesWithObservation(call,
                                    PresenceObservation::PRESENT), ',')
               << '\t' << join(leavesWithObservation(call,
                                    PresenceObservation::ABSENT), ',')
               << '\t' << join(leavesWithObservation(call,
                                    PresenceObservation::MISSING), ',')
               << '\t' << (call.presence.isAmbiguous() ? 1 : 0)
               << '\t' << (call.orientationConflict ? 1 : 0)
               << '\t' << (call.endpointConflict ? 1 : 0)
               << '\t' << (call.supportedWithoutLocalConflict ? 1 : 0)
               << '\t' << call.decision << '\t' << report.inferenceScope
               << "\tNOT_ATTEMPTED\t" << report.immutableCoreHash << '\n';
    }
}

void AncestralAdjacencyWriters::writeIssues(
    const AncestralAdjacencyReport &report, std::ostream &output) {
    output << "issue_id\tissue_kind\ttaxon\toccurrence_path\tblock_ids"
              "\tprojection_ids\tdetail\tinference_scope\tcore_hash\n";
    std::vector<BlockProjectionIssue> issues = report.issues;
    sortAndDeduplicateIssues(issues);
    for (const BlockProjectionIssue &issue : issues) {
        output << issue.issueId << '\t'
               << blockProjectionIssueKindName(issue.kind) << '\t'
               << (issue.taxon.empty() ? "." : issue.taxon) << '\t'
               << (issue.occurrencePath.empty() ? "." : issue.occurrencePath)
               << '\t' << (issue.blockIds.empty() ? "." : join(issue.blockIds, ','))
               << '\t'
               << (issue.projectionIds.empty() ? "." : join(issue.projectionIds, ','))
               << '\t' << cleanTsv(issue.detail) << '\t' << report.inferenceScope
               << '\t' << report.immutableCoreHash << '\n';
    }
}

}  // namespace trio
}  // namespace anchorwave
