#include "src/trio/evolution/AncestralInference.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <deque>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

unsigned int nucleotideMask(char value) {
    switch (value) {
        case 'A': case 'a': return 0x1u;
        case 'C': case 'c': return 0x2u;
        case 'G': case 'g': return 0x4u;
        case 'T': case 't': case 'U': case 'u': return 0x8u;
        case 'R': case 'r': return 0x5u;
        case 'Y': case 'y': return 0xAu;
        case 'S': case 's': return 0x6u;
        case 'W': case 'w': return 0x9u;
        case 'K': case 'k': return 0xCu;
        case 'M': case 'm': return 0x3u;
        case 'B': case 'b': return 0xEu;
        case 'D': case 'd': return 0xDu;
        case 'H': case 'h': return 0xBu;
        case 'V': case 'v': return 0x7u;
        default: return 0xFu;
    }
}

char maskToIupac(unsigned int mask) {
    switch (mask) {
        case 0x1u: return 'A'; case 0x2u: return 'C'; case 0x4u: return 'G';
        case 0x8u: return 'T'; case 0x5u: return 'R'; case 0xAu: return 'Y';
        case 0x6u: return 'S'; case 0x9u: return 'W'; case 0xCu: return 'K';
        case 0x3u: return 'M'; case 0xEu: return 'B'; case 0xDu: return 'D';
        case 0xBu: return 'H'; case 0x7u: return 'V'; default: return 'N';
    }
}

std::string taxonFromOccurrencePath(const std::string &path) {
    const std::size_t separator = path.find('|');
    return separator == std::string::npos ? path : path.substr(0, separator);
}

std::map<std::string, char> nucleotideObservations(
    const AlignmentSite &site, const std::set<std::string> &treeLeaves) {
    std::map<std::string, unsigned int> masks;
    for (const ResidueObservation &observation : site.observations) {
        if (!treeLeaves.count(observation.id.taxon)) {
            throw AncestralInferenceError("graph taxon '" + observation.id.taxon +
                                          "' is absent from the species tree");
        }
        masks[observation.id.taxon] |= nucleotideMask(observation.base);
    }
    std::map<std::string, char> result;
    for (const auto &entry : masks) result[entry.first] = maskToIupac(entry.second);
    return result;
}

std::map<std::string, PresenceObservation> presenceObservations(
    const AlignmentSite &site, const std::set<std::string> &treeLeaves) {
    std::map<std::string, PresenceObservation> result;
    for (const std::string &leaf : treeLeaves) result[leaf] = PresenceObservation::MISSING;
    for (const std::string &path : site.alignedAbsentPaths) {
        const std::string taxon = taxonFromOccurrencePath(path);
        if (!treeLeaves.count(taxon)) {
            throw AncestralInferenceError("gap evidence taxon '" + taxon +
                                          "' is absent from the species tree");
        }
        result[taxon] = PresenceObservation::ABSENT;
    }
    for (const ResidueObservation &observation : site.observations) {
        if (!treeLeaves.count(observation.id.taxon)) {
            throw AncestralInferenceError("graph taxon '" + observation.id.taxon +
                                          "' is absent from the species tree");
        }
        result[observation.id.taxon] = PresenceObservation::PRESENT;
    }
    return result;
}

struct Topology {
    std::map<std::string, std::set<std::string>> directed;
    std::map<std::string, std::set<std::string>> undirected;
};

Topology graphTopology(const AlignmentSiteGraph &graph) {
    Topology result;
    for (const auto &site : graph.sites()) {
        result.directed[site.first]; result.undirected[site.first];
    }
    for (const PathSegment &path : graph.pathSegments()) {
        for (std::size_t i = 1; i < path.siteIds.size(); ++i) {
            const std::string &a = path.siteIds[i - 1];
            const std::string &b = path.siteIds[i];
            if (a == b) continue;
            result.directed[a].insert(b);
            result.undirected[a].insert(b);
            result.undirected[b].insert(a);
        }
    }
    return result;
}

std::vector<std::vector<std::string>> connectedComponents(const Topology &topology) {
    std::set<std::string> unseen;
    for (const auto &entry : topology.undirected) unseen.insert(entry.first);
    std::vector<std::vector<std::string>> result;
    while (!unseen.empty()) {
        const std::string seed = *unseen.begin(); unseen.erase(unseen.begin());
        std::deque<std::string> queue(1, seed);
        std::vector<std::string> component;
        while (!queue.empty()) {
            const std::string current = queue.front(); queue.pop_front();
            component.push_back(current);
            for (const std::string &next : topology.undirected.at(current)) {
                auto found = unseen.find(next);
                if (found != unseen.end()) {
                    unseen.erase(found); queue.push_back(next);
                }
            }
        }
        result.push_back(component);
    }
    return result;
}

bool uniqueOrder(const std::vector<std::string> &component,
                 const Topology &topology, std::vector<std::string> &order) {
    const std::set<std::string> included(component.begin(), component.end());
    std::map<std::string, std::size_t> indegree;
    for (const std::string &site : component) indegree[site] = 0;
    for (const std::string &site : component) {
        for (const std::string &next : topology.directed.at(site)) {
            if (included.count(next)) ++indegree[next];
        }
    }
    std::set<std::string> ready;
    for (const auto &entry : indegree) if (entry.second == 0) ready.insert(entry.first);
    bool unique = true; order.clear();
    while (!ready.empty()) {
        if (ready.size() != 1) unique = false;
        const std::string current = *ready.begin(); ready.erase(ready.begin());
        order.push_back(current);
        for (const std::string &next : topology.directed.at(current)) {
            if (included.count(next) && --indegree[next] == 0) ready.insert(next);
        }
    }
    return unique && order.size() == component.size();
}

const ResidueObservation *observationForPath(const AlignmentSite &site,
                                             const std::string &path) {
    for (const ResidueObservation &observation : site.observations) {
        if (observation.id.occurrencePath == path) return &observation;
    }
    return nullptr;
}

std::set<std::string> knownIngroupPaths(const AlignmentSite &site,
                                        const AlignmentSiteGraph &graph) {
    std::set<std::string> result;
    for (const ResidueObservation &observation : site.observations) {
        if (observation.id.taxon == graph.ingroup1() ||
            observation.id.taxon == graph.ingroup2()) {
            result.insert(observation.id.occurrencePath);
        }
    }
    for (const std::string &path : site.alignedAbsentPaths) {
        const std::string taxon = taxonFromOccurrencePath(path);
        if (taxon == graph.ingroup1() || taxon == graph.ingroup2()) result.insert(path);
    }
    return result;
}

}  // namespace

AncestralOverlay AncestralInference::infer(
    const AlignmentSiteGraph &graph, const Phylogeny &tree,
    const EvolutionInferenceConfig &config) {
    if (config.targetNode.empty()) {
        throw AncestralInferenceError("target ancestor node must be named");
    }
    tree.internalNodeId(config.targetNode);
    const std::vector<std::string> leafVector = tree.leafNames();
    const std::set<std::string> leaves(leafVector.begin(), leafVector.end());
    AncestralOverlay overlay;
    overlay.immutableCoreHash = graph.coreHash();
    overlay.targetNode = config.targetNode;
    for (const auto &entry : graph.sites()) {
        const AlignmentSite &site = entry.second;
        AncestralSiteCall call;
        call.siteId = site.siteId;
        call.copyGroup = site.copyGroup;
        call.nucleotide = tree.inferNucleotide(
            config.targetNode, nucleotideObservations(site, leaves),
            config.nucleotideModel);
        call.presence = tree.inferPresence(
            config.targetNode, presenceObservations(site, leaves),
            config.presenceModel);
        if (config.requireCopyResolved && !site.copyResolved) {
            call.decision = "COPY_UNRESOLVED";
        } else if (call.presence.isAmbiguous() && !config.includeAmbiguousPresence) {
            call.decision = "PRESENCE_AMBIGUOUS";
        } else if (call.presence.isAmbiguous()) {
            // A binary tie contains both ABSENT and PRESENT. The explicit
            // opt-in means retain the present alternative, while preserving a
            // decision label that prevents it being mistaken for a resolved
            // presence call.
            call.emitted = true;
            call.emittedBase = call.nucleotide.ambiguityCode;
            call.decision = "PRESENCE_AMBIGUOUS_EMITTED";
        } else if (call.presence.selectedState == 0) {
            call.decision = "ANCESTRAL_ABSENT";
        } else {
            call.emitted = true;
            call.emittedBase = call.nucleotide.ambiguityCode;
            call.decision = call.nucleotide.isAmbiguous() ? "BASE_AMBIGUOUS" : "EMITTED";
        }
        overlay.calls[site.siteId] = call;
    }

    const Topology topology = graphTopology(graph);
    for (const std::vector<std::string> &component : connectedComponents(topology)) {
        std::vector<std::string> order;
        if (!uniqueOrder(component, topology, order)) {
            order = component;
            std::sort(order.begin(), order.end());
            // Non-unique site order is never converted into an unsupported
            // sequence order: each emitted site becomes its own local block.
            for (const std::string &siteId : order) {
                const AncestralSiteCall &call = overlay.calls.at(siteId);
                if (!call.emitted) continue;
                AncestralBlock block;
                block.copyGroup = call.copyGroup;
                block.siteIds.push_back(siteId);
                block.sequence.push_back(call.emittedBase);
                block.blockId = stableId("ancestor_block",
                                         {overlay.immutableCoreHash, config.targetNode,
                                          call.copyGroup, siteId});
                overlay.blocks.push_back(block);
            }
            continue;
        }
        std::size_t begin = 0;
        while (begin < order.size()) {
            while (begin < order.size() && !overlay.calls.at(order[begin]).emitted) ++begin;
            if (begin == order.size()) break;
            std::size_t end = begin;
            // Sequence continuity is determined by the immutable graph
            // topology, not by the copy label attached to each alignment
            // site. In unconstrained evidence each homologous column can have
            // a distinct inferred copy label; splitting on that label would
            // incorrectly emit one-base ancestral FASTA records. Structural
            // adjacency/karyotype inference remains a separate layer.
            while (end < order.size() && overlay.calls.at(order[end]).emitted) ++end;
            AncestralBlock block;
            block.copyGroup = overlay.calls.at(order[begin]).copyGroup;
            for (std::size_t i = begin; i < end; ++i) {
                if (overlay.calls.at(order[i]).copyGroup != block.copyGroup) {
                    // Empty means that this sequence block spans more than one
                    // resolved per-site copy label. The complete labels remain
                    // available in ancestor.calls.tsv.
                    block.copyGroup.clear();
                }
                block.siteIds.push_back(order[i]);
                block.sequence.push_back(overlay.calls.at(order[i]).emittedBase);
            }
            block.blockId = stableId("ancestor_block",
                                     {overlay.immutableCoreHash, config.targetNode,
                                      block.copyGroup, block.siteIds.front(),
                                      block.siteIds.back()});
            overlay.blocks.push_back(block);
            begin = end;
        }
    }
    std::sort(overlay.blocks.begin(), overlay.blocks.end(),
              [](const AncestralBlock &a, const AncestralBlock &b) {
                  return a.blockId < b.blockId;
              });
    overlay.validateAgainst(graph);
    return overlay;
}

void AncestralOverlay::validateAgainst(const AlignmentSiteGraph &graph) const {
    if (immutableCoreHash != graph.coreHash()) {
        throw AncestralInferenceError("ancestral overlay core hash does not match graph");
    }
    std::set<std::string> used;
    for (const AncestralBlock &block : blocks) {
        if (block.siteIds.size() != block.sequence.size() || block.sequence.empty()) {
            throw AncestralInferenceError("ancestral block spelling invariant failed");
        }
        for (std::size_t i = 0; i < block.siteIds.size(); ++i) {
            if (!graph.sites().count(block.siteIds[i]) ||
                !calls.at(block.siteIds[i]).emitted ||
                calls.at(block.siteIds[i]).emittedBase != block.sequence[i]) {
                throw AncestralInferenceError("ancestral block references an invalid call");
            }
            if (!used.insert(block.siteIds[i]).second) {
                throw AncestralInferenceError("ancestral site appears in multiple blocks");
            }
        }
    }
}

std::vector<AncestorChildAlignmentBlock>
AncestralInference::projectToIngroupChildren(
    const AlignmentSiteGraph &graph, const AncestralOverlay &overlay) {
    overlay.validateAgainst(graph);
    std::map<std::string, const ResidueObservation *> metadata;
    for (const auto &site : graph.sites()) {
        for (const ResidueObservation &observation : site.second.observations) {
            metadata.insert(std::make_pair(observation.id.occurrencePath, &observation));
        }
    }
    std::vector<AncestorChildAlignmentBlock> result;
    for (const AncestralBlock &block : overlay.blocks) {
        std::size_t begin = 0;
        while (begin < block.siteIds.size()) {
            const std::set<std::string> paths =
                knownIngroupPaths(graph.sites().at(block.siteIds[begin]), graph);
            std::size_t end = begin + 1;
            while (end < block.siteIds.size() &&
                   knownIngroupPaths(graph.sites().at(block.siteIds[end]), graph) == paths) {
                ++end;
            }
            if (paths.empty()) { begin = end; continue; }
            AncestorChildAlignmentBlock projection;
            projection.ancestorBlockId = block.blockId;
            projection.ancestorStart0 = static_cast<int64_t>(begin);
            projection.siteIds.assign(block.siteIds.begin() + static_cast<std::ptrdiff_t>(begin),
                                      block.siteIds.begin() + static_cast<std::ptrdiff_t>(end));
            projection.ancestorText.assign(block.sequence.begin() + static_cast<std::ptrdiff_t>(begin),
                                           block.sequence.begin() + static_cast<std::ptrdiff_t>(end));
            for (const std::string &path : paths) {
                const auto info = metadata.find(path);
                if (info == metadata.end()) continue;
                ChildProjectionRow row;
                row.occurrencePath = path;
                row.taxon = info->second->id.taxon;
                row.source = info->second->id.sequence;
                row.sourceSize = info->second->sourceSize;
                int64_t previous = -1;
                for (const std::string &siteId : projection.siteIds) {
                    const AlignmentSite &site = graph.sites().at(siteId);
                    const ResidueObservation *observation = observationForPath(site, path);
                    if (observation != nullptr) {
                        if (row.size == 0) row.start0 = observation->id.forwardPosition0;
                        if (previous >= 0 && observation->id.forwardPosition0 != previous + 1) {
                            row.text.clear(); break;
                        }
                        previous = observation->id.forwardPosition0;
                        ++row.size; row.text.push_back(observation->base);
                    } else if (site.alignedAbsentPaths.count(path)) {
                        row.text.push_back('-');
                    } else {
                        row.text.clear(); break;
                    }
                }
                if (row.text.size() != projection.siteIds.size()) continue;
                row.emptyComponent = row.size == 0;
                projection.children.push_back(row);
            }
            if (!projection.children.empty()) {
                std::sort(projection.children.begin(), projection.children.end(),
                          [](const ChildProjectionRow &a, const ChildProjectionRow &b) {
                              return a.occurrencePath < b.occurrencePath;
                          });
                std::vector<std::string> idFields = projection.siteIds;
                idFields.push_back(block.blockId);
                for (const ChildProjectionRow &row : projection.children)
                    idFields.push_back(row.occurrencePath);
                projection.projectionId = stableId("ancestor_child", idFields);
                result.push_back(projection);
            }
            begin = end;
        }
    }
    return result;
}

void AncestralWriters::writeCalls(const AncestralOverlay &overlay,
                                  std::ostream &output) {
    output << "site_id\tcopy_group\ttarget_node\temitted\temitted_base\tdecision"
              "\tnucleotide_states\tnucleotide_best_cost\tnucleotide_margin"
              "\tpresence_states\tpresence_best_cost\tpresence_margin\tcore_hash\n";
    for (const auto &entry : overlay.calls) {
        const AncestralSiteCall &call = entry.second;
        std::string bases(call.nucleotide.optimalStates.begin(),
                          call.nucleotide.optimalStates.end());
        std::ostringstream presence;
        for (int state : call.presence.optimalStates) presence << state;
        output << call.siteId << '\t' << (call.copyGroup.empty() ? "." : call.copyGroup)
               << '\t' << overlay.targetNode << '\t' << (call.emitted ? 1 : 0)
               << '\t' << (call.emitted ? std::string(1, call.emittedBase) : ".")
               << '\t' << call.decision << '\t' << bases << '\t'
               << call.nucleotide.bestCost << '\t' << call.nucleotide.confidenceMargin
               << '\t' << presence.str() << '\t' << call.presence.bestCost << '\t'
               << call.presence.confidenceMargin << '\t' << overlay.immutableCoreHash << '\n';
    }
}

void AncestralWriters::writeFasta(const AncestralOverlay &overlay,
                                  std::ostream &output) {
    for (const AncestralBlock &block : overlay.blocks) {
        output << '>' << block.blockId << " target=" << overlay.targetNode
               << " copy_group=" << (block.copyGroup.empty() ? "." : block.copyGroup)
               << " core_hash=" << overlay.immutableCoreHash << '\n';
        for (std::size_t i = 0; i < block.sequence.size(); i += 80) {
            output << block.sequence.substr(i, 80) << '\n';
        }
    }
}

void AncestralWriters::writeChildMaf(
    const AncestralOverlay &overlay,
    const std::vector<AncestorChildAlignmentBlock> &projections,
    std::ostream &output) {
    output << "##maf version=1 scoring=TrioAnchorGraph-ancestor-projection\n"
           << "# immutable_core_hash=" << overlay.immutableCoreHash << "\n\n";
    std::map<std::string, std::size_t> blockSizes;
    for (const AncestralBlock &block : overlay.blocks)
        blockSizes[block.blockId] = block.sequence.size();
    for (const AncestorChildAlignmentBlock &projection : projections) {
        output << "a projection_id=" << projection.projectionId << '\n'
               << "s\t" << overlay.targetNode << '.' << projection.ancestorBlockId
               << '\t' << projection.ancestorStart0 << '\t'
               << projection.ancestorText.size() << "\t+\t"
               << blockSizes.at(projection.ancestorBlockId) << '\t'
               << projection.ancestorText << '\n';
        for (const ChildProjectionRow &row : projection.children) {
            if (row.emptyComponent) {
                output << "e\t" << row.taxon << '.' << row.source << '\t'
                       << row.start0 << "\t0\t+\t" << row.sourceSize << "\tC\n";
            } else {
                output << "s\t" << row.taxon << '.' << row.source << '\t'
                       << row.start0 << '\t' << row.size << "\t+\t"
                       << row.sourceSize << '\t' << row.text << '\n';
            }
        }
        output << '\n';
    }
}

void AncestralWriters::writeChildMap(
    const std::vector<AncestorChildAlignmentBlock> &projections,
    std::ostream &output) {
    output << "projection_id\tancestor_block\tancestor_start0\tancestor_end0"
              "\tchild_taxon\tchild_path\tchild_source\tchild_start0\tchild_end0"
              "\tstate\tsite_ids\n";
    for (const AncestorChildAlignmentBlock &projection : projections) {
        std::ostringstream sites;
        for (std::size_t i = 0; i < projection.siteIds.size(); ++i) {
            if (i) sites << ','; sites << projection.siteIds[i];
        }
        for (const ChildProjectionRow &row : projection.children) {
            output << projection.projectionId << '\t' << projection.ancestorBlockId
                   << '\t' << projection.ancestorStart0 << '\t'
                   << projection.ancestorStart0 +
                          static_cast<int64_t>(projection.ancestorText.size())
                   << '\t' << row.taxon << '\t' << row.occurrencePath << '\t'
                   << row.source << '\t' << row.start0 << '\t' << row.start0 + row.size
                   << '\t' << (row.emptyComponent ? "ALIGNED_ABSENCE" : "ALIGNED")
                   << '\t' << sites.str() << '\n';
        }
    }
}

}  // namespace trio
}  // namespace anchorwave
