#include "src/trio/io/GraphWriters.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <deque>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

std::string alleleNodeId(const std::string &siteId, char base) {
    return stableId("node", {siteId, std::string(1, base)});
}

std::string sanitizeTagValue(std::string value) {
    for (char &character : value) {
        if (character == '\t' || character == '\n' || character == '\r') {
            character = '_';
        }
    }
    return value;
}

const ResidueObservation *observationForPath(const AlignmentSite &site,
                                             const std::string &occurrencePath) {
    for (const ResidueObservation &observation : site.observations) {
        if (observation.id.occurrencePath == occurrencePath) {
            return &observation;
        }
    }
    return nullptr;
}

std::map<std::string, const ResidueObservation *> pathMetadata(
    const AlignmentSiteGraph &graph) {
    std::map<std::string, const ResidueObservation *> result;
    for (const auto &site : graph.sites()) {
        for (const ResidueObservation &observation : site.second.observations) {
            result.insert(std::make_pair(observation.id.occurrencePath, &observation));
        }
    }
    return result;
}

struct SiteTopology {
    std::map<std::string, std::set<std::string>> directed;
    std::map<std::string, std::set<std::string>> undirected;
};

SiteTopology topology(const AlignmentSiteGraph &graph) {
    SiteTopology result;
    for (const auto &site : graph.sites()) {
        result.directed[site.first];
        result.undirected[site.first];
    }
    for (const PathSegment &path : graph.pathSegments()) {
        for (std::size_t i = 1; i < path.siteIds.size(); ++i) {
            const std::string &left = path.siteIds[i - 1];
            const std::string &right = path.siteIds[i];
            if (left == right) {
                continue;
            }
            result.directed[left].insert(right);
            result.undirected[left].insert(right);
            result.undirected[right].insert(left);
        }
    }
    return result;
}

std::vector<std::vector<std::string>> components(const SiteTopology &topology) {
    std::set<std::string> unseen;
    for (const auto &entry : topology.undirected) {
        unseen.insert(entry.first);
    }
    std::vector<std::vector<std::string>> result;
    while (!unseen.empty()) {
        const std::string seed = *unseen.begin();
        unseen.erase(unseen.begin());
        std::deque<std::string> queue(1, seed);
        std::vector<std::string> component;
        while (!queue.empty()) {
            const std::string current = queue.front();
            queue.pop_front();
            component.push_back(current);
            const auto found = topology.undirected.find(current);
            if (found == topology.undirected.end()) {
                continue;
            }
            for (const std::string &next : found->second) {
                const auto present = unseen.find(next);
                if (present != unseen.end()) {
                    unseen.erase(present);
                    queue.push_back(next);
                }
            }
        }
        std::sort(component.begin(), component.end());
        result.push_back(component);
    }
    return result;
}

bool uniqueTopologicalOrder(const std::vector<std::string> &component,
                            const SiteTopology &topology,
                            std::vector<std::string> &order) {
    const std::set<std::string> included(component.begin(), component.end());
    std::map<std::string, std::size_t> indegree;
    for (const std::string &site : component) {
        indegree[site] = 0;
    }
    for (const std::string &site : component) {
        const auto found = topology.directed.find(site);
        if (found == topology.directed.end()) continue;
        for (const std::string &next : found->second) {
            if (included.count(next)) ++indegree[next];
        }
    }
    std::set<std::string> ready;
    for (const auto &entry : indegree) {
        if (entry.second == 0) ready.insert(entry.first);
    }
    order.clear();
    bool unique = true;
    while (!ready.empty()) {
        if (ready.size() != 1) unique = false;
        const std::string current = *ready.begin();
        ready.erase(ready.begin());
        order.push_back(current);
        const auto found = topology.directed.find(current);
        if (found == topology.directed.end()) continue;
        for (const std::string &next : found->second) {
            if (!included.count(next)) continue;
            if (--indegree[next] == 0) ready.insert(next);
        }
    }
    return unique && order.size() == component.size();
}

std::set<std::string> knownPaths(const AlignmentSite &site) {
    std::set<std::string> result = site.alignedAbsentPaths;
    for (const ResidueObservation &observation : site.observations) {
        result.insert(observation.id.occurrencePath);
    }
    return result;
}

struct MafRowProjection {
    std::string source;
    int64_t start0 = -1;
    int64_t size = 0;
    int64_t sourceSize = 0;
    std::string text;
};

bool projectRow(const std::vector<std::string> &sites,
                const std::string &path,
                const AlignmentSiteGraph &graph,
                const ResidueObservation &metadata,
                MafRowProjection &row,
                std::string &reason) {
    row.source = metadata.id.taxon + "." + metadata.id.sequence;
    row.sourceSize = metadata.sourceSize;
    int64_t previous = -1;
    for (const std::string &siteId : sites) {
        const AlignmentSite &site = graph.sites().at(siteId);
        const ResidueObservation *observation = observationForPath(site, path);
        if (observation != nullptr) {
            if (row.start0 < 0) row.start0 = observation->id.forwardPosition0;
            if (previous >= 0 && observation->id.forwardPosition0 != previous + 1) {
                reason = "noncontiguous_forward_coordinates";
                return false;
            }
            previous = observation->id.forwardPosition0;
            ++row.size;
            row.text.push_back(observation->base);
        } else if (site.alignedAbsentPaths.count(path)) {
            row.text.push_back('-');
        } else {
            reason = "missing_coverage_state";
            return false;
        }
    }
    if (row.size == 0) {
        reason = "all_gap_row_requires_MAF_empty_component";
        return false;
    }
    return true;
}

}  // namespace

void GfaWriter::write(const AlignmentSiteGraph &graph, std::ostream &output) {
    output << "H\tVN:Z:1.1\tAW:Z:TrioAnchorGraph-v1\tCH:Z:"
           << graph.coreHash() << '\n';
    for (const auto &entry : graph.sites()) {
        std::set<char> alleles;
        for (const ResidueObservation &observation : entry.second.observations) {
            alleles.insert(observation.base);
        }
        for (char allele : alleles) {
            output << "S\t" << alleleNodeId(entry.first, allele) << '\t' << allele
                   << "\tSI:Z:" << entry.first
                   << "\tCL:Z:" << consistencyClassName(entry.second.consistency);
            if (entry.second.copyResolved) {
                output << "\tCG:Z:" << sanitizeTagValue(entry.second.copyGroup);
            }
            output << '\n';
        }
    }
    std::set<std::pair<std::string, std::string>> links;
    for (const PathSegment &path : graph.pathSegments()) {
        std::vector<std::string> nodes;
        for (const std::string &siteId : path.siteIds) {
            const ResidueObservation *observation =
                observationForPath(graph.sites().at(siteId), path.occurrencePath);
            if (observation == nullptr) {
                throw AlignmentGraphError("GFA path cannot find its site observation");
            }
            nodes.push_back(alleleNodeId(siteId, observation->base));
        }
        for (std::size_t i = 1; i < nodes.size(); ++i) {
            links.insert(std::make_pair(nodes[i - 1], nodes[i]));
        }
        output << "P\t" << path.pathId << '\t';
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (i) output << ',';
            output << nodes[i] << '+';
        }
        output << "\t*\tTX:Z:" << sanitizeTagValue(path.taxon)
               << "\tSN:Z:" << sanitizeTagValue(path.sequence)
               << "\tSO:i:" << path.start0 << "\tEO:i:" << path.end0 << '\n';
    }
    // Emit links after paths in a deterministic section. GFA record order does
    // not define truth; all consumers must use the named handles.
    for (const auto &link : links) {
        output << "L\t" << link.first << "\t+\t" << link.second << "\t+\t0M\n";
    }
    if (!output) throw AlignmentGraphError("failed while writing GFA output");
}

void GfaWriter::writeMetadata(const AlignmentSiteGraph &graph, std::ostream &output) {
    output << "site_id\tallele_node_id\ttaxon\toccurrence_path\tsequence"
              "\tposition0\tsource_size\tstate\tbase\tcopy_group\tconsistency"
              "\tevidence_ids\n";
    for (const auto &entry : graph.sites()) {
        const AlignmentSite &site = entry.second;
        std::ostringstream evidence;
        for (std::size_t i = 0; i < site.evidenceIds.size(); ++i) {
            if (i) evidence << ',';
            evidence << site.evidenceIds[i];
        }
        for (const ResidueObservation &observation : site.observations) {
            output << site.siteId << '\t' << alleleNodeId(site.siteId, observation.base)
                   << '\t' << observation.id.taxon << '\t'
                   << observation.id.occurrencePath << '\t' << observation.id.sequence
                   << '\t' << observation.id.forwardPosition0 << '\t'
                   << observation.sourceSize << "\tOBSERVED_BASE\t" << observation.base
                   << '\t' << (site.copyResolved ? site.copyGroup : ".")
                   << '\t' << consistencyClassName(site.consistency) << '\t'
                   << evidence.str() << '\n';
        }
        for (const std::string &path : site.alignedAbsentPaths) {
            output << site.siteId << "\t.\t.\t" << path
                   << "\t.\t.\t.\tALIGNED_ABSENCE\t-\t"
                   << (site.copyResolved ? site.copyGroup : ".") << '\t'
                   << consistencyClassName(site.consistency) << '\t';
            for (std::size_t i = 0; i < site.absenceEvidenceIds.size(); ++i) {
                if (i) output << ',';
                output << site.absenceEvidenceIds[i];
            }
            output << '\n';
        }
    }
    if (!output) throw AlignmentGraphError("failed while writing graph metadata");
}

MafExportResult MafGraphExporter::write(const AlignmentSiteGraph &graph,
                                        std::ostream &output) {
    MafExportResult result;
    output << "##maf version=1 scoring=TrioAnchorGraph-projection\n"
           << "# authoritative_core_hash=" << graph.coreHash() << '\n'
           << "# missing coverage is split/omitted, never encoded as a gap\n\n";
    const SiteTopology graphTopology = topology(graph);
    const std::vector<std::vector<std::string>> graphComponents = components(graphTopology);
    const std::map<std::string, const ResidueObservation *> metadata = pathMetadata(graph);

    for (const std::vector<std::string> &component : graphComponents) {
        const std::string componentId = stableId("component", component);
        std::vector<std::string> order;
        if (!uniqueTopologicalOrder(component, graphTopology, order)) {
            result.omissions.push_back({componentId, "non_unique_or_cyclic_site_order", component});
            continue;
        }

        std::size_t begin = 0;
        while (begin < order.size()) {
            const std::set<std::string> paths = knownPaths(graph.sites().at(order[begin]));
            std::size_t end = begin + 1;
            while (end < order.size() && knownPaths(graph.sites().at(order[end])) == paths) {
                ++end;
            }
            const std::vector<std::string> run(order.begin() + static_cast<std::ptrdiff_t>(begin),
                                               order.begin() + static_cast<std::ptrdiff_t>(end));
            std::set<std::string> taxa;
            for (const std::string &path : paths) {
                const auto found = metadata.find(path);
                if (found != metadata.end()) taxa.insert(found->second->id.taxon);
            }
            if (paths.size() < 2 || taxa.size() < 2) {
                result.omissions.push_back({componentId, "fewer_than_two_explicit_paths", run});
                begin = end;
                continue;
            }

            std::vector<MafRowProjection> rows;
            std::string failure;
            for (const std::string &path : paths) {
                const auto found = metadata.find(path);
                if (found == metadata.end()) {
                    failure = "gap_only_path_has_no_source_metadata";
                    break;
                }
                MafRowProjection row;
                if (!projectRow(run, path, graph, *found->second, row, failure)) break;
                rows.push_back(row);
            }
            if (!failure.empty()) {
                result.omissions.push_back({componentId, failure, run});
                begin = end;
                continue;
            }
            output << "a graph_component=" << componentId
                   << " core_hash=" << graph.coreHash() << '\n';
            for (const MafRowProjection &row : rows) {
                output << "s\t" << row.source << '\t' << row.start0 << '\t'
                       << row.size << "\t+\t" << row.sourceSize << '\t'
                       << row.text << '\n';
            }
            output << '\n';
            ++result.blocksWritten;
            begin = end;
        }
    }
    if (!output) throw AlignmentGraphError("failed while writing MAF projection");
    return result;
}

void MafGraphExporter::writeOmissions(const MafExportResult &result,
                                      const std::string &coreHash,
                                      std::ostream &output) {
    output << "component_id\treason\tsite_count\tsite_ids\tcore_hash\n";
    for (const MafExportOmission &omission : result.omissions) {
        output << omission.componentId << '\t' << omission.reason << '\t'
               << omission.siteIds.size() << '\t';
        for (std::size_t index = 0; index < omission.siteIds.size(); ++index) {
            if (index) output << ',';
            output << omission.siteIds[index];
        }
        output << '\t' << coreHash << '\n';
    }
    if (!output) throw AlignmentGraphError("failed while writing MAF omission audit");
}

}  // namespace trio
}  // namespace anchorwave
