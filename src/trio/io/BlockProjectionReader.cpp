#include "src/trio/io/BlockProjectionReader.h"

#include "src/trio/model/StableId.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

const char *const SCHEMA_VERSION = "#anchorwave-block-projections-v1";

const std::vector<std::string> &header() {
    static const std::vector<std::string> value = {
        "projection_id", "taxon_id", "sequence", "block_id",
        "start0", "end0", "orientation", "source"};
    return value;
}

void requireValue(const std::string &value, const SourceLocation &location,
                  const std::string &name, bool dotAllowed = false) {
    if (value.empty() || (!dotAllowed && value == ".")) {
        throw ManifestError(location, name + " must not be empty or '.'");
    }
}

BlockOrientation parseOrientation(const std::string &value,
                                  const SourceLocation &location) {
    if (value == "+") return BlockOrientation::FORWARD;
    if (value == "-") return BlockOrientation::REVERSE;
    if (value == "?") return BlockOrientation::UNKNOWN;
    throw ManifestError(location, "orientation must be '+', '-', or '?'");
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

}  // namespace

BlockProjectionSet BlockProjectionReader::read(
    const std::string &path, const TaxonManifest &taxa,
    const AlignmentSiteGraph &graph) {
    std::ifstream input(path.c_str());
    if (!input.good()) {
        throw ManifestError(SourceLocation(path, 0),
                            "cannot open block projection manifest");
    }
    bool versionSeen = false;
    bool headerSeen = false;
    std::size_t lineNumber = 0;
    std::string line;
    std::set<std::string> identifiers;
    BlockProjectionSet result;
    result.immutableCoreHash = graph.coreHash();
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string trimmed = detail::trimAscii(line);
        if (trimmed.empty()) continue;
        const SourceLocation location(path, lineNumber);
        if (!versionSeen) {
            if (trimmed != SCHEMA_VERSION) {
                throw ManifestError(
                    location, std::string("first non-empty line must be schema version '") +
                                  SCHEMA_VERSION + "'");
            }
            versionSeen = true;
            continue;
        }
        if (trimmed[0] == '#') continue;
        const std::vector<std::string> fields = detail::splitTsv(line);
        if (!headerSeen) {
            detail::requireExactHeader(fields, header(), location);
            headerSeen = true;
            continue;
        }
        if (fields.size() != header().size()) {
            std::ostringstream message;
            message << "expected " << header().size() << " TSV columns, observed "
                    << fields.size();
            throw ManifestError(location, message.str());
        }
        requireValue(fields[0], location, "projection_id", true);
        requireValue(fields[1], location, "taxon_id");
        requireValue(fields[2], location, "sequence");
        requireValue(fields[3], location, "block_id");
        requireValue(fields[7], location, "source");
        if (!taxa.contains(fields[1])) {
            throw ManifestError(location, "taxon_id '" + fields[1] +
                                             "' is not declared in the taxon manifest");
        }
        const int64_t start0 = detail::parseInt64(fields[4], location, "start0");
        const int64_t end0 = detail::parseInt64(fields[5], location, "end0");
        if (start0 < 0 || end0 <= start0) {
            throw ManifestError(location,
                                "block coordinates must be a non-empty, zero-based, "
                                "half-open interval");
        }
        ExtantBlockProjection projection;
        projection.taxon = fields[1];
        projection.occurrencePath = fields[1] + "|" + fields[2];
        projection.blockId = fields[3];
        projection.start0 = start0;
        projection.end0 = end0;
        projection.orientation = parseOrientation(fields[6], location);
        projection.projectionId = fields[0] == "."
                                      ? stableId("block_projection",
                                                 {projection.taxon,
                                                  projection.occurrencePath,
                                                  projection.blockId,
                                                  std::to_string(start0),
                                                  std::to_string(end0), fields[6]})
                                      : fields[0];
        if (!identifiers.insert(projection.projectionId).second) {
            throw ManifestError(location, "duplicate projection_id '" +
                                             projection.projectionId + "'");
        }
        result.projections.push_back(projection);
    }
    if (input.bad()) {
        throw ManifestError(SourceLocation(path, lineNumber),
                            "I/O failure while reading block projection manifest");
    }
    if (!versionSeen) {
        throw ManifestError(SourceLocation(path, 0),
                            std::string("missing schema version '") + SCHEMA_VERSION + "'");
    }
    if (!headerSeen) {
        throw ManifestError(SourceLocation(path, 0),
                            "missing block projection TSV header");
    }
    if (result.projections.empty()) {
        throw ManifestError(SourceLocation(path, 0),
                            "block projection manifest contains no data rows");
    }
    std::sort(result.projections.begin(), result.projections.end(), projectionLess);
    try {
        result.validateAgainst(graph);
    } catch (const AncestralAdjacencyError &error) {
        throw ManifestError(SourceLocation(path, 0), error.what());
    }
    return result;
}

}  // namespace trio
}  // namespace anchorwave
