#include "PairwiseManifestReader.h"

#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

const std::vector<std::string> &pairwiseHeader() {
    static const std::vector<std::string> header = {
        "taxon_a", "taxon_b", "maf", "anchor_map", "score_profile", "weight"};
    return header;
}

void requireNonMissing(const std::string &value,
                       const SourceLocation &location,
                       const std::string &fieldName) {
    if (value.empty() || value == ".") {
        throw ManifestError(location, fieldName + " must not be empty or '.'");
    }
}

void requirePathOrDot(const std::string &value,
                      const SourceLocation &location,
                      const std::string &fieldName) {
    if (value.empty()) {
        throw ManifestError(location, fieldName + " must be a path or '.'");
    }
}

}  // namespace

PairwiseManifestReader::PairwiseManifestReader(
    const PairwiseManifestReaderOptions &options)
    : options_(options) {}

PairwiseManifest PairwiseManifestReader::read(const std::string &path,
                                              const TaxonManifest &taxa) const {
    std::ifstream input(path.c_str());
    if (!input.good()) {
        throw ManifestError(SourceLocation(path, 0), "cannot open pairwise manifest");
    }

    PairwiseManifest manifest;
    std::map<CanonicalTaxonPair, SourceLocation> firstDefinition;
    bool headerSeen = false;
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;
        if (detail::isBlankOrComment(line)) {
            continue;
        }

        const SourceLocation location(path, lineNumber);
        const std::vector<std::string> fields = detail::splitTsv(line);
        if (!headerSeen) {
            detail::requireExactHeader(fields, pairwiseHeader(), location);
            headerSeen = true;
            continue;
        }
        if (fields.size() != pairwiseHeader().size()) {
            std::ostringstream message;
            message << "expected " << pairwiseHeader().size()
                    << " TSV columns, observed " << fields.size();
            throw ManifestError(location, message.str());
        }

        requireNonMissing(fields[0], location, "taxon_a");
        requireNonMissing(fields[1], location, "taxon_b");
        if (fields[0] == fields[1]) {
            throw ManifestError(location, "taxon_a and taxon_b must be distinct");
        }
        if (!taxa.contains(fields[0])) {
            throw ManifestError(location,
                                "taxon_a '" + fields[0] +
                                    "' is not declared in the taxon manifest");
        }
        if (!taxa.contains(fields[1])) {
            throw ManifestError(location,
                                "taxon_b '" + fields[1] +
                                    "' is not declared in the taxon manifest");
        }

        PairwiseManifestRecord record;
        record.sourceLocation = location;
        record.sourceTaxonA = fields[0];
        record.sourceTaxonB = fields[1];
        record.pair = CanonicalTaxonPair(fields[0], fields[1]);
        record.maf = fields[2];
        record.anchorMap = fields[3];
        record.scoreProfile = fields[4];
        requireNonMissing(record.maf, location, "maf");
        requirePathOrDot(record.anchorMap, location, "anchor_map");
        requireNonMissing(record.scoreProfile, location, "score_profile");
        record.maf = detail::resolveManifestEntryPath(path, record.maf, location);
        record.anchorMap =
            detail::resolveManifestEntryPath(path, record.anchorMap, location);
        record.weight = detail::parseFiniteDouble(fields[5], location, "weight");
        if (record.weight < 0.0) {
            throw ManifestError(location, "weight must be non-negative");
        }

        const std::map<CanonicalTaxonPair, SourceLocation>::const_iterator duplicate =
            firstDefinition.find(record.pair);
        if (duplicate != firstDefinition.end()) {
            std::ostringstream message;
            message << "duplicate unordered taxon pair " << record.pair.stableId()
                    << "; first defined at " << duplicate->second.path << ':'
                    << duplicate->second.line;
            throw ManifestError(location, message.str());
        }

        if (options_.validatePaths) {
            detail::validateExistingPath(record.maf, location, "maf");
            if (record.anchorMap != ".") {
                detail::validateExistingPath(record.anchorMap, location, "anchor_map");
            }
        }

        firstDefinition[record.pair] = location;
        manifest.records.push_back(record);
    }

    if (input.bad()) {
        throw ManifestError(SourceLocation(path, lineNumber),
                            "I/O failure while reading pairwise manifest");
    }
    if (!headerSeen) {
        throw ManifestError(SourceLocation(path, 0), "missing pairwise TSV header");
    }
    if (manifest.records.empty()) {
        throw ManifestError(SourceLocation(path, 0), "pairwise manifest has no data rows");
    }

    validateRequiredPairs(manifest, taxa, options_.scope, path);
    return manifest;
}

void PairwiseManifestReader::validateRequiredPairs(const PairwiseManifest &manifest,
                                                   const TaxonManifest &taxa,
                                                   const PairwiseScope scope,
                                                   const std::string &manifestPath) {
    const std::vector<TaxonId> ingroups = taxa.ingroupIds();
    const std::vector<TaxonId> outgroups = taxa.outgroupIds();
    if (ingroups.size() != 2 || outgroups.empty() ||
        taxa.primaryOutgroupId().empty()) {
        throw ManifestError(
            SourceLocation(manifestPath, 0),
            "taxon manifest must contain exactly two ingroups and one primary outgroup "
            "before pairwise edges can be validated");
    }

    std::set<CanonicalTaxonPair> required;
    required.insert(CanonicalTaxonPair(ingroups[0], ingroups[1]));
    if (scope == PairwiseScope::TRIANGLES) {
        for (const TaxonId &outgroup : outgroups) {
            required.insert(CanonicalTaxonPair(ingroups[0], outgroup));
            required.insert(CanonicalTaxonPair(ingroups[1], outgroup));
        }
    } else {
        const std::vector<TaxonId> allTaxa = taxa.allTaxonIds();
        for (std::size_t left = 0; left < allTaxa.size(); ++left) {
            for (std::size_t right = left + 1; right < allTaxa.size(); ++right) {
                required.insert(CanonicalTaxonPair(allTaxa[left], allTaxa[right]));
            }
        }
    }

    std::vector<std::string> missing;
    for (const CanonicalTaxonPair &pair : required) {
        if (!manifest.contains(pair)) {
            missing.push_back(pair.stableId());
        }
    }
    if (!missing.empty()) {
        std::ostringstream message;
        message << "missing required " << toString(scope) << " pairwise edge";
        if (missing.size() != 1) {
            message << 's';
        }
        message << ": ";
        for (std::size_t index = 0; index < missing.size(); ++index) {
            if (index != 0) {
                message << ", ";
            }
            message << missing[index];
        }
        throw ManifestError(SourceLocation(manifestPath, 0), message.str());
    }
}

}  // namespace trio
}  // namespace anchorwave
