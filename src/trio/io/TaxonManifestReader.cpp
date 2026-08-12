#include "TaxonManifestReader.h"

#include <cctype>
#include <fstream>
#include <map>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

const std::vector<std::string> &taxonHeader() {
    static const std::vector<std::string> header = {
        "taxon_id",       "role",          "fasta",          "gff",
        "anchor_sam",     "anchor_fasta",  "callability_bed", "quality_weight"};
    return header;
}

void requireToken(const std::string &value,
                  const SourceLocation &location,
                  const std::string &fieldName,
                  const bool allowDot) {
    if (value.empty() || (!allowDot && value == ".")) {
        throw ManifestError(location,
                            fieldName + (allowDot ? " must be a path or '.'"
                                                  : " must not be empty or '.'"));
    }
}

void validateTaxonId(const std::string &value, const SourceLocation &location) {
    for (std::string::const_iterator character = value.begin();
         character != value.end(); ++character) {
        const unsigned char byte = static_cast<unsigned char>(*character);
        if (*character == '|' || std::isspace(byte) || std::iscntrl(byte)) {
            throw ManifestError(
                location,
                "taxon_id must not contain '|', whitespace, or control characters");
        }
    }
}

}  // namespace

TaxonManifestReader::TaxonManifestReader(const ManifestReaderOptions &options)
    : options_(options) {}

TaxonManifest TaxonManifestReader::read(const std::string &path) const {
    std::ifstream input(path.c_str());
    if (!input.good()) {
        throw ManifestError(SourceLocation(path, 0), "cannot open taxon manifest");
    }

    TaxonManifest manifest;
    std::map<TaxonId, SourceLocation> firstDefinition;
    bool headerSeen = false;
    std::string line;
    std::size_t lineNumber = 0;
    std::size_t ingroupCount = 0;
    std::size_t preferredIngroupCount = 0;
    std::size_t primaryOutgroupCount = 0;

    while (std::getline(input, line)) {
        ++lineNumber;
        if (detail::isBlankOrComment(line)) {
            continue;
        }

        const SourceLocation location(path, lineNumber);
        const std::vector<std::string> fields = detail::splitTsv(line);
        if (!headerSeen) {
            detail::requireExactHeader(fields, taxonHeader(), location);
            headerSeen = true;
            continue;
        }
        if (fields.size() != taxonHeader().size()) {
            std::ostringstream message;
            message << "expected " << taxonHeader().size() << " TSV columns, observed "
                    << fields.size();
            throw ManifestError(location, message.str());
        }

        TaxonRecord record;
        record.sourceLocation = location;
        record.taxonId = fields[0];
        requireToken(record.taxonId, location, "taxon_id", false);
        validateTaxonId(record.taxonId, location);
        if (firstDefinition.find(record.taxonId) != firstDefinition.end()) {
            const SourceLocation &original = firstDefinition[record.taxonId];
            std::ostringstream message;
            message << "duplicate taxon_id '" << record.taxonId << "'; first defined at "
                    << original.path << ':' << original.line;
            throw ManifestError(location, message.str());
        }

        record.role = parseTaxonRole(fields[1], location);
        record.fasta = fields[2];
        record.gff = fields[3];
        record.anchorSam = fields[4];
        record.anchorFasta = fields[5];
        record.callabilityBed = fields[6];
        requireToken(record.fasta, location, "fasta", false);
        requireToken(record.gff, location, "gff", true);
        requireToken(record.anchorSam, location, "anchor_sam", true);
        requireToken(record.anchorFasta, location, "anchor_fasta", true);
        requireToken(record.callabilityBed, location, "callability_bed", true);

        record.fasta = detail::resolveManifestEntryPath(path, record.fasta, location);
        record.gff = detail::resolveManifestEntryPath(path, record.gff, location);
        record.anchorSam =
            detail::resolveManifestEntryPath(path, record.anchorSam, location);
        record.anchorFasta =
            detail::resolveManifestEntryPath(path, record.anchorFasta, location);
        record.callabilityBed =
            detail::resolveManifestEntryPath(path, record.callabilityBed, location);

        record.qualityWeight =
            detail::parseFiniteDouble(fields[7], location, "quality_weight");
        if (record.qualityWeight < 0.0) {
            throw ManifestError(location, "quality_weight must be non-negative");
        }

        if (options_.validatePaths) {
            detail::validateExistingPath(record.fasta, location, "fasta");
            if (record.gff != ".") {
                detail::validateExistingPath(record.gff, location, "gff");
            }
            if (record.anchorSam != ".") {
                detail::validateExistingPath(record.anchorSam, location, "anchor_sam");
            }
            if (record.anchorFasta != ".") {
                detail::validateExistingPath(record.anchorFasta, location, "anchor_fasta");
            }
            if (record.callabilityBed != ".") {
                detail::validateExistingPath(
                    record.callabilityBed, location, "callability_bed");
            }
        }

        if (record.role == TaxonRole::INGROUP_REFERENCE ||
            record.role == TaxonRole::INGROUP) {
            ++ingroupCount;
        }
        if (record.role == TaxonRole::INGROUP_REFERENCE) {
            ++preferredIngroupCount;
        }
        if (record.role == TaxonRole::PRIMARY_OUTGROUP) {
            ++primaryOutgroupCount;
        }

        firstDefinition[record.taxonId] = location;
        manifest.records.push_back(record);
    }

    if (input.bad()) {
        throw ManifestError(SourceLocation(path, lineNumber),
                            "I/O failure while reading taxon manifest");
    }
    if (!headerSeen) {
        throw ManifestError(SourceLocation(path, 0), "missing taxon TSV header");
    }
    if (manifest.records.empty()) {
        throw ManifestError(SourceLocation(path, 0), "taxon manifest has no data rows");
    }
    if (ingroupCount != 2) {
        std::ostringstream message;
        message << "expected exactly two ingroup taxa, observed " << ingroupCount;
        throw ManifestError(SourceLocation(path, 0), message.str());
    }
    if (preferredIngroupCount > 1) {
        throw ManifestError(SourceLocation(path, 0),
                            "at most one ingroup may have role ingroup_reference");
    }
    if (primaryOutgroupCount != 1) {
        std::ostringstream message;
        message << "expected exactly one primary_outgroup, observed "
                << primaryOutgroupCount;
        throw ManifestError(SourceLocation(path, 0), message.str());
    }

    return manifest;
}

}  // namespace trio
}  // namespace anchorwave
