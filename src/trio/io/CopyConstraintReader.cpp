#include "CopyConstraintReader.h"

#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

const char *const COPY_SCHEMA_VERSION = "#anchorwave-copy-relations-v1";

const std::vector<std::string> &copyHeader() {
    static const std::vector<std::string> header = {
        "record_type",      "family_id",       "ancestral_copy_id",
        "parent_copy_id",   "born_on_edge",    "taxon_id",
        "member_type",      "member_id",       "seq",
        "start0",           "end0",            "strand",
        "expected_count",   "relation",        "confidence",
        "constraint",       "source"};
    return header;
}

std::string optionalValue(const std::string &text) {
    return detail::isMissingValue(text) ? std::string() : text;
}

void requireValue(const std::string &value,
                  const SourceLocation &location,
                  const std::string &fieldName) {
    if (detail::isMissingValue(value)) {
        throw ManifestError(location, fieldName + " must not be empty or '.'");
    }
}

void requireDot(const std::string &value,
                const SourceLocation &location,
                const std::string &fieldName,
                const CopyRecordType type) {
    if (value != ".") {
        throw ManifestError(location,
                            fieldName + " must be '.' for " + toString(type) +
                                " records");
    }
}

CopyRecordType parseRecordType(const std::string &value,
                               const SourceLocation &location) {
    if (value == "GROUP") {
        return CopyRecordType::GROUP;
    }
    if (value == "MEMBER") {
        return CopyRecordType::MEMBER;
    }
    if (value == "COUNT") {
        return CopyRecordType::COUNT;
    }
    throw ManifestError(location,
                        "record_type must be GROUP, MEMBER, or COUNT; observed '" +
                            value + "'");
}

CopyMemberType parseMemberType(const std::string &value,
                               const SourceLocation &location) {
    if (value == "anchor") {
        return CopyMemberType::ANCHOR;
    }
    if (value == "interval") {
        return CopyMemberType::INTERVAL;
    }
    throw ManifestError(location,
                        "member_type must be anchor or interval for MEMBER records; "
                        "observed '" +
                            value + "'");
}

CopyRelation parseRelation(const std::string &value,
                           const SourceLocation &location) {
    if (value == "ortholog") {
        return CopyRelation::ORTHOLOG;
    }
    if (value == "coortholog") {
        return CopyRelation::COORTHOLOG;
    }
    if (value == "excluded") {
        return CopyRelation::EXCLUDED;
    }
    throw ManifestError(location,
                        "relation must be ortholog, coortholog, or excluded for MEMBER "
                        "records; observed '" +
                            value + "'");
}

ConstraintStrength parseStrength(const std::string &value,
                                 const SourceLocation &location) {
    if (value == "hard") {
        return ConstraintStrength::HARD;
    }
    if (value == "soft") {
        return ConstraintStrength::SOFT;
    }
    throw ManifestError(location,
                        "constraint must be explicitly 'hard' or 'soft'; observed '" +
                            value + "'");
}

std::string memberOccurrenceKey(const CopyConstraint &constraint) {
    std::ostringstream key;
    key << constraint.familyId << '\t' << constraint.taxonId << '\t'
        << toString(constraint.memberType) << '\t';
    if (constraint.memberType == CopyMemberType::ANCHOR) {
        key << constraint.memberId;
    } else {
        key << constraint.interval.sequence << ':' << constraint.interval.start0 << '-'
            << constraint.interval.end0 << ':' << constraint.interval.strand;
    }
    return key.str();
}

void validateTaxon(const TaxonManifest &taxa,
                   const TaxonId &taxonId,
                   const SourceLocation &location,
                   const std::string &fieldName) {
    if (!taxa.contains(taxonId)) {
        throw ManifestError(location,
                            fieldName + " '" + taxonId +
                                "' is not declared in the taxon manifest");
    }
}

}  // namespace

CopyConstraintReader::CopyConstraintReader(const ManifestReaderOptions &options)
    : options_(options) {}

CopyConstraintSet CopyConstraintReader::read(
    const std::string &path,
    const TaxonManifest &taxa,
    const std::vector<std::string> &defaultCapacitySpecs) const {
    CopyConstraintSet result;

    for (std::size_t index = 0; index < defaultCapacitySpecs.size(); ++index) {
        const SourceLocation location("<copy-number>", index + 1);
        CopyCapacity capacity =
            parseDefaultCapacitySpec(defaultCapacitySpecs[index], location);
        validateTaxon(taxa, capacity.key.taxonId, location, "copy-number taxon");
        const CopyCapacityMap::const_iterator duplicate =
            result.defaultCapacities.find(capacity.key);
        if (duplicate != result.defaultCapacities.end()) {
            std::ostringstream message;
            message << "duplicate default copy capacity for ";
            if (!capacity.key.familyId.empty()) {
                message << capacity.key.familyId << ':';
            }
            message << capacity.key.taxonId << "; first defined at "
                    << duplicate->second.sourceLocation.path << ':'
                    << duplicate->second.sourceLocation.line;
            throw ManifestError(location, message.str());
        }
        result.defaultCapacities[capacity.key] = capacity;
    }

    if (path.empty() || path == ".") {
        return result;
    }
    if (options_.validatePaths) {
        detail::validateExistingPath(
            path, SourceLocation(path, 0), "copy relations manifest");
    }

    std::ifstream input(path.c_str());
    if (!input.good()) {
        throw ManifestError(SourceLocation(path, 0),
                            "cannot open copy relations manifest");
    }

    bool versionSeen = false;
    bool headerSeen = false;
    std::string line;
    std::size_t lineNumber = 0;
    std::map<std::string, SourceLocation> groups;
    std::map<std::string, SourceLocation> counts;
    std::map<std::string, SourceLocation> exactMembers;
    std::map<std::string, std::pair<AncestralCopyId, SourceLocation>> hardAssignments;
    std::map<std::string, std::pair<CopyRelation, SourceLocation>> hardGroupRelations;

    while (std::getline(input, line)) {
        ++lineNumber;
        const std::string trimmed = detail::trimAscii(line);
        if (trimmed.empty()) {
            continue;
        }

        const SourceLocation location(path, lineNumber);
        if (!versionSeen) {
            if (trimmed != COPY_SCHEMA_VERSION) {
                throw ManifestError(
                    location,
                    std::string("first non-empty line must be schema version '") +
                        COPY_SCHEMA_VERSION + "'");
            }
            versionSeen = true;
            continue;
        }
        if (trimmed[0] == '#') {
            continue;
        }

        const std::vector<std::string> fields = detail::splitTsv(line);
        if (!headerSeen) {
            detail::requireExactHeader(fields, copyHeader(), location);
            headerSeen = true;
            continue;
        }
        if (fields.size() != copyHeader().size()) {
            std::ostringstream message;
            message << "expected " << copyHeader().size() << " TSV columns, observed "
                    << fields.size();
            throw ManifestError(location, message.str());
        }

        CopyConstraint constraint;
        constraint.sourceLocation = location;
        constraint.recordType = parseRecordType(fields[0], location);
        constraint.familyId = optionalValue(fields[1]);
        constraint.ancestralCopyId = optionalValue(fields[2]);
        constraint.parentCopyId = optionalValue(fields[3]);
        constraint.bornOnEdge = optionalValue(fields[4]);
        constraint.taxonId = optionalValue(fields[5]);
        constraint.memberId = optionalValue(fields[7]);
        constraint.confidence =
            detail::parseFiniteDouble(fields[14], location, "confidence");
        if (constraint.confidence < 0.0 || constraint.confidence > 1.0) {
            throw ManifestError(location, "confidence must be in the closed interval [0,1]");
        }
        constraint.strength = parseStrength(fields[15], location);
        constraint.source = fields[16];
        requireValue(constraint.source, location, "source");

        requireValue(fields[1], location, "family_id");

        if (constraint.recordType == CopyRecordType::GROUP) {
            requireValue(fields[2], location, "ancestral_copy_id");
            requireDot(fields[5], location, "taxon_id", constraint.recordType);
            requireDot(fields[6], location, "member_type", constraint.recordType);
            requireDot(fields[7], location, "member_id", constraint.recordType);
            requireDot(fields[8], location, "seq", constraint.recordType);
            requireDot(fields[9], location, "start0", constraint.recordType);
            requireDot(fields[10], location, "end0", constraint.recordType);
            requireDot(fields[11], location, "strand", constraint.recordType);
            requireDot(fields[12], location, "expected_count", constraint.recordType);
            requireDot(fields[13], location, "relation", constraint.recordType);

            const std::string key =
                constraint.familyId + '\t' + constraint.ancestralCopyId;
            if (groups.find(key) != groups.end()) {
                const SourceLocation &original = groups[key];
                std::ostringstream message;
                message << "duplicate GROUP for family_id '" << constraint.familyId
                        << "' and ancestral_copy_id '" << constraint.ancestralCopyId
                        << "'; first defined at " << original.path << ':' << original.line;
                throw ManifestError(location, message.str());
            }
            groups[key] = location;
        } else if (constraint.recordType == CopyRecordType::COUNT) {
            requireValue(fields[5], location, "taxon_id");
            validateTaxon(taxa, constraint.taxonId, location, "taxon_id");
            requireDot(fields[2], location, "ancestral_copy_id", constraint.recordType);
            requireDot(fields[3], location, "parent_copy_id", constraint.recordType);
            requireDot(fields[4], location, "born_on_edge", constraint.recordType);
            requireDot(fields[6], location, "member_type", constraint.recordType);
            requireDot(fields[7], location, "member_id", constraint.recordType);
            requireDot(fields[8], location, "seq", constraint.recordType);
            requireDot(fields[9], location, "start0", constraint.recordType);
            requireDot(fields[10], location, "end0", constraint.recordType);
            requireDot(fields[11], location, "strand", constraint.recordType);
            requireValue(fields[12], location, "expected_count");
            requireDot(fields[13], location, "relation", constraint.recordType);
            constraint.hasExpectedCount = true;
            constraint.expectedCount =
                detail::parseInt64(fields[12], location, "expected_count");
            if (constraint.expectedCount < 0) {
                throw ManifestError(location, "expected_count must be non-negative");
            }

            const std::string key = constraint.familyId + '\t' + constraint.taxonId;
            if (counts.find(key) != counts.end()) {
                const SourceLocation &original = counts[key];
                std::ostringstream message;
                message << "duplicate COUNT for family_id '" << constraint.familyId
                        << "' and taxon_id '" << constraint.taxonId
                        << "'; first defined at " << original.path << ':' << original.line;
                throw ManifestError(location, message.str());
            }
            counts[key] = location;
        } else {
            requireValue(fields[2], location, "ancestral_copy_id");
            requireDot(fields[3], location, "parent_copy_id", constraint.recordType);
            requireDot(fields[4], location, "born_on_edge", constraint.recordType);
            requireValue(fields[5], location, "taxon_id");
            validateTaxon(taxa, constraint.taxonId, location, "taxon_id");
            constraint.memberType = parseMemberType(fields[6], location);
            constraint.relation = parseRelation(fields[13], location);
            requireDot(fields[12], location, "expected_count", constraint.recordType);

            if (constraint.memberType == CopyMemberType::ANCHOR) {
                requireValue(fields[7], location, "member_id");
                requireDot(fields[8], location, "seq", constraint.recordType);
                requireDot(fields[9], location, "start0", constraint.recordType);
                requireDot(fields[10], location, "end0", constraint.recordType);
                requireDot(fields[11], location, "strand", constraint.recordType);
            } else {
                requireDot(fields[7], location, "member_id", constraint.recordType);
                requireValue(fields[8], location, "seq");
                requireValue(fields[9], location, "start0");
                requireValue(fields[10], location, "end0");
                requireValue(fields[11], location, "strand");
                if (fields[11].size() != 1) {
                    throw ManifestError(location, "strand must be '+' or '-'");
                }
                constraint.hasInterval = true;
                constraint.interval = GenomicInterval(
                    fields[8], detail::parseInt64(fields[9], location, "start0"),
                    detail::parseInt64(fields[10], location, "end0"), fields[11][0]);
                validateHalfOpenInterval(constraint.interval, location, "member interval");
            }

            const std::string occurrence = memberOccurrenceKey(constraint);
            const std::string exactKey = occurrence + '\t' +
                                         constraint.ancestralCopyId + '\t' +
                                         toString(constraint.relation) + '\t' +
                                         toString(constraint.strength);
            if (exactMembers.find(exactKey) != exactMembers.end()) {
                const SourceLocation &original = exactMembers[exactKey];
                std::ostringstream message;
                message << "duplicate MEMBER constraint; first defined at "
                        << original.path << ':' << original.line;
                throw ManifestError(location, message.str());
            }
            exactMembers[exactKey] = location;

            if (constraint.strength == ConstraintStrength::HARD) {
                const std::string groupKey =
                    occurrence + '\t' + constraint.ancestralCopyId;
                const std::map<std::string,
                               std::pair<CopyRelation, SourceLocation>>::const_iterator
                    sameGroup = hardGroupRelations.find(groupKey);
                if (sameGroup != hardGroupRelations.end() &&
                    sameGroup->second.first != constraint.relation) {
                    std::ostringstream message;
                    message << "contradictory hard MEMBER relations for the same occurrence "
                               "and ancestral copy; first defined at "
                            << sameGroup->second.second.path << ':'
                            << sameGroup->second.second.line;
                    throw ManifestError(location, message.str());
                }
                hardGroupRelations[groupKey] =
                    std::make_pair(constraint.relation, location);

                if (constraint.relation != CopyRelation::EXCLUDED) {
                    const std::map<
                        std::string,
                        std::pair<AncestralCopyId, SourceLocation>>::const_iterator prior =
                        hardAssignments.find(occurrence);
                    if (prior != hardAssignments.end() &&
                        prior->second.first != constraint.ancestralCopyId) {
                        std::ostringstream message;
                        message << "one occurrence is assigned to incompatible hard ancestral "
                                   "copies '"
                                << prior->second.first << "' and '"
                                << constraint.ancestralCopyId << "'; first assignment at "
                                << prior->second.second.path << ':'
                                << prior->second.second.line;
                        throw ManifestError(location, message.str());
                    }
                    hardAssignments[occurrence] =
                        std::make_pair(constraint.ancestralCopyId, location);
                }
            }
        }

        result.records.push_back(constraint);
    }

    if (input.bad()) {
        throw ManifestError(SourceLocation(path, lineNumber),
                            "I/O failure while reading copy relations manifest");
    }
    if (!versionSeen) {
        throw ManifestError(SourceLocation(path, 0),
                            std::string("missing schema version '") +
                                COPY_SCHEMA_VERSION + "'");
    }
    if (!headerSeen) {
        throw ManifestError(SourceLocation(path, 0),
                            "missing copy relations TSV header");
    }

    return result;
}

CopyCapacity CopyConstraintReader::parseDefaultCapacitySpec(
    const std::string &specification,
    const SourceLocation &location) {
    const std::string value = detail::trimAscii(specification);
    const std::size_t equal = value.find('=');
    if (equal == std::string::npos || equal == 0 || equal + 1 >= value.size() ||
        value.find('=', equal + 1) != std::string::npos) {
        throw ManifestError(
            location,
            "copy-number must use TAXON=N or FAMILY:TAXON=N syntax; observed '" +
                value + "'");
    }

    const std::string left = value.substr(0, equal);
    const std::string countText = value.substr(equal + 1);
    const std::size_t colon = left.find(':');
    CopyCapacity capacity;
    capacity.sourceLocation = location;
    if (colon == std::string::npos) {
        capacity.key.taxonId = left;
    } else {
        if (colon == 0 || colon + 1 >= left.size() ||
            left.find(':', colon + 1) != std::string::npos) {
            throw ManifestError(
                location,
                "copy-number must use TAXON=N or FAMILY:TAXON=N syntax; observed '" +
                    value + "'");
        }
        capacity.key.familyId = left.substr(0, colon);
        capacity.key.taxonId = left.substr(colon + 1);
    }
    requireValue(capacity.key.taxonId, location, "copy-number taxon");
    capacity.count = detail::parseInt64(countText, location, "copy-number count");
    if (capacity.count < 0) {
        throw ManifestError(location, "copy-number count must be non-negative");
    }
    return capacity;
}

}  // namespace trio
}  // namespace anchorwave
