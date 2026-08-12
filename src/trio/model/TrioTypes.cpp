#include "TrioTypes.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

namespace anchorwave {
namespace trio {
namespace {

std::string formatManifestError(const SourceLocation &location,
                                const std::string &message) {
    std::ostringstream stream;
    stream << (location.path.empty() ? std::string("<manifest>") : location.path);
    if (location.line != 0) {
        stream << ':' << location.line;
    }
    stream << ": " << message;
    return stream.str();
}

std::string joinHeader(const std::vector<std::string> &columns) {
    std::ostringstream stream;
    for (std::size_t index = 0; index < columns.size(); ++index) {
        if (index != 0) {
            stream << '\t';
        }
        stream << columns[index];
    }
    return stream.str();
}

}  // namespace

SourceLocation::SourceLocation() : line(0) {}

SourceLocation::SourceLocation(const std::string &pathValue, std::size_t lineValue)
    : path(pathValue), line(lineValue) {}

ManifestError::ManifestError(const SourceLocation &location,
                             const std::string &message)
    : std::runtime_error(formatManifestError(location, message)), location_(location) {}

const SourceLocation &ManifestError::location() const {
    return location_;
}

std::string toString(const TaxonRole role) {
    switch (role) {
        case TaxonRole::INGROUP_REFERENCE:
            return "ingroup_reference";
        case TaxonRole::INGROUP:
            return "ingroup";
        case TaxonRole::PRIMARY_OUTGROUP:
            return "primary_outgroup";
        case TaxonRole::OUTGROUP:
            return "outgroup";
    }
    return "unknown";
}

TaxonRole parseTaxonRole(const std::string &text,
                         const SourceLocation &location) {
    if (text == "ingroup_reference") {
        return TaxonRole::INGROUP_REFERENCE;
    }
    if (text == "ingroup") {
        return TaxonRole::INGROUP;
    }
    if (text == "primary_outgroup") {
        return TaxonRole::PRIMARY_OUTGROUP;
    }
    if (text == "outgroup") {
        return TaxonRole::OUTGROUP;
    }
    throw ManifestError(
        location,
        "invalid role '" + text +
            "'; expected ingroup_reference, ingroup, primary_outgroup, or outgroup");
}

GenomicInterval::GenomicInterval() : start0(0), end0(0), strand('+') {}

GenomicInterval::GenomicInterval(const std::string &sequenceValue,
                                 const std::int64_t startValue,
                                 const std::int64_t endValue,
                                 const char strandValue)
    : sequence(sequenceValue), start0(startValue), end0(endValue), strand(strandValue) {}

std::int64_t GenomicInterval::length() const {
    return end0 - start0;
}

void validateHalfOpenInterval(const GenomicInterval &interval,
                              const SourceLocation &location,
                              const std::string &fieldName) {
    if (interval.sequence.empty() || interval.sequence == ".") {
        throw ManifestError(location, fieldName + " sequence must not be empty or '.'");
    }
    if (interval.start0 < 0) {
        throw ManifestError(location, fieldName + " start0 must be non-negative");
    }
    if (interval.end0 <= interval.start0) {
        throw ManifestError(
            location,
            fieldName + " must be a non-empty 0-based half-open interval (end0 > start0)");
    }
    if (interval.strand != '+' && interval.strand != '-') {
        throw ManifestError(location, fieldName + " strand must be '+' or '-'");
    }
}

TaxonRecord::TaxonRecord()
    : role(TaxonRole::INGROUP), qualityWeight(1.0) {}

bool TaxonManifest::contains(const TaxonId &taxonId) const {
    return find(taxonId) != nullptr;
}

const TaxonRecord *TaxonManifest::find(const TaxonId &taxonId) const {
    for (const TaxonRecord &record : records) {
        if (record.taxonId == taxonId) {
            return &record;
        }
    }
    return nullptr;
}

std::vector<TaxonId> TaxonManifest::ingroupIds() const {
    std::vector<TaxonId> preferred;
    std::vector<TaxonId> remaining;
    for (const TaxonRecord &record : records) {
        if (record.role == TaxonRole::INGROUP_REFERENCE) {
            preferred.push_back(record.taxonId);
        } else if (record.role == TaxonRole::INGROUP) {
            remaining.push_back(record.taxonId);
        }
    }
    std::sort(preferred.begin(), preferred.end());
    std::sort(remaining.begin(), remaining.end());
    preferred.insert(preferred.end(), remaining.begin(), remaining.end());
    return preferred;
}

TaxonId TaxonManifest::primaryOutgroupId() const {
    for (const TaxonRecord &record : records) {
        if (record.role == TaxonRole::PRIMARY_OUTGROUP) {
            return record.taxonId;
        }
    }
    return TaxonId();
}

std::vector<TaxonId> TaxonManifest::outgroupIds() const {
    std::vector<TaxonId> primary;
    std::vector<TaxonId> additional;
    for (const TaxonRecord &record : records) {
        if (record.role == TaxonRole::PRIMARY_OUTGROUP) {
            primary.push_back(record.taxonId);
        } else if (record.role == TaxonRole::OUTGROUP) {
            additional.push_back(record.taxonId);
        }
    }
    std::sort(primary.begin(), primary.end());
    std::sort(additional.begin(), additional.end());
    primary.insert(primary.end(), additional.begin(), additional.end());
    return primary;
}

std::vector<TaxonId> TaxonManifest::allTaxonIds() const {
    std::vector<TaxonId> result;
    result.reserve(records.size());
    for (const TaxonRecord &record : records) {
        result.push_back(record.taxonId);
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::string toString(const PairwiseScope scope) {
    return scope == PairwiseScope::TRIANGLES ? "triangles" : "complete";
}

CanonicalTaxonPair::CanonicalTaxonPair() {}

CanonicalTaxonPair::CanonicalTaxonPair(const TaxonId &taxonA,
                                       const TaxonId &taxonB) {
    if (taxonA.empty() || taxonB.empty()) {
        throw std::invalid_argument("canonical taxon pair requires two non-empty IDs");
    }
    if (taxonA == taxonB) {
        throw std::invalid_argument("canonical taxon pair requires two distinct taxa");
    }
    if (taxonA < taxonB) {
        first = taxonA;
        second = taxonB;
    } else {
        first = taxonB;
        second = taxonA;
    }
}

std::string CanonicalTaxonPair::stableId() const {
    return first + "--" + second;
}

bool CanonicalTaxonPair::operator<(const CanonicalTaxonPair &other) const {
    if (first != other.first) {
        return first < other.first;
    }
    return second < other.second;
}

bool CanonicalTaxonPair::operator==(const CanonicalTaxonPair &other) const {
    return first == other.first && second == other.second;
}

bool CanonicalTaxonPair::operator!=(const CanonicalTaxonPair &other) const {
    return !(*this == other);
}

PairwiseManifestRecord::PairwiseManifestRecord() : weight(1.0) {}

bool PairwiseManifest::contains(const CanonicalTaxonPair &pair) const {
    return find(pair) != nullptr;
}

const PairwiseManifestRecord *PairwiseManifest::find(
    const CanonicalTaxonPair &pair) const {
    for (const PairwiseManifestRecord &record : records) {
        if (record.pair == pair) {
            return &record;
        }
    }
    return nullptr;
}

std::string toString(const CopyRecordType type) {
    switch (type) {
        case CopyRecordType::GROUP:
            return "GROUP";
        case CopyRecordType::MEMBER:
            return "MEMBER";
        case CopyRecordType::COUNT:
            return "COUNT";
    }
    return "UNKNOWN";
}

std::string toString(const CopyMemberType type) {
    switch (type) {
        case CopyMemberType::NONE:
            return ".";
        case CopyMemberType::ANCHOR:
            return "anchor";
        case CopyMemberType::INTERVAL:
            return "interval";
    }
    return "unknown";
}

std::string toString(const CopyRelation relation) {
    switch (relation) {
        case CopyRelation::UNSPECIFIED:
            return ".";
        case CopyRelation::ORTHOLOG:
            return "ortholog";
        case CopyRelation::COORTHOLOG:
            return "coortholog";
        case CopyRelation::EXCLUDED:
            return "excluded";
    }
    return "unknown";
}

std::string toString(const ConstraintStrength strength) {
    return strength == ConstraintStrength::HARD ? "hard" : "soft";
}

CopyConstraint::CopyConstraint()
    : recordType(CopyRecordType::GROUP),
      memberType(CopyMemberType::NONE),
      hasInterval(false),
      hasExpectedCount(false),
      expectedCount(0),
      relation(CopyRelation::UNSPECIFIED),
      confidence(1.0),
      strength(ConstraintStrength::HARD) {}

bool CopyCapacityKey::operator<(const CopyCapacityKey &other) const {
    if (familyId != other.familyId) {
        return familyId < other.familyId;
    }
    return taxonId < other.taxonId;
}

bool CopyCapacityKey::operator==(const CopyCapacityKey &other) const {
    return familyId == other.familyId && taxonId == other.taxonId;
}

CopyCapacity::CopyCapacity() : count(0) {}

bool CopyConstraintSet::findDefaultCapacity(const FamilyId &familyId,
                                            const TaxonId &taxonId,
                                            std::int64_t &count) const {
    CopyCapacityKey familyKey;
    familyKey.familyId = familyId;
    familyKey.taxonId = taxonId;
    CopyCapacityMap::const_iterator found = defaultCapacities.find(familyKey);
    if (found != defaultCapacities.end()) {
        count = found->second.count;
        return true;
    }

    CopyCapacityKey taxonKey;
    taxonKey.taxonId = taxonId;
    found = defaultCapacities.find(taxonKey);
    if (found == defaultCapacities.end()) {
        return false;
    }
    count = found->second.count;
    return true;
}

namespace detail {

namespace {

bool isAbsolutePath(const std::string &path) {
    return !path.empty() && path[0] == '/';
}

std::string currentWorkingDirectory(const SourceLocation &location) {
    std::vector<char> buffer(256, '\0');
    for (;;) {
        errno = 0;
        if (getcwd(buffer.data(), buffer.size()) != nullptr) {
            return std::string(buffer.data());
        }
        if (errno != ERANGE) {
            throw ManifestError(location,
                                "cannot determine the current working directory "
                                "while resolving a manifest-relative path");
        }
        if (buffer.size() > std::numeric_limits<std::size_t>::max() / 2) {
            throw ManifestError(location,
                                "current working directory path is too long");
        }
        buffer.resize(buffer.size() * 2, '\0');
    }
}

std::string parentDirectory(const std::string &path) {
    const std::size_t separator = path.find_last_of('/');
    if (separator == std::string::npos) {
        return std::string(".");
    }
    if (separator == 0) {
        return std::string("/");
    }
    return path.substr(0, separator);
}

std::string joinPath(const std::string &directory, const std::string &entry) {
    if (directory.empty() || directory == "/") {
        return directory + entry;
    }
    return directory + "/" + entry;
}

}  // namespace

std::string trimAscii(const std::string &text) {
    std::size_t begin = 0;
    while (begin < text.size() &&
           (text[begin] == ' ' || text[begin] == '\t' || text[begin] == '\r' ||
            text[begin] == '\n')) {
        ++begin;
    }
    std::size_t end = text.size();
    while (end > begin &&
           (text[end - 1] == ' ' || text[end - 1] == '\t' || text[end - 1] == '\r' ||
            text[end - 1] == '\n')) {
        --end;
    }
    return text.substr(begin, end - begin);
}

bool isMissingValue(const std::string &text) {
    const std::string trimmed = trimAscii(text);
    return trimmed.empty() || trimmed == ".";
}

bool isBlankOrComment(const std::string &line) {
    const std::string trimmed = trimAscii(line);
    return trimmed.empty() || trimmed[0] == '#';
}

std::vector<std::string> splitTsv(const std::string &line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.push_back(trimAscii(line.substr(begin)));
            break;
        }
        fields.push_back(trimAscii(line.substr(begin, tab - begin)));
        begin = tab + 1;
    }
    return fields;
}

void requireExactHeader(const std::vector<std::string> &observed,
                        const std::vector<std::string> &expected,
                        const SourceLocation &location) {
    if (observed != expected) {
        throw ManifestError(location,
                            "invalid TSV header; expected exactly: " + joinHeader(expected));
    }
}

std::string resolveManifestEntryPath(const std::string &manifestPath,
                                     const std::string &entryPath,
                                     const SourceLocation &location) {
    if (entryPath == "." || isAbsolutePath(entryPath)) {
        return entryPath;
    }

    std::string absoluteManifestPath = manifestPath;
    if (!isAbsolutePath(absoluteManifestPath)) {
        absoluteManifestPath =
            joinPath(currentWorkingDirectory(location), absoluteManifestPath);
    }
    return joinPath(parentDirectory(absoluteManifestPath), entryPath);
}

std::int64_t parseInt64(const std::string &text,
                        const SourceLocation &location,
                        const std::string &fieldName) {
    const std::string value = trimAscii(text);
    if (value.empty()) {
        throw ManifestError(location, fieldName + " must be an integer");
    }
    errno = 0;
    char *end = nullptr;
    const long long parsed = std::strtoll(value.c_str(), &end, 10);
    if (errno == ERANGE || end == value.c_str() || *end != '\0') {
        throw ManifestError(location,
                            fieldName + " must be a signed 64-bit integer; observed '" +
                                value + "'");
    }
    if (parsed < std::numeric_limits<std::int64_t>::min() ||
        parsed > std::numeric_limits<std::int64_t>::max()) {
        throw ManifestError(location, fieldName + " is outside signed 64-bit range");
    }
    return static_cast<std::int64_t>(parsed);
}

double parseFiniteDouble(const std::string &text,
                         const SourceLocation &location,
                         const std::string &fieldName) {
    const std::string value = trimAscii(text);
    if (value.empty()) {
        throw ManifestError(location, fieldName + " must be numeric");
    }
    errno = 0;
    char *end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (errno == ERANGE || end == value.c_str() || *end != '\0' ||
        !std::isfinite(parsed)) {
        throw ManifestError(location,
                            fieldName + " must be a finite number; observed '" + value +
                                "'");
    }
    return parsed;
}

void validateExistingPath(const std::string &path,
                          const SourceLocation &location,
                          const std::string &fieldName) {
    struct stat status;
    if (path.empty() || path == "." || stat(path.c_str(), &status) != 0) {
        throw ManifestError(location,
                            fieldName + " path does not exist: '" + path + "'");
    }
}

}  // namespace detail

}  // namespace trio
}  // namespace anchorwave
