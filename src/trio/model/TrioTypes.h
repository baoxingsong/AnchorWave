#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

using TaxonId = std::string;
using FamilyId = std::string;
using AncestralCopyId = std::string;
using OccurrencePathId = std::string;

struct SourceLocation {
    std::string path;
    std::size_t line;

    SourceLocation();
    SourceLocation(const std::string &pathValue, std::size_t lineValue);
};

class ManifestError : public std::runtime_error {
public:
    ManifestError(const SourceLocation &location, const std::string &message);

    const SourceLocation &location() const;

private:
    SourceLocation location_;
};

enum class TaxonRole {
    INGROUP_REFERENCE,
    INGROUP,
    PRIMARY_OUTGROUP,
    OUTGROUP
};

std::string toString(TaxonRole role);
TaxonRole parseTaxonRole(const std::string &text, const SourceLocation &location);

struct GenomicInterval {
    std::string sequence;
    std::int64_t start0;
    std::int64_t end0;
    char strand;

    GenomicInterval();
    GenomicInterval(const std::string &sequenceValue,
                    std::int64_t startValue,
                    std::int64_t endValue,
                    char strandValue);

    std::int64_t length() const;
};

void validateHalfOpenInterval(const GenomicInterval &interval,
                              const SourceLocation &location,
                              const std::string &fieldName);

struct TaxonRecord {
    TaxonId taxonId;
    TaxonRole role;
    std::string fasta;
    std::string gff;
    std::string anchorSam;
    std::string anchorFasta;
    std::string callabilityBed;
    double qualityWeight;
    SourceLocation sourceLocation;

    TaxonRecord();
};

struct TaxonManifest {
    std::vector<TaxonRecord> records;

    bool contains(const TaxonId &taxonId) const;
    const TaxonRecord *find(const TaxonId &taxonId) const;
    std::vector<TaxonId> ingroupIds() const;
    TaxonId primaryOutgroupId() const;
    std::vector<TaxonId> outgroupIds() const;
    std::vector<TaxonId> allTaxonIds() const;
};

enum class PairwiseScope {
    TRIANGLES,
    COMPLETE
};

std::string toString(PairwiseScope scope);

struct CanonicalTaxonPair {
    TaxonId first;
    TaxonId second;

    CanonicalTaxonPair();
    CanonicalTaxonPair(const TaxonId &taxonA, const TaxonId &taxonB);

    std::string stableId() const;
    bool operator<(const CanonicalTaxonPair &other) const;
    bool operator==(const CanonicalTaxonPair &other) const;
    bool operator!=(const CanonicalTaxonPair &other) const;
};

struct PairwiseManifestRecord {
    CanonicalTaxonPair pair;
    TaxonId sourceTaxonA;
    TaxonId sourceTaxonB;
    std::string maf;
    std::string anchorMap;
    std::string scoreProfile;
    double weight;
    SourceLocation sourceLocation;

    PairwiseManifestRecord();
};

struct PairwiseManifest {
    std::vector<PairwiseManifestRecord> records;

    bool contains(const CanonicalTaxonPair &pair) const;
    const PairwiseManifestRecord *find(const CanonicalTaxonPair &pair) const;
};

enum class CopyRecordType {
    GROUP,
    MEMBER,
    COUNT
};

enum class CopyMemberType {
    NONE,
    ANCHOR,
    INTERVAL
};

enum class CopyRelation {
    UNSPECIFIED,
    ORTHOLOG,
    COORTHOLOG,
    EXCLUDED
};

enum class ConstraintStrength {
    HARD,
    SOFT
};

std::string toString(CopyRecordType type);
std::string toString(CopyMemberType type);
std::string toString(CopyRelation relation);
std::string toString(ConstraintStrength strength);

struct CopyConstraint {
    CopyRecordType recordType;
    FamilyId familyId;
    AncestralCopyId ancestralCopyId;
    std::string parentCopyId;
    std::string bornOnEdge;
    TaxonId taxonId;
    CopyMemberType memberType;
    std::string memberId;
    bool hasInterval;
    GenomicInterval interval;
    bool hasExpectedCount;
    std::int64_t expectedCount;
    CopyRelation relation;
    double confidence;
    ConstraintStrength strength;
    std::string source;
    SourceLocation sourceLocation;

    CopyConstraint();
};

struct CopyCapacityKey {
    FamilyId familyId;  // Empty means a taxon-wide fallback capacity.
    TaxonId taxonId;

    bool operator<(const CopyCapacityKey &other) const;
    bool operator==(const CopyCapacityKey &other) const;
};

struct CopyCapacity {
    CopyCapacityKey key;
    std::int64_t count;
    SourceLocation sourceLocation;

    CopyCapacity();
};

using CopyCapacityMap = std::map<CopyCapacityKey, CopyCapacity>;

struct CopyConstraintSet {
    std::vector<CopyConstraint> records;
    CopyCapacityMap defaultCapacities;

    bool findDefaultCapacity(const FamilyId &familyId,
                             const TaxonId &taxonId,
                             std::int64_t &count) const;
};

namespace detail {

std::string trimAscii(const std::string &text);
bool isMissingValue(const std::string &text);
bool isBlankOrComment(const std::string &line);
std::vector<std::string> splitTsv(const std::string &line);
void requireExactHeader(const std::vector<std::string> &observed,
                        const std::vector<std::string> &expected,
                        const SourceLocation &location);
// Resolve a path stored in a manifest against the directory containing that
// manifest.  Relative results are made absolute so downstream consumers do not
// accidentally resolve them a second time against the process working directory.
// The sentinel "." and already-absolute paths are returned unchanged.
std::string resolveManifestEntryPath(const std::string &manifestPath,
                                     const std::string &entryPath,
                                     const SourceLocation &location);
std::int64_t parseInt64(const std::string &text,
                        const SourceLocation &location,
                        const std::string &fieldName);
double parseFiniteDouble(const std::string &text,
                         const SourceLocation &location,
                         const std::string &fieldName);
void validateExistingPath(const std::string &path,
                          const SourceLocation &location,
                          const std::string &fieldName);

}  // namespace detail

}  // namespace trio
}  // namespace anchorwave
