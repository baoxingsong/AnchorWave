#include "src/trio/io/FastaEvidenceValidator.h"

#include <zlib.h>

#include <algorithm>
#include <cctype>
#include <climits>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace anchorwave {
namespace trio {
namespace {

struct RequiredObservation {
    int64_t position0 = -1;
    char base = 'N';
    std::string identity;
};

struct SourceRequirement {
    bool hasExpectedSize = false;
    int64_t expectedSize = 0;
    std::vector<RequiredObservation> observations;
};

using TaxonRequirements = std::map<std::string, SourceRequirement>;
using Requirements = std::map<std::string, TaxonRequirements>;

unsigned int nucleotideMask(char value, bool allowEvidenceMissingSymbols) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(value)))) {
        case 'A': return 0x1u;
        case 'C': return 0x2u;
        case 'G': return 0x4u;
        case 'T': case 'U': return 0x8u;
        case 'R': return 0x5u;
        case 'Y': return 0xAu;
        case 'S': return 0x6u;
        case 'W': return 0x9u;
        case 'K': return 0xCu;
        case 'M': return 0x3u;
        case 'B': return 0xEu;
        case 'D': return 0xDu;
        case 'H': return 0xBu;
        case 'V': return 0x7u;
        case 'N': case 'X': return 0xFu;
        case '.': case '*': case '?':
            return allowEvidenceMissingSymbols ? 0xFu : 0u;
        default: return 0u;
    }
}

std::string quotedChar(char value) {
    const unsigned char byte = static_cast<unsigned char>(value);
    if (std::isprint(byte)) {
        return std::string("'") + value + "'";
    }
    static const char digits[] = "0123456789ABCDEF";
    std::string result("0x00");
    result[2] = digits[(byte >> 4) & 0xFu];
    result[3] = digits[byte & 0xFu];
    return result;
}

void addDiagnostic(std::vector<FastaEvidenceDiagnostic> &diagnostics,
                   const std::string &code, const std::string &taxon,
                   const std::string &fastaPath, const std::string &sequence,
                   int64_t position0, std::size_t line,
                   const std::string &message) {
    FastaEvidenceDiagnostic diagnostic;
    diagnostic.code = code;
    diagnostic.taxon = taxon;
    diagnostic.fastaPath = fastaPath;
    diagnostic.sequence = sequence;
    diagnostic.position0 = position0;
    diagnostic.line = line;
    diagnostic.message = message;
    diagnostics.push_back(diagnostic);
}

bool sameResidueId(const ResidueId &left, const ResidueId &right) {
    return left == right;
}

Requirements deriveRequirements(
    const AlignmentEvidenceSet &evidence,
    FastaEvidenceValidationReport &report) {
    Requirements requirements;
    for (const auto &entry : evidence.sourceSizes) {
        const std::string &taxon = entry.first.first;
        const std::string &sequence = entry.first.second;
        if (taxon.empty() || sequence.empty()) {
            addDiagnostic(report.diagnostics, "EVIDENCE_EMPTY_SOURCE_ID", taxon,
                          "", sequence, -1, 0,
                          "declared source size has an empty taxon or sequence ID");
            continue;
        }
        SourceRequirement &requirement = requirements[taxon][sequence];
        if (entry.second < 0) {
            addDiagnostic(report.diagnostics, "EVIDENCE_INVALID_SOURCE_SIZE",
                          taxon, "", sequence, -1, 0,
                          "declared source size is negative");
            continue;
        }
        requirement.hasExpectedSize = true;
        requirement.expectedSize = entry.second;
    }

    for (const auto &entry : evidence.residues) {
        ++report.counts.residueObservationsRequired;
        const ResidueId &key = entry.first;
        const ResidueObservation &observation = entry.second;
        if (!sameResidueId(key, observation.id)) {
            addDiagnostic(report.diagnostics, "EVIDENCE_RESIDUE_KEY_MISMATCH",
                          key.taxon, "", key.sequence, key.forwardPosition0, 0,
                          "residue map key differs from the stored observation ID");
        }
        if (key.taxon.empty() || key.sequence.empty()) {
            addDiagnostic(report.diagnostics, "EVIDENCE_EMPTY_SOURCE_ID",
                          key.taxon, "", key.sequence, key.forwardPosition0, 0,
                          "residue has an empty taxon or sequence ID");
            continue;
        }

        SourceRequirement &requirement = requirements[key.taxon][key.sequence];
        if (observation.sourceSize < 0) {
            addDiagnostic(report.diagnostics, "EVIDENCE_INVALID_SOURCE_SIZE",
                          key.taxon, "", key.sequence, key.forwardPosition0, 0,
                          "residue observation has a negative source size");
        } else if (!requirement.hasExpectedSize) {
            requirement.hasExpectedSize = true;
            requirement.expectedSize = observation.sourceSize;
        } else if (requirement.expectedSize != observation.sourceSize) {
            std::ostringstream message;
            message << "residue source size " << observation.sourceSize
                    << " differs from declared size " << requirement.expectedSize;
            addDiagnostic(report.diagnostics, "EVIDENCE_SOURCE_SIZE_CONFLICT",
                          key.taxon, "", key.sequence, key.forwardPosition0, 0,
                          message.str());
        }
        if (key.forwardPosition0 < 0 ||
            (requirement.hasExpectedSize &&
             key.forwardPosition0 >= requirement.expectedSize)) {
            addDiagnostic(report.diagnostics, "EVIDENCE_POSITION_OUT_OF_RANGE",
                          key.taxon, "", key.sequence, key.forwardPosition0, 0,
                          "residue position is outside the declared source interval");
            continue;
        }
        if (nucleotideMask(observation.base, true) == 0u) {
            addDiagnostic(report.diagnostics, "EVIDENCE_INVALID_BASE",
                          key.taxon, "", key.sequence, key.forwardPosition0, 0,
                          "residue uses unsupported symbol " +
                              quotedChar(observation.base));
            continue;
        }
        RequiredObservation required;
        required.position0 = key.forwardPosition0;
        required.base = observation.base;
        required.identity = key.canonicalString();
        requirement.observations.push_back(required);
    }

    for (auto &taxon : requirements) {
        for (auto &source : taxon.second) {
            std::sort(source.second.observations.begin(),
                      source.second.observations.end(),
                      [](const RequiredObservation &left,
                         const RequiredObservation &right) {
                          return std::tie(left.position0, left.identity, left.base) <
                                 std::tie(right.position0, right.identity, right.base);
                      });
        }
    }
    report.counts.taxaRequired = requirements.size();
    for (const auto &taxon : requirements) {
        report.counts.sourcesRequired += taxon.second.size();
    }
    return requirements;
}

class FastaStreamInspector {
public:
    FastaStreamInspector(const std::string &taxon, const std::string &path,
                         const TaxonRequirements &requirements,
                         FastaEvidenceValidationReport &report)
        : taxon_(taxon), path_(path), requirements_(requirements),
          report_(report) {}

    void run() {
        gzFile input = gzopen(path_.c_str(), "rb");
        if (input == nullptr) {
            addDiagnostic(report_.diagnostics, "FASTA_OPEN_ERROR", taxon_, path_,
                          "", -1, 0, "unable to open FASTA");
            markMissingSources();
            return;
        }
        ++report_.counts.fastaFilesScanned;
        (void)gzbuffer(input, 1U << 20);

        char buffer[1U << 16];
        bool ioFailed = false;
        while (true) {
            const int count = gzread(input, buffer, static_cast<unsigned int>(sizeof(buffer)));
            if (count > 0) {
                for (int i = 0; i < count; ++i) consume(buffer[i]);
                continue;
            }
            if (count < 0) {
                int errorNumber = Z_OK;
                const char *errorText = gzerror(input, &errorNumber);
                addDiagnostic(report_.diagnostics, "FASTA_IO_ERROR", taxon_, path_,
                              "", -1, line_,
                              std::string("error while reading FASTA: ") +
                                  (errorText == nullptr ? "unknown zlib error" : errorText));
                ioFailed = true;
            }
            break;
        }
        finishInput();
        const int closeStatus = gzclose(input);
        if (closeStatus != Z_OK && !ioFailed) {
            addDiagnostic(report_.diagnostics, "FASTA_IO_ERROR", taxon_, path_,
                          "", -1, line_, "error while closing FASTA stream");
        }
        markMissingSources();
        countValidatedSources();
    }

private:
    enum class HeaderStage { LEADING_SPACE, ID, DESCRIPTION };

    struct RequiredSourceStatus {
        std::size_t occurrences = 0;
        bool firstOccurrenceChecksPassed = false;
    };

    void consume(char value) {
        if (value == '\r') return;  // supports CRLF without spelling CR as a base
        if (value == '\n') {
            if (inHeader_) commitHeader();
            ++line_;
            atLineStart_ = true;
            return;
        }
        if (atLineStart_) {
            atLineStart_ = false;
            if (value == '>') {
                finalizeRecord();
                beginHeader();
                return;
            }
        }
        if (inHeader_) {
            consumeHeader(value);
        } else {
            consumeSequence(value);
        }
    }

    void beginHeader() {
        inHeader_ = true;
        headerStage_ = HeaderStage::LEADING_SPACE;
        headerId_.clear();
        headerHasInvalidControl_ = false;
        headerLine_ = line_;
    }

    void consumeHeader(char value) {
        const unsigned char byte = static_cast<unsigned char>(value);
        if (std::iscntrl(byte)) headerHasInvalidControl_ = true;
        if (headerStage_ == HeaderStage::DESCRIPTION) return;
        if (std::isspace(byte)) {
            if (headerStage_ == HeaderStage::ID) {
                headerStage_ = HeaderStage::DESCRIPTION;
            }
            return;
        }
        if (headerStage_ == HeaderStage::LEADING_SPACE) {
            headerStage_ = HeaderStage::ID;
        }
        headerId_.push_back(value);
    }

    void commitHeader() {
        inHeader_ = false;
        haveRecord_ = true;
        currentSequence_ = headerId_;
        currentLength_ = 0;
        firstInvalidBaseCount_ = 0;
        firstInvalidBasePosition0_ = -1;
        firstInvalidBaseLine_ = 0;
        firstInvalidBase_ = '\0';
        mismatchCount_ = 0;
        firstMismatchPosition0_ = -1;
        firstMismatchFastaBase_ = '\0';
        firstMismatchEvidenceBase_ = '\0';
        currentRequired_ = false;
        compareCurrentRecord_ = false;
        currentObservationIndex_ = 0;

        ++report_.counts.sequencesScanned;
        if (headerId_.empty()) {
            haveNamedRecord_ = false;
            addDiagnostic(report_.diagnostics, "FASTA_EMPTY_SEQUENCE_ID", taxon_,
                          path_, "", -1, headerLine_,
                          "FASTA header has no sequence ID");
            return;
        }
        haveNamedRecord_ = true;
        if (headerHasInvalidControl_) {
            addDiagnostic(report_.diagnostics, "FASTA_INVALID_HEADER", taxon_,
                          path_, headerId_, -1, headerLine_,
                          "FASTA header contains an unsupported control character");
        }

        const auto inserted = firstHeaderLine_.insert(
            std::make_pair(headerId_, headerLine_));
        const auto required = requirements_.find(headerId_);
        currentRequired_ = required != requirements_.end();
        if (!inserted.second) {
            std::ostringstream message;
            message << "duplicate sequence ID; first header is at line "
                    << inserted.first->second;
            addDiagnostic(report_.diagnostics, "FASTA_DUPLICATE_SEQUENCE_ID",
                          taxon_, path_, headerId_, -1, headerLine_, message.str());
        }
        if (currentRequired_) {
            RequiredSourceStatus &status = sourceStatus_[headerId_];
            ++status.occurrences;
            if (status.occurrences == 1) {
                ++report_.counts.sourcesFound;
                compareCurrentRecord_ = true;
            }
        }
    }

    void consumeSequence(char value) {
        const unsigned char byte = static_cast<unsigned char>(value);
        if (std::isspace(byte)) return;
        if (!haveRecord_) {
            if (!reportedSequenceBeforeHeader_) {
                addDiagnostic(report_.diagnostics, "FASTA_SEQUENCE_BEFORE_HEADER",
                              taxon_, path_, "", -1, line_,
                              "sequence content appears before the first FASTA header");
                reportedSequenceBeforeHeader_ = true;
            }
            return;
        }

        const int64_t position0 =
            currentLength_ > static_cast<std::uint64_t>(
                                 std::numeric_limits<int64_t>::max())
                ? std::numeric_limits<int64_t>::max()
                : static_cast<int64_t>(currentLength_);
        ++currentLength_;
        ++report_.counts.basesScanned;
        const unsigned int fastaMask = nucleotideMask(value, false);
        if (fastaMask == 0u) {
            if (firstInvalidBaseCount_ == 0) {
                firstInvalidBasePosition0_ = position0;
                firstInvalidBaseLine_ = line_;
                firstInvalidBase_ = value;
            }
            ++firstInvalidBaseCount_;
        }

        if (!compareCurrentRecord_) return;
        const std::vector<RequiredObservation> &observations =
            requirements_.at(currentSequence_).observations;
        while (currentObservationIndex_ < observations.size() &&
               observations[currentObservationIndex_].position0 < position0) {
            ++currentObservationIndex_;
        }
        while (currentObservationIndex_ < observations.size() &&
               observations[currentObservationIndex_].position0 == position0) {
            const RequiredObservation &observation =
                observations[currentObservationIndex_];
            const unsigned int evidenceMask = nucleotideMask(observation.base, true);
            if (fastaMask != 0u && (fastaMask & evidenceMask) != 0u) {
                ++report_.counts.residueObservationsValidated;
            } else if (fastaMask != 0u) {
                if (mismatchCount_ == 0) {
                    firstMismatchPosition0_ = position0;
                    firstMismatchFastaBase_ = value;
                    firstMismatchEvidenceBase_ = observation.base;
                }
                ++mismatchCount_;
            }
            ++currentObservationIndex_;
        }
    }

    void finalizeRecord() {
        if (!haveRecord_) return;
        if (firstInvalidBaseCount_ != 0) {
            std::ostringstream message;
            message << "sequence contains " << firstInvalidBaseCount_
                    << " non-IUPAC symbol(s); first is "
                    << quotedChar(firstInvalidBase_);
            addDiagnostic(report_.diagnostics, "FASTA_INVALID_BASE", taxon_, path_,
                          haveNamedRecord_ ? currentSequence_ : "",
                          firstInvalidBasePosition0_, firstInvalidBaseLine_,
                          message.str());
        }
        if (mismatchCount_ != 0) {
            std::ostringstream message;
            message << mismatchCount_
                    << " evidence observation(s) are incompatible with the forward "
                       "FASTA spelling; first has evidence "
                    << quotedChar(firstMismatchEvidenceBase_) << " and FASTA "
                    << quotedChar(firstMismatchFastaBase_);
            addDiagnostic(report_.diagnostics, "FASTA_EVIDENCE_BASE_MISMATCH",
                          taxon_, path_, currentSequence_,
                          firstMismatchPosition0_, 0, message.str());
        }
        if (compareCurrentRecord_) {
            const SourceRequirement &requirement =
                requirements_.at(currentSequence_);
            bool passed = firstInvalidBaseCount_ == 0 && mismatchCount_ == 0;
            if (currentLength_ > static_cast<std::uint64_t>(
                                     std::numeric_limits<int64_t>::max())) {
                addDiagnostic(report_.diagnostics, "FASTA_SOURCE_TOO_LONG", taxon_,
                              path_, currentSequence_, -1, headerLine_,
                              "FASTA source length exceeds signed 64-bit coordinates");
                passed = false;
            } else if (requirement.hasExpectedSize &&
                       static_cast<int64_t>(currentLength_) !=
                           requirement.expectedSize) {
                std::ostringstream message;
                message << "FASTA length " << currentLength_
                        << " differs from evidence source size "
                        << requirement.expectedSize;
                addDiagnostic(report_.diagnostics, "FASTA_SOURCE_SIZE_MISMATCH",
                              taxon_, path_, currentSequence_, -1, headerLine_,
                              message.str());
                passed = false;
            }
            sourceStatus_[currentSequence_].firstOccurrenceChecksPassed = passed;
        }
        haveRecord_ = false;
        haveNamedRecord_ = false;
    }

    void finishInput() {
        if (inHeader_) commitHeader();
        finalizeRecord();
    }

    void markMissingSources() {
        for (const auto &requirement : requirements_) {
            const auto status = sourceStatus_.find(requirement.first);
            if (status == sourceStatus_.end() || status->second.occurrences == 0) {
                addDiagnostic(report_.diagnostics, "FASTA_MISSING_SEQUENCE_ID",
                              taxon_, path_, requirement.first, -1, 0,
                              "sequence required by alignment evidence is absent from FASTA");
            }
        }
    }

    void countValidatedSources() {
        for (const auto &status : sourceStatus_) {
            if (status.second.occurrences == 1 &&
                status.second.firstOccurrenceChecksPassed) {
                ++report_.counts.sourcesValidated;
            }
        }
    }

    const std::string &taxon_;
    const std::string &path_;
    const TaxonRequirements &requirements_;
    FastaEvidenceValidationReport &report_;
    std::map<std::string, std::size_t> firstHeaderLine_;
    std::map<std::string, RequiredSourceStatus> sourceStatus_;

    std::size_t line_ = 1;
    bool atLineStart_ = true;
    bool inHeader_ = false;
    HeaderStage headerStage_ = HeaderStage::LEADING_SPACE;
    std::string headerId_;
    bool headerHasInvalidControl_ = false;
    std::size_t headerLine_ = 0;
    bool haveRecord_ = false;
    bool haveNamedRecord_ = false;
    std::string currentSequence_;
    std::uint64_t currentLength_ = 0;
    bool currentRequired_ = false;
    bool compareCurrentRecord_ = false;
    std::size_t currentObservationIndex_ = 0;
    bool reportedSequenceBeforeHeader_ = false;

    std::uint64_t firstInvalidBaseCount_ = 0;
    int64_t firstInvalidBasePosition0_ = -1;
    std::size_t firstInvalidBaseLine_ = 0;
    char firstInvalidBase_ = '\0';
    std::uint64_t mismatchCount_ = 0;
    int64_t firstMismatchPosition0_ = -1;
    char firstMismatchFastaBase_ = '\0';
    char firstMismatchEvidenceBase_ = '\0';
};

std::string errorMessage(const FastaEvidenceValidationReport &report) {
    std::ostringstream message;
    message << "FASTA evidence validation failed with "
            << report.diagnostics.size() << " diagnostic(s)";
    for (const FastaEvidenceDiagnostic &diagnostic : report.diagnostics) {
        message << '\n' << diagnostic.canonicalString();
    }
    return message.str();
}

}  // namespace

bool FastaEvidenceDiagnostic::operator<(
    const FastaEvidenceDiagnostic &other) const {
    return std::tie(taxon, code, sequence, position0, line, fastaPath, message) <
           std::tie(other.taxon, other.code, other.sequence, other.position0,
                    other.line, other.fastaPath, other.message);
}

std::string FastaEvidenceDiagnostic::canonicalString() const {
    std::ostringstream output;
    output << '[' << code << ']';
    if (!taxon.empty()) output << " taxon=" << taxon;
    if (!fastaPath.empty()) output << " fasta=" << fastaPath;
    if (!sequence.empty()) output << " sequence=" << sequence;
    if (position0 >= 0) output << " position0=" << position0;
    if (line != 0) output << " line=" << line;
    output << ": " << message;
    return output.str();
}

FastaEvidenceValidationError::FastaEvidenceValidationError(
    const FastaEvidenceValidationReport &report)
    : std::runtime_error(errorMessage(report)), report_(report) {}

FastaEvidenceValidationReport FastaEvidenceValidator::inspect(
    const AlignmentEvidenceSet &evidence,
    const TaxonFastaMap &fastaByTaxon) {
    FastaEvidenceValidationReport report;
    const Requirements requirements = deriveRequirements(evidence, report);
    for (const auto &taxon : requirements) {
        const auto fasta = fastaByTaxon.find(taxon.first);
        if (fasta == fastaByTaxon.end() || fasta->second.empty()) {
            addDiagnostic(report.diagnostics, "MISSING_TAXON_FASTA", taxon.first,
                          "", "", -1, 0,
                          "no FASTA path was supplied for evidence taxon");
            for (const auto &source : taxon.second) {
                addDiagnostic(report.diagnostics, "FASTA_MISSING_SEQUENCE_ID",
                              taxon.first, "", source.first, -1, 0,
                              "sequence cannot be checked because taxon FASTA is missing");
            }
            continue;
        }
        FastaStreamInspector inspector(taxon.first, fasta->second, taxon.second,
                                       report);
        inspector.run();
    }
    std::sort(report.diagnostics.begin(), report.diagnostics.end());
    return report;
}

FastaEvidenceValidationCounts FastaEvidenceValidator::validate(
    const AlignmentEvidenceSet &evidence,
    const TaxonFastaMap &fastaByTaxon) {
    const FastaEvidenceValidationReport report = inspect(evidence, fastaByTaxon);
    if (!report.valid()) throw FastaEvidenceValidationError(report);
    return report.counts;
}

}  // namespace trio
}  // namespace anchorwave
