#include "src/trio/io/FastaEvidenceValidator.h"

#include "src/trio/io/MafEvidenceReader.h"

#include "gtest/gtest.h"

#include <zlib.h>

#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {
namespace {

class TempFastaDirectory {
public:
    TempFastaDirectory() {
        char pattern[] = "/tmp/anchorwave-fasta-validator-XXXXXX";
        char *created = mkdtemp(pattern);
        if (created == nullptr) throw std::runtime_error("mkdtemp failed");
        path_ = created;
    }

    ~TempFastaDirectory() {
        for (std::vector<std::string>::const_reverse_iterator file = files_.rbegin();
             file != files_.rend(); ++file) {
            unlink(file->c_str());
        }
        rmdir(path_.c_str());
    }

    std::string plain(const std::string &name, const std::string &contents) {
        const std::string path = path_ + "/" + name;
        std::ofstream output(path.c_str(), std::ios::binary);
        if (!output.good()) throw std::runtime_error("cannot write plain FASTA fixture");
        output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
        if (!output.good()) throw std::runtime_error("cannot finish plain FASTA fixture");
        output.close();
        files_.push_back(path);
        return path;
    }

    std::string gzip(const std::string &name, const std::string &contents) {
        const std::string path = path_ + "/" + name;
        gzFile output = gzopen(path.c_str(), "wb");
        if (output == nullptr) throw std::runtime_error("cannot create gzip fixture");
        const int written = gzwrite(output, contents.data(),
                                    static_cast<unsigned int>(contents.size()));
        const int closed = gzclose(output);
        if (written != static_cast<int>(contents.size()) || closed != Z_OK) {
            throw std::runtime_error("cannot finish gzip fixture");
        }
        files_.push_back(path);
        return path;
    }

private:
    std::string path_;
    std::vector<std::string> files_;
};

ResidueId residueId(const std::string &taxon, const std::string &sequence,
                    int64_t position0, const std::string &occurrenceSuffix = "") {
    ResidueId id;
    id.taxon = taxon;
    id.sequence = sequence;
    id.occurrencePath = taxon + "|" + sequence + occurrenceSuffix;
    id.forwardPosition0 = position0;
    return id;
}

void addObservation(AlignmentEvidenceSet &evidence, const std::string &taxon,
                    const std::string &sequence, int64_t sourceSize,
                    int64_t position0, char base,
                    const std::string &occurrenceSuffix = "") {
    ResidueObservation observation;
    observation.id = residueId(taxon, sequence, position0, occurrenceSuffix);
    observation.sourceSize = sourceSize;
    observation.base = base;
    evidence.residues[observation.id] = observation;
}

std::vector<std::string> codes(const FastaEvidenceValidationReport &report) {
    std::vector<std::string> result;
    for (const FastaEvidenceDiagnostic &diagnostic : report.diagnostics) {
        result.push_back(diagnostic.code);
    }
    return result;
}

TEST(FastaEvidenceValidator, StreamsPlainFastaAndAcceptsIupacIntersections) {
    TempFastaDirectory files;
    const std::string fasta = files.plain(
        "T.fa", ">chr1 description\nAcRyN\n>chr2\ntt\n>unused\nG\n");

    AlignmentEvidenceSet evidence;
    evidence.sourceSizes[std::make_pair("T", "chr1")] = 5;
    // chr2 intentionally has no sourceSizes entry: its requirement is derived
    // from the residue observation itself.
    addObservation(evidence, "T", "chr1", 5, 0, 'A');
    addObservation(evidence, "T", "chr1", 5, 2, 'G');  // G intersects R
    addObservation(evidence, "T", "chr1", 5, 3, 'C');  // C intersects Y
    addObservation(evidence, "T", "chr1", 5, 4, 'A');  // N accepts A
    addObservation(evidence, "T", "chr2", 2, 1, 'T');

    const FastaEvidenceValidationCounts counts = FastaEvidenceValidator::validate(
        evidence, FastaEvidenceValidator::TaxonFastaMap{{"T", fasta}});
    EXPECT_EQ(1u, counts.taxaRequired);
    EXPECT_EQ(1u, counts.fastaFilesScanned);
    EXPECT_EQ(3u, counts.sequencesScanned);
    EXPECT_EQ(8u, counts.basesScanned);
    EXPECT_EQ(2u, counts.sourcesRequired);
    EXPECT_EQ(2u, counts.sourcesFound);
    EXPECT_EQ(2u, counts.sourcesValidated);
    EXPECT_EQ(5u, counts.residueObservationsRequired);
    EXPECT_EQ(5u, counts.residueObservationsValidated);
}

TEST(FastaEvidenceValidator, ReadsGzipAndValidatesNegativeStrandNormalizedBases) {
    TempFastaDirectory files;
    const std::string aFasta = files.plain("I1.fa", ">a\nAAANNNNN\n");
    const std::string bFasta = files.gzip("I2.fa.gz", ">b\nNNNNNGCA\n");

    std::istringstream maf(
        "a score=3\n"
        "s a 0 3 + 8 AAA\n"
        "s b 0 3 - 8 TGC\n");
    const std::vector<MafBlock> blocks =
        MafEvidenceReader::read(maf, "negative.maf", "I1", "I2");
    PairwiseMafInput input;
    input.leftTaxon = "I1";
    input.rightTaxon = "I2";
    input.mafPath = "negative.maf";
    input.runId = "negative";
    AlignmentEvidenceSet evidence;
    PairwiseEvidenceNormalizer::appendBlocks(input, blocks, evidence);

    const FastaEvidenceValidationCounts counts = FastaEvidenceValidator::validate(
        evidence, {{"I1", aFasta}, {"I2", bFasta}});
    EXPECT_EQ(2u, counts.fastaFilesScanned);
    EXPECT_EQ(2u, counts.sourcesValidated);
    EXPECT_EQ(6u, counts.residueObservationsValidated);
}

TEST(FastaEvidenceValidator, ReportsMissingDuplicateSizeAndBaseErrors) {
    TempFastaDirectory files;
    const std::string fasta = files.plain(
        "T.fa", ">chr1\nACGT\n>chr1 duplicate\nACGT\n>short\nAAA\n");

    AlignmentEvidenceSet evidence;
    evidence.sourceSizes[std::make_pair("T", "chr1")] = 4;
    evidence.sourceSizes[std::make_pair("T", "short")] = 4;
    evidence.sourceSizes[std::make_pair("T", "missing")] = 2;
    addObservation(evidence, "T", "chr1", 4, 0, 'C');

    const FastaEvidenceValidationReport report = FastaEvidenceValidator::inspect(
        evidence, {{"T", fasta}});
    EXPECT_FALSE(report.valid());
    EXPECT_EQ((std::vector<std::string>{
                  "FASTA_DUPLICATE_SEQUENCE_ID",
                  "FASTA_EVIDENCE_BASE_MISMATCH",
                  "FASTA_MISSING_SEQUENCE_ID",
                  "FASTA_SOURCE_SIZE_MISMATCH"}),
              codes(report));
    EXPECT_TRUE(std::is_sorted(report.diagnostics.begin(),
                               report.diagnostics.end()));
    EXPECT_EQ(3u, report.counts.sourcesRequired);
    EXPECT_EQ(2u, report.counts.sourcesFound);
    EXPECT_EQ(0u, report.counts.sourcesValidated);

    // Inspection is deterministic for the same bytes and evidence.
    const FastaEvidenceValidationReport repeated = FastaEvidenceValidator::inspect(
        evidence, {{"T", fasta}});
    ASSERT_EQ(report.diagnostics.size(), repeated.diagnostics.size());
    for (std::size_t i = 0; i < report.diagnostics.size(); ++i) {
        EXPECT_EQ(report.diagnostics[i].canonicalString(),
                  repeated.diagnostics[i].canonicalString());
    }
}

TEST(FastaEvidenceValidator, DiagnosesMalformedFastaWithoutUnboundedPerBaseErrors) {
    TempFastaDirectory files;
    const std::string fasta = files.plain(
        "bad.fa", "AC\n>\nA!Z%\n>chr1\nA\n");
    AlignmentEvidenceSet evidence;
    evidence.sourceSizes[std::make_pair("T", "chr1")] = 1;

    const FastaEvidenceValidationReport report = FastaEvidenceValidator::inspect(
        evidence, {{"T", fasta}});
    EXPECT_EQ((std::vector<std::string>{
                  "FASTA_EMPTY_SEQUENCE_ID", "FASTA_INVALID_BASE",
                  "FASTA_SEQUENCE_BEFORE_HEADER"}),
              codes(report));
    ASSERT_EQ(3u, report.diagnostics.size());
    EXPECT_NE(std::string::npos,
              report.diagnostics[1].message.find("3 non-IUPAC symbol(s)"));
    EXPECT_EQ(1u, report.counts.sourcesValidated);
}

TEST(FastaEvidenceValidator, MissingTaxonFastaIsAValidationError) {
    AlignmentEvidenceSet evidence;
    evidence.sourceSizes[std::make_pair("T", "chr1")] = 1;
    const FastaEvidenceValidationReport report =
        FastaEvidenceValidator::inspect(evidence, {});
    EXPECT_EQ((std::vector<std::string>{"FASTA_MISSING_SEQUENCE_ID",
                                        "MISSING_TAXON_FASTA"}),
              codes(report));
    EXPECT_EQ(0u, report.counts.fastaFilesScanned);
}

TEST(FastaEvidenceValidator, ValidateThrowsAndPreservesReportAndCounts) {
    TempFastaDirectory files;
    const std::string fasta = files.gzip("T.fa.gz", ">chr1\nA\n");
    AlignmentEvidenceSet evidence;
    evidence.sourceSizes[std::make_pair("T", "chr1")] = 1;
    addObservation(evidence, "T", "chr1", 1, 0, 'C');

    try {
        (void)FastaEvidenceValidator::validate(evidence, {{"T", fasta}});
        FAIL() << "expected FastaEvidenceValidationError";
    } catch (const FastaEvidenceValidationError &error) {
        ASSERT_EQ(1u, error.report().diagnostics.size());
        EXPECT_EQ("FASTA_EVIDENCE_BASE_MISMATCH",
                  error.report().diagnostics[0].code);
        EXPECT_EQ(1u, error.report().counts.fastaFilesScanned);
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("position0=0"));
    }
}

TEST(FastaEvidenceValidator, RejectsInvalidEvidenceCoordinatesAndSymbols) {
    TempFastaDirectory files;
    const std::string fasta = files.plain("T.fa", ">chr1\nA\n");
    AlignmentEvidenceSet evidence;
    evidence.sourceSizes[std::make_pair("T", "chr1")] = 1;
    addObservation(evidence, "T", "chr1", 1, 2, 'A');
    addObservation(evidence, "T", "chr1", 1, 0, '!');

    const FastaEvidenceValidationReport report = FastaEvidenceValidator::inspect(
        evidence, {{"T", fasta}});
    EXPECT_EQ((std::vector<std::string>{"EVIDENCE_INVALID_BASE",
                                        "EVIDENCE_POSITION_OUT_OF_RANGE"}),
              codes(report));
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
