#include "../../src/trio/io/CopyConstraintReader.h"
#include "../../src/trio/io/PairwiseManifestReader.h"
#include "../../src/trio/io/TaxonManifestReader.h"

#include <gtest/gtest.h>

#include <fcntl.h>
#include <unistd.h>

#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using anchorwave::trio::CanonicalTaxonPair;
using anchorwave::trio::ConstraintStrength;
using anchorwave::trio::CopyConstraintReader;
using anchorwave::trio::CopyMemberType;
using anchorwave::trio::CopyRecordType;
using anchorwave::trio::ManifestError;
using anchorwave::trio::ManifestReaderOptions;
using anchorwave::trio::PairwiseManifestReader;
using anchorwave::trio::PairwiseManifestReaderOptions;
using anchorwave::trio::PairwiseScope;
using anchorwave::trio::TaxonManifest;
using anchorwave::trio::TaxonManifestReader;

class TempFile {
public:
    explicit TempFile(const std::string &contents) {
        char pattern[] = "/tmp/anchorwave-trio-manifest-XXXXXX";
        const int descriptor = mkstemp(pattern);
        if (descriptor < 0) {
            throw std::runtime_error("mkstemp failed");
        }
        close(descriptor);
        path_ = pattern;
        std::ofstream output(path_.c_str(), std::ios::binary);
        if (!output.good()) {
            throw std::runtime_error("cannot create temporary manifest");
        }
        output << contents;
        output.close();
    }

    ~TempFile() {
        unlink(path_.c_str());
    }

    const std::string &path() const {
        return path_;
    }

private:
    std::string path_;
};

class TempDirectory {
public:
    TempDirectory() {
        char pattern[] = "/tmp/anchorwave-trio-manifest-dir-XXXXXX";
        char *created = mkdtemp(pattern);
        if (created == nullptr) {
            throw std::runtime_error("mkdtemp failed");
        }
        path_ = created;
    }

    ~TempDirectory() {
        for (std::vector<std::string>::const_reverse_iterator file = files_.rbegin();
             file != files_.rend(); ++file) {
            unlink(file->c_str());
        }
        rmdir(path_.c_str());
    }

    std::string write(const std::string &name, const std::string &contents = "") {
        const std::string filePath = path_ + "/" + name;
        std::ofstream output(filePath.c_str(), std::ios::binary);
        if (!output.good()) {
            throw std::runtime_error("cannot create temporary fixture");
        }
        output << contents;
        output.close();
        files_.push_back(filePath);
        return filePath;
    }

    const std::string &path() const {
        return path_;
    }

private:
    std::string path_;
    std::vector<std::string> files_;
};

class ScopedWorkingDirectory {
public:
    explicit ScopedWorkingDirectory(const std::string &path) {
        char *current = getcwd(nullptr, 0);
        if (current == nullptr) {
            throw std::runtime_error("getcwd failed");
        }
        original_ = current;
        std::free(current);
        if (chdir(path.c_str()) != 0) {
            throw std::runtime_error("chdir failed");
        }
        current = getcwd(nullptr, 0);
        if (current == nullptr) {
            chdir(original_.c_str());
            throw std::runtime_error("getcwd after chdir failed");
        }
        current_ = current;
        std::free(current);
    }

    ~ScopedWorkingDirectory() {
        chdir(original_.c_str());
    }

    const std::string &path() const {
        return current_;
    }

private:
    std::string original_;
    std::string current_;
};

std::string taxonHeader() {
    return "taxon_id\trole\tfasta\tgff\tanchor_sam\tanchor_fasta\t"
           "callability_bed\tquality_weight\n";
}

std::string fourTaxa() {
    return std::string("# manifest comment\n") + taxonHeader() +
           "R\tingroup_reference\tR.fa\t.\t.\t.\t.\t1.0\n"
           "S\tingroup\tS.fa\t.\t.\t.\t.\t0.95\n"
           "O1\tprimary_outgroup\tO1.fa\t.\t.\t.\t.\t0.9\n"
           "O2\toutgroup\tO2.fa\t.\t.\t.\t.\t0.8\n";
}

std::string pairwiseHeader() {
    return "taxon_a\ttaxon_b\tmaf\tanchor_map\tscore_profile\tweight\n";
}

std::string requiredTrianglePairs() {
    return pairwiseHeader() +
           "S\tR\tR_S.maf\t.\tingroup\t1.0\n"
           "R\tO1\tR_O1.maf\t.\toutgroup\t0.9\n"
           "O1\tS\tS_O1.maf\t.\toutgroup\t0.9\n"
           "R\tO2\tR_O2.maf\t.\toutgroup\t0.8\n"
           "S\tO2\tS_O2.maf\t.\toutgroup\t0.8\n";
}

std::string copyHeader() {
    return "record_type\tfamily_id\tancestral_copy_id\tparent_copy_id\t"
           "born_on_edge\ttaxon_id\tmember_type\tmember_id\tseq\tstart0\t"
           "end0\tstrand\texpected_count\trelation\tconfidence\tconstraint\t"
           "source\n";
}

TaxonManifest readFourTaxa() {
    TempFile file(fourTaxa());
    return TaxonManifestReader().read(file.path());
}

TEST(TaxonManifestReaderTest, ParsesCommentsRolesAndOptionalDotPaths) {
    TempFile file(fourTaxa());
    const TaxonManifest manifest = TaxonManifestReader().read(file.path());

    ASSERT_EQ(4U, manifest.records.size());
    ASSERT_EQ(2U, manifest.ingroupIds().size());
    EXPECT_EQ("R", manifest.ingroupIds()[0]);
    EXPECT_EQ("S", manifest.ingroupIds()[1]);
    EXPECT_EQ("O1", manifest.primaryOutgroupId());
    ASSERT_EQ(2U, manifest.outgroupIds().size());
    EXPECT_EQ("O1", manifest.outgroupIds()[0]);
    EXPECT_EQ("O2", manifest.outgroupIds()[1]);
    EXPECT_DOUBLE_EQ(0.95, manifest.find("S")->qualityWeight);
}

TEST(TaxonManifestReaderTest, RejectsDuplicateTaxonWithLineNumber) {
    TempFile file(taxonHeader() +
                  "R\tingroup_reference\tR.fa\t.\t.\t.\t.\t1\n"
                  "S\tingroup\tS.fa\t.\t.\t.\t.\t1\n"
                  "R\tprimary_outgroup\tO.fa\t.\t.\t.\t.\t1\n");
    try {
        TaxonManifestReader().read(file.path());
        FAIL() << "expected ManifestError";
    } catch (const ManifestError &error) {
        EXPECT_EQ(4U, error.location().line);
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("duplicate taxon_id 'R'"));
    }
}

TEST(TaxonManifestReaderTest, RejectsTaxonIdsThatBreakOccurrencePathEncoding) {
    const std::vector<std::string> invalidIds = {"bad|taxon", "bad taxon",
                                                  std::string("bad\vcontrol")};
    for (std::vector<std::string>::const_iterator taxonId = invalidIds.begin();
         taxonId != invalidIds.end(); ++taxonId) {
        TempFile file(taxonHeader() +
                      *taxonId + "\tingroup_reference\tR.fa\t.\t.\t.\t.\t1\n"
                      "S\tingroup\tS.fa\t.\t.\t.\t.\t1\n"
                      "O\tprimary_outgroup\tO.fa\t.\t.\t.\t.\t1\n");
        try {
            TaxonManifestReader().read(file.path());
            FAIL() << "expected invalid taxon_id to be rejected: " << *taxonId;
        } catch (const ManifestError &error) {
            EXPECT_EQ(2U, error.location().line);
            EXPECT_NE(std::string::npos,
                      std::string(error.what()).find("taxon_id must not contain"));
        }
    }
}

TEST(TaxonManifestReaderTest, OptionalExistenceValidationIsStrict) {
    ManifestReaderOptions options;
    options.validatePaths = true;
    TempFile file(taxonHeader() +
                  "R\tingroup_reference\t/no/such/R.fa\t.\t.\t.\t.\t1\n"
                  "S\tingroup\t/no/such/S.fa\t.\t.\t.\t.\t1\n"
                  "O\tprimary_outgroup\t/no/such/O.fa\t.\t.\t.\t.\t1\n");
    EXPECT_THROW(TaxonManifestReader(options).read(file.path()), ManifestError);
}

TEST(ManifestReaderPathTest, ResolvesRelativeTaxonAndPairwisePathsAtManifestDirectory) {
    TempDirectory directory;
    const std::string absoluteReference = directory.write("absolute-R.fa", ">chr1\nA\n");
    directory.write("S.fa", ">chr1\nA\n");
    directory.write("O1.fa", ">chr1\nA\n");
    directory.write("O2.fa", ">chr1\nA\n");
    directory.write("R.gff", "##gff-version 3\n");
    directory.write("R.sam", "@HD\tVN:1.6\n");
    directory.write("R.anchors.fa", ">a\nA\n");
    directory.write("R.callable.bed", "chr1\t0\t1\n");

    directory.write(
        "taxa.tsv", taxonHeader() +
                        "R\tingroup_reference\t" + absoluteReference +
                        "\tR.gff\tR.sam\tR.anchors.fa\tR.callable.bed\t1\n"
                        "S\tingroup\tS.fa\t.\t.\t.\t.\t1\n"
                        "O1\tprimary_outgroup\tO1.fa\t.\t.\t.\t.\t1\n"
                        "O2\toutgroup\tO2.fa\t.\t.\t.\t.\t1\n");
    const std::size_t directorySeparator = directory.path().find_last_of('/');
    ASSERT_NE(std::string::npos, directorySeparator);
    const std::string parentDirectory = directory.path().substr(0, directorySeparator);
    const std::string relativeDirectory = directory.path().substr(directorySeparator + 1);
    ScopedWorkingDirectory workingDirectory(parentDirectory);
    const std::string resolvedDirectory =
        workingDirectory.path() + "/" + relativeDirectory;
    ManifestReaderOptions taxonOptions;
    taxonOptions.validatePaths = true;
    const TaxonManifest taxa = TaxonManifestReader(taxonOptions).read(
        relativeDirectory + "/taxa.tsv");

    ASSERT_NE(nullptr, taxa.find("R"));
    EXPECT_EQ(absoluteReference, taxa.find("R")->fasta);
    EXPECT_EQ(resolvedDirectory + "/R.gff", taxa.find("R")->gff);
    EXPECT_EQ(resolvedDirectory + "/R.sam", taxa.find("R")->anchorSam);
    EXPECT_EQ(resolvedDirectory + "/R.anchors.fa", taxa.find("R")->anchorFasta);
    EXPECT_EQ(resolvedDirectory + "/R.callable.bed", taxa.find("R")->callabilityBed);
    EXPECT_EQ(resolvedDirectory + "/S.fa", taxa.find("S")->fasta);
    EXPECT_EQ(".", taxa.find("S")->gff);

    directory.write("R_S.maf");
    directory.write("R_O1.maf");
    directory.write("S_O1.maf");
    directory.write("R_O2.maf");
    const std::string absoluteMaf = directory.write("absolute-S_O2.maf");
    directory.write("R_S.anchors.tsv");
    directory.write(
        "pairwise.tsv", pairwiseHeader() +
                            "S\tR\tR_S.maf\tR_S.anchors.tsv\tingroup\t1\n"
                            "R\tO1\tR_O1.maf\t.\toutgroup\t1\n"
                            "O1\tS\tS_O1.maf\t.\toutgroup\t1\n"
                            "R\tO2\tR_O2.maf\t.\toutgroup\t1\n"
                            "S\tO2\t" + absoluteMaf + "\t.\toutgroup\t1\n");
    PairwiseManifestReaderOptions pairOptions;
    pairOptions.validatePaths = true;
    const anchorwave::trio::PairwiseManifest pairs =
        PairwiseManifestReader(pairOptions).read(
            relativeDirectory + "/pairwise.tsv", taxa);

    const anchorwave::trio::PairwiseManifestRecord *rs =
        pairs.find(CanonicalTaxonPair("R", "S"));
    ASSERT_NE(nullptr, rs);
    EXPECT_EQ(resolvedDirectory + "/R_S.maf", rs->maf);
    EXPECT_EQ(resolvedDirectory + "/R_S.anchors.tsv", rs->anchorMap);
    const anchorwave::trio::PairwiseManifestRecord *so2 =
        pairs.find(CanonicalTaxonPair("S", "O2"));
    ASSERT_NE(nullptr, so2);
    EXPECT_EQ(absoluteMaf, so2->maf);
    EXPECT_EQ(".", so2->anchorMap);
}

TEST(ManifestReaderNumericTest, RejectsFloatingPointOverflow) {
    TempFile taxaFile(taxonHeader() +
                       "R\tingroup_reference\tR.fa\t.\t.\t.\t.\t1e309\n"
                       "S\tingroup\tS.fa\t.\t.\t.\t.\t1\n"
                       "O\tprimary_outgroup\tO.fa\t.\t.\t.\t.\t1\n");
    EXPECT_THROW(TaxonManifestReader().read(taxaFile.path()), ManifestError);

    const TaxonManifest taxa = readFourTaxa();
    TempFile pairFile(pairwiseHeader() +
                      "R\tS\tR_S.maf\t.\tingroup\t1e309\n");
    EXPECT_THROW(PairwiseManifestReader().read(pairFile.path(), taxa), ManifestError);
}

TEST(PairwiseManifestReaderTest, CanonicalizesPairsAndRequiresEveryOutgroupTriangle) {
    const TaxonManifest taxa = readFourTaxa();
    TempFile file(requiredTrianglePairs());
    const anchorwave::trio::PairwiseManifest manifest =
        PairwiseManifestReader().read(file.path(), taxa);

    ASSERT_EQ(5U, manifest.records.size());
    const CanonicalTaxonPair rs("R", "S");
    ASSERT_TRUE(manifest.contains(rs));
    EXPECT_EQ("R--S", rs.stableId());
    EXPECT_EQ("S", manifest.find(rs)->sourceTaxonA);
    EXPECT_EQ("R", manifest.find(rs)->sourceTaxonB);
}

TEST(PairwiseManifestReaderTest, ReportsMissingExtraOutgroupEdge) {
    const TaxonManifest taxa = readFourTaxa();
    const std::string missing = pairwiseHeader() +
                                "R\tS\tR_S.maf\t.\tingroup\t1\n"
                                "R\tO1\tR_O1.maf\t.\toutgroup\t1\n"
                                "S\tO1\tS_O1.maf\t.\toutgroup\t1\n"
                                "R\tO2\tR_O2.maf\t.\toutgroup\t1\n";
    TempFile file(missing);
    try {
        PairwiseManifestReader().read(file.path(), taxa);
        FAIL() << "expected ManifestError";
    } catch (const ManifestError &error) {
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("O2--S"));
    }
}

TEST(PairwiseManifestReaderTest, CompleteScopeRequiresOutgroupPair) {
    const TaxonManifest taxa = readFourTaxa();
    TempFile file(requiredTrianglePairs());
    PairwiseManifestReaderOptions options;
    options.scope = PairwiseScope::COMPLETE;
    try {
        PairwiseManifestReader(options).read(file.path(), taxa);
        FAIL() << "expected ManifestError";
    } catch (const ManifestError &error) {
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("O1--O2"));
    }
}

TEST(PairwiseManifestReaderTest, RejectsReversedDuplicatePair) {
    const TaxonManifest taxa = readFourTaxa();
    TempFile file(requiredTrianglePairs() +
                  "R\tS\tduplicate.maf\t.\tingroup\t1\n");
    try {
        PairwiseManifestReader().read(file.path(), taxa);
        FAIL() << "expected ManifestError";
    } catch (const ManifestError &error) {
        EXPECT_EQ(7U, error.location().line);
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("duplicate unordered taxon pair"));
    }
}

TEST(CopyConstraintReaderTest, ParsesHardSoftIntervalsCountsAndFallbackCapacities) {
    const TaxonManifest taxa = readFourTaxa();
    const std::string contents =
        std::string("#anchorwave-copy-relations-v1\n# provenance comment\n") +
        copyHeader() +
        "GROUP\tF1\tA1\t.\tpre_root\t.\t.\t.\t.\t.\t.\t.\t.\t.\t1.0\thard\tcurated\n"
        "MEMBER\tF1\tA1\t.\t.\tR\tanchor\trGene1\t.\t.\t.\t.\t.\tortholog\t1.0\thard\tcurated\n"
        "MEMBER\tF1\tA1\t.\t.\tS\tinterval\t.\tchr1\t10\t20\t-\t.\tcoortholog\t0.8\tsoft\tsynteny\n"
        "COUNT\tF1\t.\t.\t.\tR\t.\t.\t.\t.\t.\t.\t2\t.\t1.0\thard\tcytogenetics\n";
    TempFile file(contents);
    const std::vector<std::string> defaults = {"R=3", "F1:S=4"};
    const anchorwave::trio::CopyConstraintSet result =
        CopyConstraintReader().read(file.path(), taxa, defaults);

    ASSERT_EQ(4U, result.records.size());
    EXPECT_EQ(CopyRecordType::GROUP, result.records[0].recordType);
    EXPECT_EQ(CopyMemberType::ANCHOR, result.records[1].memberType);
    EXPECT_EQ(ConstraintStrength::HARD, result.records[1].strength);
    ASSERT_TRUE(result.records[2].hasInterval);
    EXPECT_EQ(10, result.records[2].interval.start0);
    EXPECT_EQ(20, result.records[2].interval.end0);
    EXPECT_EQ(10, result.records[2].interval.length());
    EXPECT_EQ('-', result.records[2].interval.strand);
    ASSERT_TRUE(result.records[3].hasExpectedCount);
    EXPECT_EQ(2, result.records[3].expectedCount);

    std::int64_t count = -1;
    ASSERT_TRUE(result.findDefaultCapacity("F1", "S", count));
    EXPECT_EQ(4, count);
    ASSERT_TRUE(result.findDefaultCapacity("another-family", "R", count));
    EXPECT_EQ(3, count);
}

TEST(CopyConstraintReaderTest, RejectsNegativeHalfOpenIntervalWithLineNumber) {
    const TaxonManifest taxa = readFourTaxa();
    const std::string contents =
        std::string("#anchorwave-copy-relations-v1\n") + copyHeader() +
        "MEMBER\tF1\tA1\t.\t.\tR\tinterval\t.\tchr1\t-1\t20\t+\t.\tortholog\t1.0\thard\tcurated\n";
    TempFile file(contents);
    try {
        CopyConstraintReader().read(file.path(), taxa);
        FAIL() << "expected ManifestError";
    } catch (const ManifestError &error) {
        EXPECT_EQ(3U, error.location().line);
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("start0 must be non-negative"));
    }
}

TEST(CopyConstraintReaderTest, RejectsIncompatibleHardAssignments) {
    const TaxonManifest taxa = readFourTaxa();
    const std::string contents =
        std::string("#anchorwave-copy-relations-v1\n") + copyHeader() +
        "MEMBER\tF1\tA1\t.\t.\tR\tanchor\trGene1\t.\t.\t.\t.\t.\tortholog\t1.0\thard\tone\n"
        "MEMBER\tF1\tA2\t.\t.\tR\tanchor\trGene1\t.\t.\t.\t.\t.\tortholog\t1.0\thard\ttwo\n";
    TempFile file(contents);
    EXPECT_THROW(CopyConstraintReader().read(file.path(), taxa), ManifestError);
}

TEST(CopyConstraintReaderTest, RequiresVersionBeforeCommentsAndHeader) {
    const TaxonManifest taxa = readFourTaxa();
    TempFile file(std::string("# ordinary comment\n") + copyHeader());
    try {
        CopyConstraintReader().read(file.path(), taxa);
        FAIL() << "expected ManifestError";
    } catch (const ManifestError &error) {
        EXPECT_EQ(1U, error.location().line);
        EXPECT_NE(std::string::npos,
                  std::string(error.what()).find("schema version"));
    }
}

TEST(CopyConstraintReaderTest, ParsesAndValidatesDefaultCapacitySyntax) {
    const anchorwave::trio::CopyCapacity family =
        CopyConstraintReader::parseDefaultCapacitySpec("F1:R=2");
    EXPECT_EQ("F1", family.key.familyId);
    EXPECT_EQ("R", family.key.taxonId);
    EXPECT_EQ(2, family.count);

    EXPECT_THROW(CopyConstraintReader::parseDefaultCapacitySpec("F1:R=-1"),
                 ManifestError);
    EXPECT_THROW(CopyConstraintReader::parseDefaultCapacitySpec("F1:R"),
                 ManifestError);
}

}  // namespace
