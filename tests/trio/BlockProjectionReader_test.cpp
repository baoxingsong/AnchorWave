#include "src/trio/io/BlockProjectionReader.h"

#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <unistd.h>

namespace anchorwave {
namespace trio {
namespace {

class ProjectionTempFile {
public:
    explicit ProjectionTempFile(const std::string &contents) {
        char pattern[] = "/tmp/anchorwave-block-projections-XXXXXX";
        const int descriptor = mkstemp(pattern);
        if (descriptor < 0) throw std::runtime_error("mkstemp failed");
        close(descriptor);
        path_ = pattern;
        std::ofstream output(path_.c_str(), std::ios::binary);
        output << contents;
        if (!output) throw std::runtime_error("temporary projection write failed");
    }
    ~ProjectionTempFile() { std::remove(path_.c_str()); }
    const std::string &path() const { return path_; }
private:
    std::string path_;
};

AlignmentSiteGraph projectionGraph() {
    AlignmentEvidenceSet evidence;
    for (const auto &pair : std::vector<std::pair<std::string, std::string>>{
             {"I1", "I2"}, {"I1", "O1"}, {"I2", "O1"}}) {
        std::istringstream stream(
            "a score=1\n"
            "s chr1 0 4 + 4 ACGT\n"
            "s chr1 0 4 + 4 ACGT\n");
        PairwiseMafInput input;
        input.leftTaxon = pair.first;
        input.rightTaxon = pair.second;
        input.runId = pair.first + pair.second;
        input.mafPath = input.runId + ".maf";
        PairwiseEvidenceNormalizer::appendBlocks(
            input, MafEvidenceReader::read(stream, input.mafPath,
                                           pair.first, pair.second), evidence);
    }
    AlignmentGraphBuildOptions options;
    options.ingroup1 = "I1";
    options.ingroup2 = "I2";
    options.primaryOutgroup = "O1";
    return AlignmentSiteGraphBuilder::build(evidence, options);
}

TaxonManifest projectionTaxa() {
    TaxonManifest result;
    for (const auto &entry : std::vector<std::pair<std::string, TaxonRole>>{
             {"I1", TaxonRole::INGROUP_REFERENCE},
             {"I2", TaxonRole::INGROUP},
             {"O1", TaxonRole::PRIMARY_OUTGROUP}}) {
        TaxonRecord record;
        record.taxonId = entry.first;
        record.role = entry.second;
        result.records.push_back(record);
    }
    return result;
}

std::string manifestHeader() {
    return "#anchorwave-block-projections-v1\n"
           "projection_id\ttaxon_id\tsequence\tblock_id\tstart0\tend0\t"
           "orientation\tsource\n";
}

TEST(BlockProjectionReader, ParsesAndDeterministicallyIdsForwardCoordinates) {
    ProjectionTempFile file(
        manifestHeader() +
        ".\tI1\tchr1\tB1\t0\t2\t+\tmacro-synteny\n"
        ".\tI1\tchr1\tB2\t2\t4\t-\tmacro-synteny\n");
    const AlignmentSiteGraph graph = projectionGraph();
    const BlockProjectionSet result =
        BlockProjectionReader::read(file.path(), projectionTaxa(), graph);
    ASSERT_EQ(2u, result.projections.size());
    EXPECT_EQ(graph.coreHash(), result.immutableCoreHash);
    EXPECT_EQ("I1|chr1", result.projections[0].occurrencePath);
    EXPECT_FALSE(result.projections[0].projectionId.empty());
    EXPECT_EQ(BlockOrientation::REVERSE, result.projections[1].orientation);
}

TEST(BlockProjectionReader, RejectsUnknownTaxonAndInvalidCoordinates) {
    const AlignmentSiteGraph graph = projectionGraph();
    ProjectionTempFile unknown(
        manifestHeader() +
        ".\tBAD\tchr1\tB1\t0\t2\t+\tx\n");
    EXPECT_THROW(BlockProjectionReader::read(unknown.path(), projectionTaxa(), graph),
                 ManifestError);
    ProjectionTempFile invalid(
        manifestHeader() +
        ".\tI1\tchr1\tB1\t2\t2\t+\tx\n");
    EXPECT_THROW(BlockProjectionReader::read(invalid.path(), projectionTaxa(), graph),
                 ManifestError);
}

TEST(BlockProjectionReader, RejectsASequencePathAbsentFromGraph) {
    ProjectionTempFile file(
        manifestHeader() +
        ".\tI1\tmissing_chr\tB1\t0\t2\t+\tx\n");
    EXPECT_THROW(BlockProjectionReader::read(file.path(), projectionTaxa(),
                                             projectionGraph()),
                 ManifestError);
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
