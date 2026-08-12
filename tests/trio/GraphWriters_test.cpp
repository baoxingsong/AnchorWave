#include "src/trio/io/GraphWriters.h"

#include "gtest/gtest.h"

#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

void addWriterPair(const std::string &leftTaxon, const std::string &rightTaxon,
                   const std::string &leftText, const std::string &rightText,
                   AlignmentEvidenceSet &evidence) {
    std::size_t leftSize = 0, rightSize = 0;
    for (char c : leftText) leftSize += c != '-';
    for (char c : rightText) rightSize += c != '-';
    std::ostringstream maf;
    maf << "a score=2\n"
        << "s " << leftTaxon << "chr 0 " << leftSize << " + " << leftSize
        << ' ' << leftText << "\n"
        << "s " << rightTaxon << "chr 0 " << rightSize << " + " << rightSize
        << ' ' << rightText << "\n";
    std::istringstream stream(maf.str());
    PairwiseMafInput input;
    input.leftTaxon = leftTaxon;
    input.rightTaxon = rightTaxon;
    input.mafPath = leftTaxon + rightTaxon + ".maf";
    input.runId = leftTaxon + rightTaxon;
    PairwiseEvidenceNormalizer::appendBlocks(
        input,
        MafEvidenceReader::read(stream, input.mafPath, leftTaxon, rightTaxon),
        evidence);
}

AlignmentSiteGraph graphFixture(bool missingOutgroupBase) {
    AlignmentEvidenceSet evidence;
    if (missingOutgroupBase) {
        addWriterPair("I1", "I2", "ACG", "ATG", evidence);
        addWriterPair("I1", "O1", "ACG", "A-G", evidence);
        addWriterPair("I2", "O1", "ATG", "A-G", evidence);
    } else {
        addWriterPair("I1", "I2", "AC", "AT", evidence);
        addWriterPair("I1", "O1", "AC", "AC", evidence);
        addWriterPair("I2", "O1", "AT", "AC", evidence);
    }
    AlignmentGraphBuildOptions options;
    options.ingroup1 = "I1";
    options.ingroup2 = "I2";
    options.primaryOutgroup = "O1";
    return AlignmentSiteGraphBuilder::build(evidence, options);
}

TEST(GraphWriters, GfaContainsStableSitesAllelesAndExactPaths) {
    const AlignmentSiteGraph graph = graphFixture(false);
    std::ostringstream output;
    GfaWriter::write(graph, output);
    const std::string gfa = output.str();
    EXPECT_NE(std::string::npos, gfa.find("H\tVN:Z:1.1"));
    EXPECT_NE(std::string::npos, gfa.find("\tSI:Z:site_"));
    EXPECT_NE(std::string::npos, gfa.find("P\tpath_"));
    EXPECT_NE(std::string::npos, gfa.find("\tCH:Z:" + graph.coreHash()));
}

TEST(GraphWriters, MetadataDistinguishesObservedBaseAndAlignedAbsence) {
    const AlignmentSiteGraph graph = graphFixture(true);
    std::ostringstream output;
    GfaWriter::writeMetadata(graph, output);
    EXPECT_NE(std::string::npos, output.str().find("OBSERVED_BASE"));
    EXPECT_NE(std::string::npos, output.str().find("ALIGNED_ABSENCE"));
}

TEST(GraphWriters, MafIsDerivedAndRoundTripParsableForLinearTriangle) {
    const AlignmentSiteGraph graph = graphFixture(false);
    std::ostringstream output;
    const MafExportResult result = MafGraphExporter::write(graph, output);
    EXPECT_EQ(1u, result.blocksWritten);
    std::istringstream parse(output.str());
    MafReadOptions options;
    options.requireExactlyTwoRows = false;
    const std::vector<MafBlock> blocks =
        MafEvidenceReader::read(parse, "projection.maf", "unused-left", "unused-right", options);
    ASSERT_EQ(1u, blocks.size());
    EXPECT_EQ(3u, blocks[0].rows.size());
    EXPECT_EQ(2u, blocks[0].rows[0].text.size());
}

TEST(GraphWriters, MissingCoverageChangesPathSetAndSplitsProjection) {
    const AlignmentSiteGraph graph = graphFixture(true);
    std::ostringstream output;
    const MafExportResult result = MafGraphExporter::write(graph, output);
    // The first column has all three paths; the second is an explicit gap in O1,
    // so it remains exportable and is never confused with missing evidence.
    EXPECT_GE(result.blocksWritten, 1u);
    EXPECT_NE(std::string::npos, output.str().find("-"));
}

TEST(GraphWriters, WritesEveryLossyMafProjectionOmission) {
    MafExportResult result;
    MafExportOmission omission;
    omission.componentId = "component_1";
    omission.reason = "non_unique_or_cyclic_site_order";
    omission.siteIds = {"site_a", "site_b"};
    result.omissions.push_back(omission);
    std::ostringstream output;
    MafGraphExporter::writeOmissions(result, "core_1", output);
    EXPECT_NE(std::string::npos,
              output.str().find("component_1\tnon_unique_or_cyclic_site_order\t2\t"
                                "site_a,site_b\tcore_1"));
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
