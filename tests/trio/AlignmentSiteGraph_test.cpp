#include "src/trio/model/AlignmentSiteGraph.h"

#include "gtest/gtest.h"

#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

void addPair(const std::string &leftTaxon, const std::string &rightTaxon,
             const std::string &leftText, const std::string &rightText,
             const std::string &run, AlignmentEvidenceSet &evidence) {
    std::size_t leftSize = 0;
    std::size_t rightSize = 0;
    for (char c : leftText) leftSize += c != '-';
    for (char c : rightText) rightSize += c != '-';
    std::ostringstream maf;
    maf << "a score=10\n"
        << "s " << leftTaxon << "chr 0 " << leftSize << " + 100"
        << ' ' << leftText << "\n"
        << "s " << rightTaxon << "chr 0 " << rightSize << " + 100"
        << ' ' << rightText << "\n";
    std::istringstream stream(maf.str());
    PairwiseMafInput input;
    input.leftTaxon = leftTaxon;
    input.rightTaxon = rightTaxon;
    input.mafPath = run + ".maf";
    input.runId = run;
    PairwiseEvidenceNormalizer::appendBlocks(
        input, MafEvidenceReader::read(stream, input.mafPath, leftTaxon, rightTaxon),
        evidence);
}

AlignmentGraphBuildOptions options() {
    AlignmentGraphBuildOptions result;
    result.ingroup1 = "I1";
    result.ingroup2 = "I2";
    result.primaryOutgroup = "O1";
    return result;
}

TEST(AlignmentSiteGraph, ImportsAConcordantTriangleDeterministically) {
    AlignmentEvidenceSet first;
    addPair("I1", "I2", "ACG", "ATG", "12", first);
    addPair("I1", "O1", "ACG", "ACG", "1o", first);
    addPair("I2", "O1", "ATG", "ACG", "2o", first);
    const AlignmentSiteGraph graph = AlignmentSiteGraphBuilder::build(first, options());
    ASSERT_EQ(3u, graph.sites().size());
    for (const auto &site : graph.sites()) {
        EXPECT_EQ(ConsistencyClass::CONSISTENT, site.second.consistency);
        EXPECT_EQ(3u, site.second.directPairSupport.size());
    }
    EXPECT_FALSE(graph.coreHash().empty());

    AlignmentEvidenceSet reordered;
    addPair("I2", "O1", "ATG", "ACG", "2o", reordered);
    addPair("I1", "O1", "ACG", "ACG", "1o", reordered);
    addPair("I1", "I2", "ACG", "ATG", "12", reordered);
    EXPECT_EQ(graph.coreHash(), AlignmentSiteGraphBuilder::build(reordered, options()).coreHash());
}

TEST(AlignmentSiteGraph, RequiresAllPrimaryPairRuns) {
    AlignmentEvidenceSet evidence;
    addPair("I1", "I2", "A", "A", "12", evidence);
    addPair("I1", "O1", "A", "A", "1o", evidence);
    EXPECT_THROW(AlignmentSiteGraphBuilder::build(evidence, options()), AlignmentGraphError);
}

TEST(AlignmentSiteGraph, RejectsHardCopyContradictionBeforeGraphUnion) {
    AlignmentEvidenceSet evidence;
    addPair("I1", "I2", "A", "A", "12", evidence);
    addPair("I1", "O1", "A", "A", "1o", evidence);
    addPair("I2", "O1", "A", "A", "2o", evidence);
    AlignmentGraphBuildOptions build = options();
    for (const auto &residue : evidence.residues) {
        CopyAssignment assignment;
        assignment.copyGroup = residue.first.taxon == "O1" ? "ancient_B" : "ancient_A";
        assignment.hard = true;
        build.copyAssignments[residue.first] = assignment;
    }
    const AlignmentSiteGraph graph = AlignmentSiteGraphBuilder::build(evidence, build);
    bool found = false;
    for (const GraphConflict &conflict : graph.conflicts()) {
        found |= conflict.classification == ConsistencyClass::COPY_COLLISION &&
                 !conflict.poaEligible;
    }
    EXPECT_TRUE(found);
}

TEST(AlignmentSiteGraph, MissingEvidenceIsNotEncodedAsDeletion) {
    AlignmentEvidenceSet evidence;
    addPair("I1", "I2", "AC", "AC", "12", evidence);
    addPair("I1", "O1", "AC", "A-", "1o", evidence);
    addPair("I2", "O1", "AC", "A-", "2o", evidence);
    const AlignmentSiteGraph graph = AlignmentSiteGraphBuilder::build(evidence, options());
    std::size_t incomplete = 0;
    for (const auto &site : graph.sites()) {
        if (site.second.consistency == ConsistencyClass::MISSING_EVIDENCE) {
            ++incomplete;
        }
    }
    EXPECT_EQ(1u, incomplete);
    // The O1 path has one observed base; no synthetic gap/base was inserted.
    std::size_t outgroupBases = 0;
    for (const PathSegment &path : graph.pathSegments()) {
        if (path.taxon == "O1") outgroupBases += path.sequenceText.size();
    }
    EXPECT_EQ(1u, outgroupBases);
    for (const auto &site : graph.sites()) {
        EXPECT_TRUE(site.second.alignedAbsentPaths.empty());
    }
    bool terminalGapAuditedAsMissing = false;
    for (const EvidenceDisposition &audit : graph.evidenceAudit()) {
        terminalGapAuditedAsMissing |=
            !audit.accepted &&
            audit.classification == ConsistencyClass::MISSING_EVIDENCE &&
            audit.reason == "aligned gap lacks two observed flanks";
    }
    EXPECT_TRUE(terminalGapAuditedAsMissing);
}

TEST(AlignmentSiteGraph, InternallyFlankedGapIsExplicitAbsence) {
    AlignmentEvidenceSet evidence;
    addPair("I1", "I2", "ACG", "ACG", "12", evidence);
    addPair("I1", "O1", "ACG", "A-G", "1o", evidence);
    addPair("I2", "O1", "ACG", "A-G", "2o", evidence);
    const AlignmentSiteGraph graph = AlignmentSiteGraphBuilder::build(evidence, options());
    bool explicitAbsence = false;
    for (const auto &site : graph.sites()) {
        explicitAbsence |= site.second.alignedAbsentPaths.count("O1|O1chr") != 0;
    }
    EXPECT_TRUE(explicitAbsence);
}

TEST(AlignmentSiteGraph, DetectsGapVersusTransitiveResidueConflict) {
    AlignmentEvidenceSet evidence;
    addPair("I1", "I2", "ACG", "ACG", "12", evidence);
    addPair("I1", "O1", "ACG", "A-C", "1o", evidence);
    addPair("I2", "O1", "ACG", "ACG", "2o", evidence);
    const AlignmentSiteGraph graph = AlignmentSiteGraphBuilder::build(evidence, options());
    bool poaCandidate = false;
    for (const GraphConflict &conflict : graph.conflicts()) {
        poaCandidate |= conflict.classification == ConsistencyClass::RESIDUE_TRANSITIVITY &&
                        conflict.poaEligible;
    }
    EXPECT_TRUE(poaCandidate);
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
