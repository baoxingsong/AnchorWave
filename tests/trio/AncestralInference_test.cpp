#include "src/trio/evolution/AncestralInference.h"

#include "src/trio/impl/CopyRelationshipResolver.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

void addAncestorPair(const std::string &a, const std::string &b,
                     const std::string &textA, const std::string &textB,
                     AlignmentEvidenceSet &evidence) {
    const std::size_t sizeA = static_cast<std::size_t>(
        std::count_if(textA.begin(), textA.end(), [](char base) { return base != '-'; }));
    const std::size_t sizeB = static_cast<std::size_t>(
        std::count_if(textB.begin(), textB.end(), [](char base) { return base != '-'; }));
    std::ostringstream maf;
    maf << "a score=5\n"
        << "s " << a << "chr 0 " << sizeA << " + " << sizeA
        << ' ' << textA << "\n"
        << "s " << b << "chr 0 " << sizeB << " + " << sizeB
        << ' ' << textB << "\n";
    std::istringstream stream(maf.str());
    PairwiseMafInput input;
    input.leftTaxon = a; input.rightTaxon = b;
    input.runId = a + b; input.mafPath = input.runId;
    PairwiseEvidenceNormalizer::appendBlocks(
        input, MafEvidenceReader::read(stream, input.runId, a, b), evidence);
}

AlignmentSiteGraph ambiguousPresenceGraph() {
    AlignmentEvidenceSet evidence;
    addAncestorPair("I1", "I2", "ACA", "A-A", evidence);
    // The other required primary-triangle runs exist but have no callable
    // blocks for this local component, so O1 remains missing rather than
    // becoming an asserted deletion.
    evidence.observedPairs.insert(
        PairwiseEvidenceNormalizer::canonicalPairId("I1", "O1"));
    evidence.observedPairs.insert(
        PairwiseEvidenceNormalizer::canonicalPairId("I2", "O1"));
    const CopyResolutionResult copies = CopyRelationshipResolver::resolve(
        evidence, CopyConstraintSet(), CopyResolverOptions());
    AlignmentGraphBuildOptions options;
    options.ingroup1 = "I1";
    options.ingroup2 = "I2";
    options.primaryOutgroup = "O1";
    options.copyAssignments = copies.assignments;
    options.allowedEvidenceIds = copies.selectedEvidenceIds;
    options.restrictToAllowedEvidence = true;
    options.requireResolvedCopyAssignments = true;
    return AlignmentSiteGraphBuilder::build(evidence, options);
}

AlignmentSiteGraph ancestorGraph() {
    AlignmentEvidenceSet evidence;
    addAncestorPair("I1", "I2", "AG", "GG", evidence);
    addAncestorPair("I1", "O1", "AG", "GG", evidence);
    addAncestorPair("I2", "O1", "GG", "GG", evidence);
    addAncestorPair("I1", "O2", "AG", "GG", evidence);
    addAncestorPair("I2", "O2", "GG", "GG", evidence);
    const CopyResolutionResult copies = CopyRelationshipResolver::resolve(
        evidence, CopyConstraintSet(), CopyResolverOptions());
    AlignmentGraphBuildOptions options;
    options.ingroup1 = "I1"; options.ingroup2 = "I2"; options.primaryOutgroup = "O1";
    options.copyAssignments = copies.assignments;
    options.allowedEvidenceIds = copies.selectedEvidenceIds;
    options.restrictToAllowedEvidence = true;
    options.requireResolvedCopyAssignments = true;
    return AlignmentSiteGraphBuilder::build(evidence, options);
}

TEST(AncestralInference, MultipleOutgroupsPolarizeCallsOnFrozenGraph) {
    const AlignmentSiteGraph graph = ancestorGraph();
    const std::string before = graph.coreHash();
    const Phylogeny tree = Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    EvolutionInferenceConfig config;
    config.targetNode = "A";
    const AncestralOverlay overlay = AncestralInference::infer(graph, tree, config);
    EXPECT_EQ(before, graph.coreHash());
    EXPECT_EQ(before, overlay.immutableCoreHash);
    ASSERT_EQ(2u, overlay.calls.size());
    for (const auto &call : overlay.calls) {
        EXPECT_TRUE(call.second.emitted);
        EXPECT_EQ('G', call.second.emittedBase);
    }
}

TEST(AncestralInference, LinearSequenceDoesNotSplitAtPerSiteCopyLabels) {
    const AlignmentSiteGraph graph = ancestorGraph();
    const Phylogeny tree = Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    EvolutionInferenceConfig config;
    config.targetNode = "A";
    const AncestralOverlay overlay = AncestralInference::infer(graph, tree, config);
    ASSERT_EQ(1u, overlay.blocks.size());
    EXPECT_EQ(2u, overlay.blocks.front().siteIds.size());
    EXPECT_EQ("GG", overlay.blocks.front().sequence);
    // The sequence path is continuous, while its two independently inferred
    // per-site copy labels are intentionally not collapsed into one label.
    EXPECT_TRUE(overlay.blocks.front().copyGroup.empty());
}

TEST(AncestralInference, AmbiguousPresenceRequiresExplicitEmissionOptIn) {
    const AlignmentSiteGraph graph = ambiguousPresenceGraph();
    const Phylogeny tree = Phylogeny::fromNewick("((I1,I2)A,O1)Root;");
    EvolutionInferenceConfig conservative;
    conservative.targetNode = "A";
    const AncestralOverlay excluded =
        AncestralInference::infer(graph, tree, conservative);

    EvolutionInferenceConfig permissive = conservative;
    permissive.includeAmbiguousPresence = true;
    const AncestralOverlay included =
        AncestralInference::infer(graph, tree, permissive);

    bool sawAmbiguous = false;
    for (const auto &entry : excluded.calls) {
        if (!entry.second.presence.isAmbiguous()) continue;
        sawAmbiguous = true;
        EXPECT_FALSE(entry.second.emitted);
        EXPECT_EQ("PRESENCE_AMBIGUOUS", entry.second.decision);
        const AncestralSiteCall &optedIn = included.calls.at(entry.first);
        EXPECT_TRUE(optedIn.emitted);
        EXPECT_EQ("PRESENCE_AMBIGUOUS_EMITTED", optedIn.decision);
    }
    EXPECT_TRUE(sawAmbiguous);
}

TEST(AncestralInference, EvolutionModelChangesCallsButNeverAlignmentGraph) {
    const AlignmentSiteGraph graph = ancestorGraph();
    const std::string immutable = graph.coreHash();
    const Phylogeny tree = Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    EvolutionInferenceConfig first;
    first.targetNode = "A";
    EvolutionInferenceConfig second = first;
    second.nucleotideModel.transitionCost[2][0] = 50.0;
    second.nucleotideModel.transitionCost[0][2] = 50.0;
    (void)AncestralInference::infer(graph, tree, first);
    (void)AncestralInference::infer(graph, tree, second);
    EXPECT_EQ(immutable, graph.coreHash());
}

TEST(AncestralInference, WritesUnorderedBlocksAndAlignmentsToBothChildren) {
    const AlignmentSiteGraph graph = ancestorGraph();
    const Phylogeny tree = Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    EvolutionInferenceConfig config;
    config.targetNode = "A";
    const AncestralOverlay overlay = AncestralInference::infer(graph, tree, config);
    const std::vector<AncestorChildAlignmentBlock> projections =
        AncestralInference::projectToIngroupChildren(graph, overlay);
    ASSERT_FALSE(projections.empty());
    std::ostringstream fasta, maf, map;
    AncestralWriters::writeFasta(overlay, fasta);
    AncestralWriters::writeChildMaf(overlay, projections, maf);
    AncestralWriters::writeChildMap(projections, map);
    EXPECT_NE(std::string::npos, fasta.str().find(">ancestor_block_"));
    EXPECT_NE(std::string::npos, maf.str().find("s\tA.ancestor_block_"));
    EXPECT_NE(std::string::npos, maf.str().find("s\tI1.I1chr"));
    EXPECT_NE(std::string::npos, maf.str().find("s\tI2.I2chr"));
    EXPECT_NE(std::string::npos, map.str().find("child_path"));
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
