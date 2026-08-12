#include "src/trio/evolution/AncestralAdjacency.h"

#include "src/trio/model/StableId.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

void addAdjacencyPair(const std::string &leftTaxon,
                      const std::string &rightTaxon,
                      const std::string &sequence,
                      AlignmentEvidenceSet &evidence) {
    std::ostringstream maf;
    maf << "a score=10\n"
        << "s chr 0 " << sequence.size() << " + " << sequence.size()
        << ' ' << sequence << "\n"
        << "s chr 0 " << sequence.size() << " + " << sequence.size()
        << ' ' << sequence << "\n";
    std::istringstream inputStream(maf.str());
    PairwiseMafInput input;
    input.leftTaxon = leftTaxon;
    input.rightTaxon = rightTaxon;
    input.runId = leftTaxon + "_" + rightTaxon;
    input.mafPath = input.runId + ".maf";
    PairwiseEvidenceNormalizer::appendBlocks(
        input,
        MafEvidenceReader::read(inputStream, input.mafPath,
                                leftTaxon, rightTaxon),
        evidence);
}

AlignmentSiteGraph structuralGraph() {
    AlignmentEvidenceSet evidence;
    addAdjacencyPair("I1", "I2", "ACGTGC", evidence);
    addAdjacencyPair("I1", "O1", "ACGTGC", evidence);
    addAdjacencyPair("I2", "O1", "ACGTGC", evidence);
    addAdjacencyPair("I1", "O2", "ACGTGC", evidence);
    addAdjacencyPair("I2", "O2", "ACGTGC", evidence);
    AlignmentGraphBuildOptions options;
    options.ingroup1 = "I1";
    options.ingroup2 = "I2";
    options.primaryOutgroup = "O1";
    options.requireResolvedCopyAssignments = true;
    for (const auto &entry : evidence.residues) {
        CopyAssignment assignment;
        assignment.copyGroup = "resolved_chr_copy";
        assignment.hard = true;
        assignment.confidence = 1.0;
        assignment.provenance = "test";
        options.copyAssignments[entry.first] = assignment;
    }
    return AlignmentSiteGraphBuilder::build(evidence, options);
}

std::string pathForTaxon(const AlignmentSiteGraph &graph,
                         const std::string &taxon) {
    for (const PathSegment &path : graph.pathSegments()) {
        if (path.taxon == taxon) return path.occurrencePath;
    }
    return "";
}

ExtantBlockProjection projection(const AlignmentSiteGraph &graph,
                                 const std::string &taxon,
                                 const std::string &blockId,
                                 int64_t start0,
                                 BlockOrientation orientation) {
    ExtantBlockProjection result;
    result.taxon = taxon;
    result.occurrencePath = pathForTaxon(graph, taxon);
    result.blockId = blockId;
    result.start0 = start0;
    result.end0 = start0 + 1;
    result.orientation = orientation;
    result.projectionId = stableId(
        "test_projection",
        {taxon, result.occurrencePath, blockId, std::to_string(start0),
         blockOrientationName(orientation)});
    return result;
}

const AncestralAdjacencyCall *findCall(
    const AncestralAdjacencyReport &report,
    const std::string &block1,
    BlockEnd end1,
    const std::string &block2,
    BlockEnd end2) {
    OrientedBlockEndpoint first = {block1, end1};
    OrientedBlockEndpoint second = {block2, end2};
    if (second < first) std::swap(first, second);
    for (const AncestralAdjacencyCall &call : report.calls) {
        if (call.endpoint1 == first && call.endpoint2 == second) return &call;
    }
    return nullptr;
}

TEST(AncestralAdjacency, ProjectsIndependentBlocksAndLeavesGraphImmutable) {
    const AlignmentSiteGraph graph = structuralGraph();
    const std::string before = graph.coreHash();
    std::vector<std::string> canonicalSites;
    for (const PathSegment &path : graph.pathSegments()) {
        if (path.taxon == "I1") {
            canonicalSites = path.siteIds;
            break;
        }
    }
    ASSERT_EQ(6u, canonicalSites.size());
    std::vector<StructuralBlockDefinition> definitions(3);
    definitions[0].blockId = "B1";
    definitions[0].siteIds = {canonicalSites[0], canonicalSites[1]};
    definitions[1].blockId = "B2";
    definitions[1].siteIds = {canonicalSites[2], canonicalSites[3]};
    definitions[2].blockId = "B3";
    definitions[2].siteIds = {canonicalSites[4], canonicalSites[5]};

    const BlockProjectionSet projected =
        AncestralAdjacencyInference::projectBlocks(graph, definitions);
    EXPECT_EQ(before, graph.coreHash());
    ASSERT_EQ(12u, projected.projections.size());
    for (const ExtantBlockProjection &occurrence : projected.projections) {
        EXPECT_EQ(BlockOrientation::FORWARD, occurrence.orientation);
    }

    const Phylogeny tree =
        Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    AncestralAdjacencyConfig config;
    config.targetNode = "A";
    const AncestralAdjacencyReport report =
        AncestralAdjacencyInference::infer(graph, tree, projected, config);
    EXPECT_EQ(before, graph.coreHash());
    EXPECT_EQ("CANDIDATE_BLOCK_ADJACENCIES_ONLY", report.inferenceScope);
    ASSERT_EQ(2u, report.calls.size());
    const AncestralAdjacencyCall *first =
        findCall(report, "B1", BlockEnd::TAIL, "B2", BlockEnd::HEAD);
    ASSERT_NE(nullptr, first);
    EXPECT_FALSE(first->presence.isAmbiguous());
    EXPECT_EQ(1, first->presence.selectedState);
    EXPECT_TRUE(first->supportedWithoutLocalConflict);

    std::ostringstream projectionsTsv, callsTsv;
    AncestralAdjacencyWriters::writeProjections(projected, projectionsTsv);
    AncestralAdjacencyWriters::writeCalls(report, callsTsv);
    EXPECT_NE(std::string::npos,
              projectionsTsv.str().find("orientation\tcore_hash"));
    EXPECT_NE(std::string::npos,
              callsTsv.str().find("B1\tTAIL\tB2\tHEAD"));
    EXPECT_NE(std::string::npos,
              callsTsv.str().find("CANDIDATE_BLOCK_ADJACENCIES_ONLY"));
    EXPECT_NE(std::string::npos,
              callsTsv.str().find("NOT_ATTEMPTED"));
}

TEST(AncestralAdjacency, RepresentsReverseOrientationWithExplicitEndpoints) {
    const AlignmentSiteGraph graph = structuralGraph();
    BlockProjectionSet projected;
    projected.immutableCoreHash = graph.coreHash();
    for (const std::string &taxon : {"I1", "I2", "O1", "O2"}) {
        projected.projections.push_back(
            projection(graph, taxon, "B1", 0, BlockOrientation::FORWARD));
        projected.projections.push_back(
            projection(graph, taxon, "B2", 1, BlockOrientation::REVERSE));
    }
    const Phylogeny tree =
        Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    AncestralAdjacencyConfig config;
    config.targetNode = "A";
    const AncestralAdjacencyReport report =
        AncestralAdjacencyInference::infer(graph, tree, projected, config);
    ASSERT_EQ(1u, report.calls.size());
    const AncestralAdjacencyCall *call =
        findCall(report, "B1", BlockEnd::TAIL, "B2", BlockEnd::TAIL);
    ASSERT_NE(nullptr, call);
    EXPECT_EQ("SUPPORTED_CANDIDATE_ADJACENCY", call->decision);
    EXPECT_TRUE(call->supportedWithoutLocalConflict);
}

TEST(AncestralAdjacency, FlagsOrientationAlternativesWithoutBuildingChromosomes) {
    const AlignmentSiteGraph graph = structuralGraph();
    BlockProjectionSet projected;
    projected.immutableCoreHash = graph.coreHash();
    for (const std::string &taxon : {"I1", "I2"}) {
        projected.projections.push_back(
            projection(graph, taxon, "B1", 0, BlockOrientation::FORWARD));
        projected.projections.push_back(
            projection(graph, taxon, "B2", 1, BlockOrientation::FORWARD));
    }
    for (const std::string &taxon : {"O1", "O2"}) {
        projected.projections.push_back(
            projection(graph, taxon, "B1", 0, BlockOrientation::FORWARD));
        projected.projections.push_back(
            projection(graph, taxon, "B2", 1, BlockOrientation::REVERSE));
    }
    const Phylogeny tree =
        Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    AncestralAdjacencyConfig config;
    config.targetNode = "A";
    const AncestralAdjacencyReport report =
        AncestralAdjacencyInference::infer(graph, tree, projected, config);
    ASSERT_EQ(2u, report.calls.size());
    for (const AncestralAdjacencyCall &call : report.calls) {
        EXPECT_TRUE(call.orientationConflict);
        EXPECT_FALSE(call.supportedWithoutLocalConflict);
    }
    const AncestralAdjacencyCall *ingroupCall =
        findCall(report, "B1", BlockEnd::TAIL, "B2", BlockEnd::HEAD);
    ASSERT_NE(nullptr, ingroupCall);
    EXPECT_EQ(1, ingroupCall->presence.selectedState);
    EXPECT_EQ("SUPPORTED_ORIENTATION_CONFLICT", ingroupCall->decision);
}

TEST(AncestralAdjacency, AuditsAmbiguousEndpointDegreeAndMissingDuplicates) {
    const AlignmentSiteGraph graph = structuralGraph();
    const std::string immutable = graph.coreHash();
    BlockProjectionSet projected;
    projected.immutableCoreHash = immutable;
    projected.projections.push_back(
        projection(graph, "I1", "B1", 0, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "I1", "B2", 1, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "I1", "B3", 2, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "I2", "B2", 0, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "I2", "B1", 1, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "I2", "B3", 2, BlockOrientation::FORWARD));
    // O1 has a duplicated B1, so every candidate involving B1 is MISSING in O1.
    projected.projections.push_back(
        projection(graph, "O1", "B1", 0, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "O1", "B1", 2, BlockOrientation::FORWARD));
    projected.projections.push_back(
        projection(graph, "O1", "B2", 3, BlockOrientation::FORWARD));

    const Phylogeny tree =
        Phylogeny::fromNewick("((I1,I2)A,O1,O2)Root;");
    AncestralAdjacencyConfig config;
    config.targetNode = "A";
    const AncestralAdjacencyReport report =
        AncestralAdjacencyInference::infer(graph, tree, projected, config);
    EXPECT_EQ(immutable, graph.coreHash());
    const AncestralAdjacencyCall *left =
        findCall(report, "B1", BlockEnd::TAIL, "B2", BlockEnd::HEAD);
    const AncestralAdjacencyCall *right =
        findCall(report, "B1", BlockEnd::TAIL, "B3", BlockEnd::HEAD);
    ASSERT_NE(nullptr, left);
    ASSERT_NE(nullptr, right);
    EXPECT_TRUE(left->presence.isAmbiguous());
    EXPECT_TRUE(right->presence.isAmbiguous());
    EXPECT_TRUE(left->endpointConflict);
    EXPECT_TRUE(right->endpointConflict);
    EXPECT_EQ(PresenceObservation::MISSING, left->leafObservations.at("O1"));

    bool duplicateIssue = false;
    for (const BlockProjectionIssue &issue : report.issues) {
        duplicateIssue |=
            issue.kind == BlockProjectionIssueKind::DUPLICATE_BLOCK_COPY;
    }
    EXPECT_TRUE(duplicateIssue);
    std::ostringstream issues;
    AncestralAdjacencyWriters::writeIssues(report, issues);
    EXPECT_NE(std::string::npos, issues.str().find("DUPLICATE_BLOCK_COPY"));

    std::reverse(projected.projections.begin(), projected.projections.end());
    const AncestralAdjacencyReport reordered =
        AncestralAdjacencyInference::infer(graph, tree, projected, config);
    std::ostringstream firstTsv, secondTsv;
    AncestralAdjacencyWriters::writeCalls(report, firstTsv);
    AncestralAdjacencyWriters::writeCalls(reordered, secondTsv);
    EXPECT_EQ(firstTsv.str(), secondTsv.str());
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
