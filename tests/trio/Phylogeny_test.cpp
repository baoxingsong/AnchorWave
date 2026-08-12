#include "gtest/gtest.h"

#include "../../src/trio/Phylogeny.h"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using anchorwave::trio::BinaryInferenceResult;
using anchorwave::trio::BinaryParsimonyModel;
using anchorwave::trio::NucleotideInferenceResult;
using anchorwave::trio::NucleotideParsimonyModel;
using anchorwave::trio::Phylogeny;
using anchorwave::trio::PresenceObservation;

Phylogeny::NodeId nodeNamed(const Phylogeny& tree, const std::string& name) {
    for (Phylogeny::NodeId node = 0; node < tree.nodeCount(); ++node) {
        if (tree.nodeName(node) == name) {
            return node;
        }
    }
    throw std::runtime_error("test helper did not find node: " + name);
}

TEST(PhylogenyNewick, ParsesNamesLengthsWhitespaceCommentsAndQuotes) {
    const Phylogeny tree = Phylogeny::fromNewick(
        " [tree comment] ( ('inner one':0.10, I2:2e-1) A : .30, "
        "O1:1.0, 'O''2' ) Root ; ");

    EXPECT_EQ(6u, tree.nodeCount());
    EXPECT_EQ(4u, tree.leafCount());
    EXPECT_EQ("Root", tree.nodeName(tree.root()));

    const Phylogeny::NodeId ancestor = tree.internalNodeId("A");
    EXPECT_TRUE(tree.isInternal(ancestor));
    ASSERT_TRUE(tree.hasBranchLength(ancestor));
    EXPECT_NEAR(0.30, tree.branchLength(ancestor), 1e-12);

    const Phylogeny::NodeId innerOne = nodeNamed(tree, "inner one");
    ASSERT_TRUE(tree.isLeaf(innerOne));
    EXPECT_NEAR(0.10, tree.branchLength(innerOne), 1e-12);
    EXPECT_EQ("O'2", tree.nodeName(nodeNamed(tree, "O'2")));
    EXPECT_FALSE(tree.hasBranchLength(nodeNamed(tree, "O'2")));
    EXPECT_THROW(tree.branchLength(nodeNamed(tree, "O'2")), std::logic_error);

    const std::vector<std::string> expectedLeaves = {"inner one", "I2", "O1", "O'2"};
    EXPECT_EQ(expectedLeaves, tree.leafNames());
}

TEST(PhylogenyNewick, AcceptsUnnamedInternalNodesAndOptionalSemicolon) {
    const Phylogeny tree = Phylogeny::fromNewick("((R,S),O)");
    EXPECT_EQ(5u, tree.nodeCount());
    EXPECT_TRUE(tree.nodeName(tree.root()).empty());
    EXPECT_EQ(3u, tree.leafCount());
}

TEST(PhylogenyNewick, RejectsMalformedTreesAndDuplicateLeaves) {
    const std::vector<std::string> invalidTrees = {
        "",
        "(A,B",
        "(A,,B)Root;",
        "(A,A)Root;",
        "(A:bad,B)Root;",
        "(A:-0.1,B)Root;",
        "(A,B)Root; trailing",
        "()Root;",
        "(A,B))Root;",
        "(A:[comment],B)Root;",
        "(A,B)[unterminated"
    };

    for (std::size_t i = 0; i < invalidTrees.size(); ++i) {
        EXPECT_THROW(Phylogeny::fromNewick(invalidTrees[i]), std::invalid_argument)
            << "tree was unexpectedly accepted: " << invalidTrees[i];
    }
}

TEST(PhylogenyNewick, InternalNodeLookupRejectsLeavesMissingAndAmbiguousNames) {
    const Phylogeny tree =
        Phylogeny::fromNewick("((R,S)A,(O1,O2)A)Root;");
    EXPECT_THROW(tree.internalNodeId("A"), std::invalid_argument);
    EXPECT_THROW(tree.internalNodeId("R"), std::invalid_argument);
    EXPECT_THROW(tree.internalNodeId("missing"), std::invalid_argument);
    EXPECT_THROW(tree.internalNodeId(""), std::invalid_argument);
}

TEST(PhylogenyNucleotide, MultiOutgroupsResolveTheIngroupAncestor) {
    const Phylogeny tree =
        Phylogeny::fromNewick("((R,S)A,O1,O2)Root;");
    const std::map<std::string, char> observations = {
        {"R", 'A'}, {"S", 'G'}, {"O1", 'G'}, {"O2", 'G'}
    };

    const NucleotideInferenceResult result =
        tree.inferNucleotide("A", observations);

    ASSERT_EQ(1u, result.optimalStates.size());
    EXPECT_EQ('G', result.optimalStates[0]);
    EXPECT_EQ('G', result.selectedState);
    EXPECT_EQ('G', result.ambiguityCode);
    EXPECT_DOUBLE_EQ(1.0, result.bestCost);
    EXPECT_DOUBLE_EQ(2.0, result.secondBestCost);
    EXPECT_DOUBLE_EQ(1.0, result.confidenceMargin);
    EXPECT_FALSE(result.isAmbiguous());
}

TEST(PhylogenyNucleotide, IsInvariantToTaxonAndChildOrder) {
    const Phylogeny first =
        Phylogeny::fromNewick("((R,S)A,O1,O2)Root;");
    const Phylogeny reordered =
        Phylogeny::fromNewick("(O2,(S,R)A,O1)Root;");

    std::map<std::string, char> firstOrder;
    firstOrder["R"] = 'A';
    firstOrder["S"] = 'G';
    firstOrder["O1"] = 'G';
    firstOrder["O2"] = 'G';

    std::map<std::string, char> anotherInsertionOrder;
    anotherInsertionOrder["O2"] = 'G';
    anotherInsertionOrder["O1"] = 'G';
    anotherInsertionOrder["S"] = 'G';
    anotherInsertionOrder["R"] = 'A';

    const NucleotideInferenceResult lhs =
        first.inferNucleotide("A", firstOrder);
    const NucleotideInferenceResult rhs =
        reordered.inferNucleotide("A", anotherInsertionOrder);

    EXPECT_EQ(lhs.optimalStates, rhs.optimalStates);
    EXPECT_EQ(lhs.selectedState, rhs.selectedState);
    EXPECT_EQ(lhs.ambiguityCode, rhs.ambiguityCode);
    EXPECT_EQ(lhs.stateCosts, rhs.stateCosts);
    EXPECT_DOUBLE_EQ(lhs.confidenceMargin, rhs.confidenceMargin);
}

TEST(PhylogenyNucleotide, MissingOutgroupProducesDeterministicIupacTie) {
    const Phylogeny tree = Phylogeny::fromNewick("((R,S)A,O)Root;");
    const std::map<std::string, char> observations = {
        {"R", 'a'}, {"S", 'G'}, {"O", '?'}
    };

    const NucleotideInferenceResult result =
        tree.inferNucleotide("A", observations);

    const std::vector<char> expected = {'A', 'G'};
    EXPECT_EQ(expected, result.optimalStates);
    EXPECT_EQ('A', result.selectedState);
    EXPECT_EQ('R', result.ambiguityCode);
    EXPECT_DOUBLE_EQ(0.0, result.confidenceMargin);
    EXPECT_TRUE(result.isAmbiguous());
}

TEST(PhylogenyNucleotide, SupportsIupacLeavesAndAllMissingEvidence) {
    const Phylogeny tree = Phylogeny::fromNewick("(R,S)A;");

    const NucleotideInferenceResult constrained = tree.inferNucleotide(
        "A", {{"R", 'R'}, {"S", 'A'}});
    EXPECT_EQ(std::vector<char>(1, 'A'), constrained.optimalStates);
    EXPECT_EQ('A', constrained.ambiguityCode);

    const NucleotideInferenceResult allMissing = tree.inferNucleotide(
        "A", {{"R", 'N'}, {"S", '.'}});
    const std::vector<char> allStates = {'A', 'C', 'G', 'T'};
    EXPECT_EQ(allStates, allMissing.optimalStates);
    EXPECT_EQ('A', allMissing.selectedState);
    EXPECT_EQ('N', allMissing.ambiguityCode);
    EXPECT_DOUBLE_EQ(0.0, allMissing.confidenceMargin);
}

TEST(PhylogenyNucleotide, RejectsUnknownTaxaSymbolsLeafTargetsAndInvalidCosts) {
    const Phylogeny tree = Phylogeny::fromNewick("(R,S)A;");
    EXPECT_THROW(tree.inferNucleotide("A", {{"unknown", 'A'}}),
                 std::invalid_argument);
    EXPECT_THROW(tree.inferNucleotide("A", {{"R", 'Z'}}),
                 std::invalid_argument);
    EXPECT_THROW(tree.inferNucleotide(nodeNamed(tree, "R"), {{"R", 'A'}}),
                 std::invalid_argument);

    NucleotideParsimonyModel invalid;
    invalid.transitionCost[0][1] = -1.0;
    EXPECT_THROW(tree.inferNucleotide("A", {{"R", 'A'}}, invalid),
                 std::invalid_argument);
}

TEST(PhylogenyPresence, MultiOutgroupsPolarizeADeletion) {
    const Phylogeny tree =
        Phylogeny::fromNewick("((R,S)A,O1,O2)Root;");
    const std::map<std::string, PresenceObservation> observations = {
        {"R", PresenceObservation::PRESENT},
        {"S", PresenceObservation::ABSENT},
        {"O1", PresenceObservation::PRESENT},
        {"O2", PresenceObservation::PRESENT}
    };

    const BinaryInferenceResult result =
        tree.inferPresence("A", observations);

    ASSERT_EQ(1u, result.optimalStates.size());
    EXPECT_EQ(1, result.optimalStates[0]);
    EXPECT_EQ(1, result.selectedState);
    EXPECT_DOUBLE_EQ(1.0, result.bestCost);
    EXPECT_DOUBLE_EQ(2.0, result.secondBestCost);
    EXPECT_DOUBLE_EQ(1.0, result.confidenceMargin);
    EXPECT_FALSE(result.isAmbiguous());
}

TEST(PhylogenyPresence, MissingEvidenceKeepsAnExplicitDeterministicTie) {
    const Phylogeny tree = Phylogeny::fromNewick("((R,S)A,O)Root;");
    const std::map<std::string, PresenceObservation> observations = {
        {"R", PresenceObservation::PRESENT},
        {"S", PresenceObservation::ABSENT},
        {"O", PresenceObservation::MISSING}
    };

    const BinaryInferenceResult result = tree.inferPresence("A", observations);
    const std::vector<int> expected = {0, 1};
    EXPECT_EQ(expected, result.optimalStates);
    EXPECT_EQ(0, result.selectedState);
    EXPECT_DOUBLE_EQ(0.0, result.confidenceMargin);
    EXPECT_TRUE(result.isAmbiguous());
}

TEST(PhylogenyPresence, GainAndLossCostsAreIndependentOfNucleotideSettings) {
    const Phylogeny tree = Phylogeny::fromNewick("(R,S)A;");
    const std::map<std::string, PresenceObservation> observations = {
        {"R", PresenceObservation::PRESENT},
        {"S", PresenceObservation::ABSENT}
    };

    const BinaryInferenceResult lossFavoured = tree.inferPresence(
        "A", observations, BinaryParsimonyModel::gainLoss(5.0, 1.0));
    const BinaryInferenceResult gainFavoured = tree.inferPresence(
        "A", observations, BinaryParsimonyModel::gainLoss(1.0, 5.0));

    EXPECT_EQ(1, lossFavoured.selectedState);
    EXPECT_EQ(0, gainFavoured.selectedState);
    EXPECT_FALSE(lossFavoured.isAmbiguous());
    EXPECT_FALSE(gainFavoured.isAmbiguous());

    NucleotideParsimonyModel nucleotideModel = NucleotideParsimonyModel::equalCost();
    nucleotideModel.transitionCost[0][2] = 100.0;
    nucleotideModel.transitionCost[2][0] = 100.0;
    const NucleotideInferenceResult nucleotide = tree.inferNucleotide(
        "A", {{"R", 'A'}, {"S", 'A'}}, nucleotideModel);
    EXPECT_EQ('A', nucleotide.selectedState);

    const BinaryInferenceResult repeated = tree.inferPresence(
        "A", observations, BinaryParsimonyModel::gainLoss(5.0, 1.0));
    EXPECT_EQ(lossFavoured.stateCosts, repeated.stateCosts);
    EXPECT_EQ(lossFavoured.selectedState, repeated.selectedState);
}

TEST(PhylogenyPresence, DirectedGainLossCostsRespectTheOriginalRooting) {
    const Phylogeny tree =
        Phylogeny::fromNewick("((R,S)A,O)Root;");
    const BinaryParsimonyModel model =
        BinaryParsimonyModel::gainLoss(7.0, 1.0);

    const BinaryInferenceResult result = tree.inferPresence(
        "A",
        {{"R", PresenceObservation::PRESENT},
         {"S", PresenceObservation::PRESENT},
         {"O", PresenceObservation::ABSENT}},
        model);

    EXPECT_EQ(1, result.selectedState);
    EXPECT_FALSE(result.isAmbiguous());
}

TEST(PhylogenyPresence, RejectsUnknownTaxaAndInvalidModelCosts) {
    const Phylogeny tree = Phylogeny::fromNewick("(R,S)A;");
    EXPECT_THROW(
        tree.inferPresence("A", {{"unknown", PresenceObservation::PRESENT}}),
        std::invalid_argument);

    BinaryParsimonyModel invalid;
    invalid.transitionCost[0][1] = -0.5;
    EXPECT_THROW(
        tree.inferPresence("A", {{"R", PresenceObservation::PRESENT}}, invalid),
        std::invalid_argument);
}

}  // namespace
