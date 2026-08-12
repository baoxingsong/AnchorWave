#include "src/trio/model/AlignmentEvidence.h"
#include "src/trio/model/StableId.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <map>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

std::vector<MafBlock> readFixture(const std::string &text,
                                  const std::string &left,
                                  const std::string &right) {
    std::istringstream input(text);
    return MafEvidenceReader::read(input, "fixture.maf", left, right);
}

TEST(AlignmentEvidence, NormalizesBaseAndGapEvidence) {
    PairwiseMafInput input;
    input.leftTaxon = "I1";
    input.rightTaxon = "O1";
    input.mafPath = "fixture.maf";
    input.runId = "run-i1-o1";
    const std::vector<MafBlock> blocks = readFixture(
        "a score=8\n"
        "s chr1 1 4 + 20 AC-GT\n"
        "s chrO 3 5 + 30 ACTGT\n",
        "I1", "O1");
    AlignmentEvidenceSet evidence;
    PairwiseEvidenceNormalizer::appendBlocks(input, blocks, evidence);
    EXPECT_EQ(4u, evidence.homologies.size());
    ASSERT_EQ(1u, evidence.alignedAbsences.size());
    EXPECT_EQ("I1", evidence.alignedAbsences[0].absentTaxon);
    EXPECT_EQ(2, evidence.alignedAbsences[0].absentLeftFlank0);
    EXPECT_EQ(3, evidence.alignedAbsences[0].absentRightFlank0);
    EXPECT_EQ(9u, evidence.residues.size());
}

TEST(AlignmentEvidence, EvidenceIdsAreInputOrientationIndependentAtPairLevel) {
    EXPECT_EQ(PairwiseEvidenceNormalizer::canonicalPairId("I1", "O1"),
              PairwiseEvidenceNormalizer::canonicalPairId("O1", "I1"));
}

TEST(AlignmentEvidence, NegativeStrandCoordinatesAreForwardNormalized) {
    PairwiseMafInput input;
    input.leftTaxon = "I1";
    input.rightTaxon = "I2";
    input.mafPath = "negative.maf";
    input.runId = "negative";
    const std::vector<MafBlock> blocks = readFixture(
        "a score=3\n"
        "s a 0 3 + 10 ACG\n"
        "s b 2 3 - 10 TGC\n",
        "I1", "I2");
    AlignmentEvidenceSet evidence;
    PairwiseEvidenceNormalizer::appendBlocks(input, blocks, evidence);
    ASSERT_EQ(3u, evidence.homologies.size());
    std::vector<int64_t> positions;
    for (const HomologyEvidence &edge : evidence.homologies) {
        positions.push_back(edge.right.forwardPosition0);
    }
    std::sort(positions.begin(), positions.end());
    EXPECT_EQ((std::vector<int64_t>{5, 6, 7}), positions);

    std::map<int64_t, char> forwardBases;
    for (const auto &entry : evidence.residues) {
        if (entry.first.taxon == "I2") {
            forwardBases[entry.first.forwardPosition0] = entry.second.base;
        }
    }
    EXPECT_EQ('G', forwardBases[5]);
    EXPECT_EQ('C', forwardBases[6]);
    EXPECT_EQ('A', forwardBases[7]);
}

TEST(AlignmentEvidence, SameResidueAgreesAcrossOppositeMafOrientations) {
    PairwiseMafInput forward;
    forward.leftTaxon = "I1";
    forward.rightTaxon = "I2";
    forward.mafPath = "forward.maf";
    forward.runId = "forward";
    PairwiseMafInput reverse = forward;
    reverse.rightTaxon = "O1";
    reverse.mafPath = "reverse.maf";
    reverse.runId = "reverse";

    AlignmentEvidenceSet evidence;
    PairwiseEvidenceNormalizer::appendBlocks(
        forward,
        readFixture("a\ns chr1 5 3 + 10 GCA\ns chr2 0 3 + 3 GCA\n", "I1", "I2"),
        evidence);
    EXPECT_NO_THROW(PairwiseEvidenceNormalizer::appendBlocks(
        reverse,
        readFixture("a\ns chr1 2 3 - 10 TGC\ns chrO 0 3 + 3 TGC\n", "I1", "O1"),
        evidence));

    ResidueId first;
    first.taxon = "I1";
    first.occurrencePath = "I1|chr1";
    first.sequence = "chr1";
    first.forwardPosition0 = 5;
    EXPECT_EQ('G', evidence.residues.at(first).base);
    first.forwardPosition0 = 6;
    EXPECT_EQ('C', evidence.residues.at(first).base);
    first.forwardPosition0 = 7;
    EXPECT_EQ('A', evidence.residues.at(first).base);
}

TEST(AlignmentEvidence, DetectsConflictingBaseForSameSourceResidue) {
    PairwiseMafInput first;
    first.leftTaxon = "I1";
    first.rightTaxon = "I2";
    first.mafPath = "first.maf";
    first.runId = "first";
    PairwiseMafInput second = first;
    second.rightTaxon = "O1";
    second.mafPath = "second.maf";
    second.runId = "second";
    AlignmentEvidenceSet evidence;
    PairwiseEvidenceNormalizer::appendBlocks(
        first, readFixture("a\ns chr1 0 1 + 2 A\ns chr2 0 1 + 2 A\n", "I1", "I2"),
        evidence);
    EXPECT_THROW(
        PairwiseEvidenceNormalizer::appendBlocks(
            second,
            readFixture("a\ns chr1 0 1 + 2 C\ns chrO 0 1 + 2 C\n", "I1", "O1"),
            evidence),
        AlignmentEvidenceError);
}

TEST(AlignmentEvidence, DetectsSourceSizeChangesWithoutCoordinateOverlap) {
    PairwiseMafInput first;
    first.leftTaxon = "I1";
    first.rightTaxon = "I2";
    first.mafPath = "first.maf";
    first.runId = "first";
    PairwiseMafInput second = first;
    second.rightTaxon = "O1";
    second.mafPath = "second.maf";
    second.runId = "second";
    AlignmentEvidenceSet evidence;
    PairwiseEvidenceNormalizer::appendBlocks(
        first,
        readFixture("a\ns chr1 0 1 + 10 A\ns chr2 0 1 + 2 A\n", "I1", "I2"),
        evidence);
    EXPECT_THROW(PairwiseEvidenceNormalizer::appendBlocks(
                     second,
                     readFixture("a\ns chr1 18 1 + 20 C\ns chrO 0 1 + 2 C\n",
                                 "I1", "O1"),
                     evidence),
                 AlignmentEvidenceError);
}

TEST(StableId, DelimitsFieldsAndIsDeterministic) {
    EXPECT_EQ(stableId("x", {"a", "bc"}), stableId("x", {"a", "bc"}));
    EXPECT_NE(stableId("x", {"a", "bc"}), stableId("x", {"ab", "c"}));
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
