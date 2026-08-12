#include "src/trio/impl/LocalizedPoaRepair.h"

#include "gtest/gtest.h"

#include <algorithm>

namespace anchorwave {
namespace trio {
namespace {

TEST(TwoPieceAffine, UsesDeclaredOpenPlusLengthTimesExtendConvention) {
    TwoPieceAffineScoring scoring;
    scoring.gapOpen1 = -4;
    scoring.gapExtend1 = -2;
    scoring.gapOpen2 = -10;
    scoring.gapExtend2 = -1;
    EXPECT_DOUBLE_EQ(-6, scoring.gapScore(1));
    EXPECT_DOUBLE_EQ(-10, scoring.gapScore(3));
    EXPECT_DOUBLE_EQ(-16, scoring.gapScore(6));
}

TEST(TwoPieceAffine, GlobalAlignmentRoundTripsAndRescoresExactly) {
    TwoPieceAffineScoring scoring;
    const PairwiseAlignment alignment =
        TwoPieceAffineAligner::global("ACGTACGT", "ACGTTACGT", scoring);
    std::string first = alignment.first;
    std::string second = alignment.second;
    first.erase(std::remove(first.begin(), first.end(), '-'), first.end());
    second.erase(std::remove(second.begin(), second.end(), '-'), second.end());
    EXPECT_EQ("ACGTACGT", first);
    EXPECT_EQ("ACGTTACGT", second);
    EXPECT_DOUBLE_EQ(alignment.score,
                     TwoPieceAffineAligner::scoreAlignedPair(
                         alignment.first, alignment.second, scoring));
}

LocalizedPoaRequest requestFixture() {
    LocalizedPoaRequest request;
    request.conflictId = "conflict-1";
    request.immutableLeftSite = "left-anchor";
    request.immutableRightSite = "right-anchor";
    request.sequences = {{"I1", "ACGTT"}, {"I2", "ACTT"}, {"O1", "ACGCTT"}};
    return request;
}

TEST(LocalizedPoa, IsInvariantToInputSequenceOrder) {
    const LocalizedPoaOptions options;
    LocalizedPoaRequest first = requestFixture();
    LocalizedPoaRequest second = first;
    std::reverse(second.sequences.begin(), second.sequences.end());
    const LocalizedPoaPatch a = LocalizedPoaRepair::repair(first, options);
    const LocalizedPoaPatch b = LocalizedPoaRepair::repair(second, options);
    EXPECT_EQ(a.pathIds, b.pathIds);
    EXPECT_EQ(a.alignedRows, b.alignedRows);
    EXPECT_EQ(a.alignmentHash, b.alignmentHash);
    EXPECT_DOUBLE_EQ(a.repairedScore, b.repairedScore);
}

TEST(LocalizedPoa, EveryRowSpellsItsInputAndBoundariesAreImmutableMetadata) {
    const LocalizedPoaRequest request = requestFixture();
    const LocalizedPoaPatch patch =
        LocalizedPoaRepair::repair(request, LocalizedPoaOptions());
    EXPECT_EQ(request.immutableLeftSite, patch.immutableLeftSite);
    EXPECT_EQ(request.immutableRightSite, patch.immutableRightSite);
    for (std::size_t i = 0; i < patch.pathIds.size(); ++i) {
        std::string ungapped = patch.alignedRows[i];
        ungapped.erase(std::remove(ungapped.begin(), ungapped.end(), '-'), ungapped.end());
        const auto original = std::find_if(
            request.sequences.begin(), request.sequences.end(),
            [&](const PoaSequence &value) { return value.pathId == patch.pathIds[i]; });
        ASSERT_NE(request.sequences.end(), original);
        EXPECT_EQ(original->sequence, ungapped);
    }
}

TEST(LocalizedPoa, RejectsCandidateBelowBaselineDeltaWithoutMutatingAnything) {
    LocalizedPoaRequest request = requestFixture();
    const LocalizedPoaPatch initial =
        LocalizedPoaRepair::repair(request, LocalizedPoaOptions());
    for (std::size_t i = 0; i < initial.pathIds.size(); ++i) {
        request.baselineRows[initial.pathIds[i]] = initial.alignedRows[i];
    }
    LocalizedPoaOptions strict;
    strict.minimumScoreDelta = 1.0;
    const LocalizedPoaPatch rejected = LocalizedPoaRepair::repair(request, strict);
    EXPECT_FALSE(rejected.accepted);
    EXPECT_DOUBLE_EQ(0.0, rejected.scoreDelta);
    EXPECT_EQ("candidate_patch_below_minimum_delta", rejected.disposition);
}

TEST(LocalizedPoa, EnforcesWindowResourceLimits) {
    LocalizedPoaRequest request = requestFixture();
    LocalizedPoaOptions options;
    options.maxTotalBases = 2;
    EXPECT_THROW(LocalizedPoaRepair::repair(request, options), LocalizedPoaError);
}

TEST(LocalizedPoa, EnforcesPairwiseDpCellLimitBeforeAllocating) {
    LocalizedPoaRequest request = requestFixture();
    LocalizedPoaOptions options;
    options.maxPairwiseDpCells = 4;
    EXPECT_THROW(LocalizedPoaRepair::repair(request, options), LocalizedPoaError);
}

TEST(TwoPieceAffine, RejectsOversizedMatrixBeforeAllocating) {
    TwoPieceAffineScoring scoring;
    EXPECT_THROW(TwoPieceAffineAligner::global("ACGT", "ACGT", scoring, 24),
                 TwoPieceAffineError);
    EXPECT_NO_THROW(TwoPieceAffineAligner::global("ACGT", "ACGT", scoring, 25));
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
