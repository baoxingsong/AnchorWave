#include "src/service/AlignmentGapPlanner.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <vector>

TEST(AlignmentGapPlannerTest, SelectsOnlyHighCostForwardInterAnchorGaps) {
    std::vector<anchorwave::AlignmentAnchorSpan> anchors{
            {100, 200, 1000, 1100, false},
            {40201, 40300, 41101, 41200, false},
            {40501, 40600, 41401, 41500, false}};

    const auto gaps = anchorwave::planParallelInterAnchorGaps(anchors, 6);
    ASSERT_EQ(1U, gaps.size());
    EXPECT_EQ(1U, gaps[0].anchorIndex);
    EXPECT_EQ(201U, gaps[0].referenceStart);
    EXPECT_EQ(40200U, gaps[0].referenceEnd);
    EXPECT_EQ(1101U, gaps[0].queryStart);
    EXPECT_EQ(41100U, gaps[0].queryEnd);
    EXPECT_FALSE(gaps[0].reverse);
    EXPECT_EQ(1600000000ULL, gaps[0].estimatedCost);
    EXPECT_EQ(1600000000ULL, gaps[0].geometricCost);
    EXPECT_EQ(0ULL, gaps[0].schedulingPriorityCost);
}

TEST(AlignmentGapPlannerTest, PlansReverseGapInAlignmentOrientation) {
    std::vector<anchorwave::AlignmentAnchorSpan> anchors{
            {100, 200, 90000, 90100, true},
            {40201, 40300, 49000, 50000, true}};

    const auto gaps = anchorwave::planParallelInterAnchorGaps(anchors, 4);
    ASSERT_EQ(1U, gaps.size());
    EXPECT_TRUE(gaps[0].reverse);
    EXPECT_EQ(201U, gaps[0].referenceStart);
    EXPECT_EQ(40200U, gaps[0].referenceEnd);
    EXPECT_EQ(50001U, gaps[0].queryStart);
    EXPECT_EQ(89999U, gaps[0].queryEnd);
}

TEST(AlignmentGapPlannerTest, SkipsMixedStrandsOverlapsAndSingleWorker) {
    std::vector<anchorwave::AlignmentAnchorSpan> mixed{
            {100, 200, 1000, 1100, false},
            {40201, 40300, 49000, 50000, true}};
    EXPECT_TRUE(anchorwave::planParallelInterAnchorGaps(mixed, 6).empty());

    std::vector<anchorwave::AlignmentAnchorSpan> overlap{
            {100, 300, 1000, 1200, false},
            {250, 50000, 1100, 50000, false}};
    EXPECT_TRUE(anchorwave::planParallelInterAnchorGaps(overlap, 6).empty());

    std::vector<anchorwave::AlignmentAnchorSpan> highCost{
            {100, 200, 1000, 1100, false},
            {40201, 40300, 41101, 41200, false}};
    EXPECT_TRUE(anchorwave::planParallelInterAnchorGaps(highCost, 1).empty());
}

TEST(AlignmentGapPlannerTest, UsesSaturatingCostAndDocumentedThreshold) {
    EXPECT_EQ(100000000ULL,
              anchorwave::minimumParallelAlignmentGapCost());
    EXPECT_TRUE(anchorwave::shouldParallelizeAlignmentGap(40000, 40000, 6));
    EXPECT_FALSE(anchorwave::shouldParallelizeAlignmentGap(1000, 1000, 6));
    EXPECT_FALSE(anchorwave::shouldParallelizeAlignmentGap(40000, 40000, 1));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::alignmentGapEstimatedCost(
                      std::numeric_limits<uint64_t>::max(), 2));
    EXPECT_EQ(1ULL, anchorwave::alignmentGapRuntimePriorityCost(0.0));
    EXPECT_EQ(60000000ULL,
              anchorwave::alignmentGapRuntimePriorityCost(1.0));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::alignmentGapRuntimePriorityCost(
                      std::numeric_limits<double>::max()));
}
