#include "src/service/ReadSamUtils.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

TEST(ReadSamMemory, RetainsOnlyExpectedCopiesPlusOneScores) {
    std::vector<double> scores;
    for (int value = 1; value <= 1000; ++value) {
        anchorwave::read_sam_detail::retainTopCopyScores(
                scores, static_cast<double>(value), 2);
    }

    ASSERT_EQ(3u, scores.size());
    std::sort(scores.begin(), scores.end());
    EXPECT_DOUBLE_EQ(998.0, scores[0]);
    EXPECT_DOUBLE_EQ(999.0, scores[1]);
    EXPECT_DOUBLE_EQ(1000.0, scores[2]);
}

TEST(ReadSamMemory, PreservesSecondaryCopyDecision) {
    std::vector<double> scores;
    for (double score : {0.98, 0.95, 0.91, 0.2, 0.1}) {
        anchorwave::read_sam_detail::retainTopCopyScores(scores, score, 2);
    }
    EXPECT_TRUE(anchorwave::read_sam_detail::secondaryCopyExceeds(
            scores, 2, 0.90));
    EXPECT_FALSE(anchorwave::read_sam_detail::secondaryCopyExceeds(
            scores, 2, 0.95));
}
