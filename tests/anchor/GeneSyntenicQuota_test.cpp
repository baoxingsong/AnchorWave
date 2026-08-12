#include "src/impl/geneSyntenic.h"
#include "src/service/NovelAnchorThreshold.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

TEST(NovelAnchorThresholdTest, UsesTheIntended64BitSquare) {
    EXPECT_EQ(UINT64_C(10000000000),
              anchorwave::novelAnchorMinimumArea(100000));
    EXPECT_FALSE(anchorwave::novelAnchorAreaExceeds(
            100000, 100000, UINT64_C(10000000000)));
    EXPECT_TRUE(anchorwave::novelAnchorAreaExceeds(
            100001, 100000, UINT64_C(10000000000)));
}

TEST(NovelAnchorThresholdTest, RejectsInvalidOrOverflowingValues) {
    EXPECT_THROW(anchorwave::novelAnchorMinimumArea(0),
                 std::invalid_argument);
    EXPECT_THROW(anchorwave::novelAnchorMinimumArea(-1),
                 std::invalid_argument);
    EXPECT_THROW(anchorwave::novelAnchorMinimumArea(
                         std::numeric_limits<int64_t>::max()),
                 std::overflow_error);
}

TEST(LongestPathQuotaTest, StopsWhenTheLastAnchorIsConsumed) {
    std::vector<AlignmentMatch> matches;
    matches.emplace_back("ref", "query", 10, 19, 20, 29, 1.0,
                         POSITIVE, 0, 0, "gene", "gene");

    std::map<std::string, std::map<int64_t, AlignmentMatch>> refIndex;
    std::map<std::string, std::map<int64_t, AlignmentMatch>> queryIndex;
    refIndex["ref"][0] = matches.front();
    queryIndex["query"][0] = matches.front();

    std::vector<std::vector<AlignmentMatch>> chains;
    double indelScore = -0.01;
    double gapOpenPenalty = -0.03;
    double minimumAlignmentScore = 0.5;
    const int maximumDistance = 25;
    int referenceMaximumTimes = 1;
    int queryMaximumTimes = 1;
    double indelDistance = 3.0;

    EXPECT_NO_THROW(longestPathQuotav2(
            matches, chains, refIndex, queryIndex, indelScore,
            gapOpenPenalty, minimumAlignmentScore, maximumDistance,
            referenceMaximumTimes, queryMaximumTimes, indelDistance, false, 2));
    ASSERT_EQ(1U, chains.size());
    ASSERT_EQ(1U, chains.front().size());
    EXPECT_EQ("gene", chains.front().front().getReferenceGeneName());
}
