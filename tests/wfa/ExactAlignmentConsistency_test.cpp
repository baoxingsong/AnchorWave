#include "src/myImportandFunction/WfaAlignment.h"
#include "src/myImportandFunction/AlignmentAlgorithmSelector.h"
#include "src/myImportandFunction/AlignmentResourceScheduler.h"
#include "src/myImportandFunction/alignSlidingWindow.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

// Legacy public symbol used by approximate-tier windowing. The implementation now
// performs one semiglobal traceback pass; this test compares it against the
// former two-pass definition built from independent exact prefix alignments.
int32_t alignment_minimap2(
        const std::string &query, const std::string &reference,
        std::string &queryAlignment, std::string &referenceAlignment,
        const int32_t &matchingScore, int32_t mismatchingPenalty,
        int32_t openGapPenalty1, int32_t extendGapPenalty1,
        int32_t openGapPenalty2, int32_t extendGapPenalty2,
        int32_t &endPositionQuery, int32_t &endPositionReference);

namespace {

constexpr uint64_t kTestMemoryBudget = 1024ULL * 1024ULL * 1024ULL;

std::string removeGaps(std::string sequence) {
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '-'),
                   sequence.end());
    return sequence;
}

int64_t scoreAlignment(const std::string &query,
                       const std::string &reference) {
    EXPECT_EQ(query.size(), reference.size());
    int64_t score = 0;
    std::size_t column = 0;
    while (column < query.size()) {
        const bool queryGap = query[column] == '-';
        const bool referenceGap = reference[column] == '-';
        EXPECT_FALSE(queryGap && referenceGap);
        if (queryGap || referenceGap) {
            const bool gapInQuery = queryGap;
            std::size_t end = column + 1;
            while (end < query.size() &&
                   (query[end] == '-') == gapInQuery &&
                   (reference[end] == '-') != gapInQuery) {
                ++end;
            }
            const int64_t length = static_cast<int64_t>(end - column);
            score -= std::min<int64_t>(8 + 2 * length, 75 + length);
            column = end;
            continue;
        }
        if (query[column] != reference[column]) {
            score -= 6;
        }
        ++column;
    }
    return score;
}

struct SemiglobalReference {
    int64_t score = KSW_NEG_INF;
    int32_t queryEnd = 0;
    int32_t referenceEnd = 0;
    std::string queryAlignment;
    std::string referenceAlignment;
};

SemiglobalReference twoPassSemiglobalReference(
        const std::string &query, const std::string &reference) {
    SemiglobalReference queryEndBest;
    SemiglobalReference referenceEndBest;
    for (std::size_t referenceEnd = 1;
         referenceEnd <= reference.size(); ++referenceEnd) {
        std::string queryAlignment;
        std::string referenceAlignment;
        const std::string prefix = reference.substr(0, referenceEnd);
        const int64_t score = alignSlidingWindow_minimap2(
                query, prefix, query.size(), prefix.size(),
                queryAlignment, referenceAlignment, -1,
                -6, -8, -2, -75, -1);
        // extd2 updates mqe on equality, retaining the latest endpoint.
        if (score >= queryEndBest.score) {
            queryEndBest = SemiglobalReference{
                    score, static_cast<int32_t>(query.size()),
                    static_cast<int32_t>(referenceEnd),
                    std::move(queryAlignment),
                    std::move(referenceAlignment)};
        }
    }
    for (std::size_t queryEnd = 1; queryEnd <= query.size(); ++queryEnd) {
        std::string queryAlignment;
        std::string referenceAlignment;
        const std::string prefix = query.substr(0, queryEnd);
        const int64_t score = alignSlidingWindow_minimap2(
                prefix, reference, prefix.size(), reference.size(),
                queryAlignment, referenceAlignment, -1,
                -6, -8, -2, -75, -1);
        // extd2 updates mte on equality, retaining the latest endpoint.
        if (score >= referenceEndBest.score) {
            referenceEndBest = SemiglobalReference{
                    score, static_cast<int32_t>(queryEnd),
                    static_cast<int32_t>(reference.size()),
                    std::move(queryAlignment),
                    std::move(referenceAlignment)};
        }
    }
    // The production code deliberately resolves an mqe/mte tie toward mte.
    return queryEndBest.score > referenceEndBest.score
           ? queryEndBest : referenceEndBest;
}

void expectSameExactScore(const std::string &query,
                          const std::string &reference) {
    std::string kswQuery;
    std::string kswReference;
    const int64_t kswScore = alignSlidingWindow_minimap2(
            query, reference,
            static_cast<int64_t>(query.size()),
            static_cast<int64_t>(reference.size()),
            kswQuery, kswReference, -1,
            -6, -8, -2, -75, -1);
    EXPECT_EQ(query, removeGaps(kswQuery));
    EXPECT_EQ(reference, removeGaps(kswReference));
    EXPECT_EQ(kswScore, scoreAlignment(kswQuery, kswReference));

    anchorwave::WfaAlignmentResult singletrack;
    anchorwave::WfaAlignmentResult high;
    anchorwave::WfaAlignmentResult medium;
    anchorwave::WfaAlignmentResult low;
    anchorwave::WfaAlignmentResult bidirectional;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithSingletrackWfa(
                      query, reference, kTestMemoryBudget,
                      -6, -8, -2, -75, -1, singletrack));
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      query, reference, kTestMemoryBudget,
                      -6, -8, -2, -75, -1, high));
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithMediumWfa(
                      query, reference, kTestMemoryBudget,
                      -6, -8, -2, -75, -1, medium));
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithLowWfa(
                      query, reference, kTestMemoryBudget,
                      -6, -8, -2, -75, -1, low));
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithBiWfa(
                      query, reference, kTestMemoryBudget,
                      -6, -8, -2, -75, -1, bidirectional));

    EXPECT_EQ(kswScore, singletrack.score);
    EXPECT_EQ(kswScore, high.score);
    EXPECT_EQ(kswScore, medium.score);
    EXPECT_EQ(kswScore, low.score);
    EXPECT_EQ(kswScore, bidirectional.score);
    EXPECT_EQ(singletrack.score,
              scoreAlignment(singletrack.queryAlignment,
                             singletrack.referenceAlignment));
    EXPECT_EQ(high.score,
              scoreAlignment(high.queryAlignment, high.referenceAlignment));
    EXPECT_EQ(medium.score,
              scoreAlignment(medium.queryAlignment,
                             medium.referenceAlignment));
    EXPECT_EQ(low.score,
              scoreAlignment(low.queryAlignment, low.referenceAlignment));
    EXPECT_EQ(bidirectional.score,
              scoreAlignment(bidirectional.queryAlignment,
                             bidirectional.referenceAlignment));
}

TEST(ExactAlignmentConsistency, AllCandidatesUseTheSameTwoPieceAffineScore) {
    std::mt19937 generator(20260809);
    std::uniform_int_distribution<int> base(0, 3);
    std::uniform_int_distribution<int> mutation(0, 99);
    constexpr char bases[] = "ACGT";

    for (int example = 0; example < 100; ++example) {
        const int referenceLength = 40 + example * 3;
        std::string reference;
        reference.reserve(referenceLength);
        for (int i = 0; i < referenceLength; ++i) {
            reference.push_back(bases[base(generator)]);
        }

        std::string query;
        query.reserve(reference.size() + 64);
        for (char nucleotide : reference) {
            const int action = mutation(generator);
            if (action < 7) {
                continue;
            }
            if (action < 15) {
                const int insertionLength = 1 + mutation(generator) % 12;
                for (int inserted = 0; inserted < insertionLength;
                     ++inserted) {
                    query.push_back(bases[base(generator)]);
                }
            }
            query.push_back(action >= 15 && action < 30
                            ? bases[base(generator)] : nucleotide);
        }
        expectSameExactScore(query, reference);
    }
}

TEST(ExactAlignmentConsistency,
     BandedKswReportsTheScoreOfItsEmittedEndToEndAlignment) {
    std::mt19937 generator(20260812);
    std::uniform_int_distribution<int> base(0, 3);
    constexpr char bases[] = "ACGT";
    std::string reference;
    reference.reserve(2400);
    for (int position = 0; position < 2400; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference.substr(0, 800);
    for (int position = 0; position < 500; ++position) {
        query.push_back(bases[base(generator)]);
    }
    query.append(reference, 800, 1100);

    std::string fullQuery;
    std::string fullReference;
    const int64_t fullScore = alignSlidingWindow_minimap2(
            query, reference, query.size(), reference.size(),
            fullQuery, fullReference, -1, -6, -8, -2, -75, -1);
    std::string bandedQuery;
    std::string bandedReference;
    const int64_t bandedScore = alignSlidingWindow_minimap2(
            query, reference, query.size(), reference.size(),
            bandedQuery, bandedReference, 500,
            -6, -8, -2, -75, -1);

    EXPECT_EQ(query, removeGaps(bandedQuery));
    EXPECT_EQ(reference, removeGaps(bandedReference));
    EXPECT_EQ(bandedScore, scoreAlignment(bandedQuery, bandedReference));
    EXPECT_LE(bandedScore, fullScore);
}

TEST(ExactAlignmentConsistency,
     SinglePassSemiglobalWindowMatchesTheHistoricalTwoPassDefinition) {
    std::mt19937 generator(20260810);
    std::uniform_int_distribution<int> base(0, 3);
    std::uniform_int_distribution<int> mutation(0, 9);
    constexpr char bases[] = "ACGT";

    for (int example = 0; example < 20; ++example) {
        std::string reference;
        const int referenceLength = 18 + example;
        reference.reserve(referenceLength);
        for (int position = 0; position < referenceLength; ++position) {
            reference.push_back(bases[base(generator)]);
        }
        std::string query = reference;
        if (example % 3 == 0) {
            query.insert(query.begin() + query.size() / 2,
                         1 + example % 4, 'A');
        } else if (example % 3 == 1 && query.size() > 4) {
            query.erase(query.begin() + query.size() / 3,
                        query.begin() + query.size() / 3 + 1 + example % 3);
        }
        for (char &character : query) {
            if (mutation(generator) == 0) {
                character = bases[base(generator)];
            }
        }

        const SemiglobalReference expected =
                twoPassSemiglobalReference(query, reference);
        std::string queryAlignment;
        std::string referenceAlignment;
        int32_t queryEnd = 0;
        int32_t referenceEnd = 0;
        const int32_t score = alignment_minimap2(
                query, reference, queryAlignment, referenceAlignment,
                1, -6, -8, -2, -75, -1,
                queryEnd, referenceEnd);

        EXPECT_EQ(expected.score, score) << "example " << example;
        EXPECT_EQ(expected.queryEnd, queryEnd) << "example " << example;
        EXPECT_EQ(expected.referenceEnd, referenceEnd)
                << "example " << example;
        EXPECT_EQ(expected.queryAlignment, queryAlignment)
                << "example " << example;
        EXPECT_EQ(expected.referenceAlignment, referenceAlignment)
                << "example " << example;
    }
}

TEST(ExactAlignmentConsistency, HandlesStrongLengthImbalance) {
    std::string query;
    query.reserve(1800);
    for (int i = 0; i < 450; ++i) {
        query += "ACGT";
    }
    const std::string reference =
            "ACGTACGTACGTACGTACGTACGTACGTACGT";
    expectSameExactScore(query, reference);
}

TEST(ExactAlignmentConsistency,
     AmbiguousBasesUseTheSameStrictObjectiveInEveryExactEngine) {
    const std::string reference =
            "ACGTTGCAACCGGTTAACCGGTTACGTACGTAACCGGTTA";
    const std::string query =
            "ACGTNNNNACCGGTTAACNGGTTACGTACGTNNCCGGTTA";
    expectSameExactScore(query, reference);
}

TEST(ExactAlignmentConsistency,
     ScoreCertifiedKsw2MatchesTheUnbandedGlobalOptimum) {
    std::mt19937 generator(20260810);
    std::uniform_int_distribution<int> base(0, 3);
    constexpr char bases[] = "ACGT";
    for (int example = 0; example < 64; ++example) {
        std::string reference;
        reference.reserve(320);
        for (int position = 0; position < 240 + example; ++position) {
            reference.push_back(bases[base(generator)]);
        }
        std::string query = reference;
        query.insert(40 + example % 31, 12 + example % 23, 'A');
        query.erase(150 + example % 41, 7 + example % 17);
        for (std::size_t position = 5; position < query.size(); position += 37) {
            query[position] = query[position] == 'A' ? 'C' : 'A';
        }

        std::string fullQuery;
        std::string fullReference;
        const int64_t fullScore = alignSlidingWindow_minimap2(
                query, reference,
                static_cast<int64_t>(query.size()),
                static_cast<int64_t>(reference.size()),
                fullQuery, fullReference, -1,
                -6, -8, -2, -75, -1);
        std::string certifiedQuery;
        std::string certifiedReference;
        const auto certified = alignScoreCertifiedKsw2(
                query, reference, certifiedQuery, certifiedReference,
                4, 256, 0.0, -6, -8, -2, -75, -1);
        ASSERT_EQ(anchorwave::Ksw2CertifiedStatus::Completed,
                  certified.status) << "example=" << example;
        EXPECT_EQ(fullScore, certified.optimalScore);
        EXPECT_EQ(fullScore, certified.score);
        EXPECT_EQ(fullScore,
                  scoreAlignment(certifiedQuery, certifiedReference));
        EXPECT_EQ(query, removeGaps(certifiedQuery));
        EXPECT_EQ(reference, removeGaps(certifiedReference));
    }
}

TEST(ExactAlignmentConsistency,
     ScoreCertifiedKsw2SkipsTheSmallestBandAndNeverEmitsUncertifiedAlignment) {
    const std::string reference = std::string(80, 'A') +
                                  std::string(80, 'C');
    const std::string query = std::string(40, 'G') +
                              std::string(80, 'A') +
                              std::string(40, 'C');

    std::string fullQuery;
    std::string fullReference;
    const int64_t fullScore = alignSlidingWindow_minimap2(
            query, reference,
            static_cast<int64_t>(query.size()),
            static_cast<int64_t>(reference.size()),
            fullQuery, fullReference, -1,
            -6, -8, -2, -75, -1);

    std::string certifiedQuery;
    std::string certifiedReference;
    const auto expanded = alignScoreCertifiedKsw2(
            query, reference, certifiedQuery, certifiedReference,
            2, 64, 0.0, -6, -8, -2, -75, -1);
    ASSERT_EQ(anchorwave::Ksw2CertifiedStatus::Completed, expanded.status);
    EXPECT_EQ(fullScore, expanded.score);
    EXPECT_GT(expanded.tracebackAttempts, 1ULL);
    EXPECT_GE(expanded.finalBandWidth, 40);

    certifiedQuery = "must be cleared";
    certifiedReference = "must be cleared";
    const auto insufficient = alignScoreCertifiedKsw2(
            query, reference, certifiedQuery, certifiedReference,
            2, 20, 0.0, -6, -8, -2, -75, -1);
    EXPECT_EQ(anchorwave::Ksw2CertifiedStatus::NotCertified,
              insufficient.status);
    EXPECT_GT(insufficient.tracebackAttempts, 1ULL);
    EXPECT_EQ(20, insufficient.finalBandWidth);
    EXPECT_TRUE(certifiedQuery.empty());
    EXPECT_TRUE(certifiedReference.empty());
    EXPECT_EQ(query,
              removeGaps(insufficient.bestEffortQueryAlignment));
    EXPECT_EQ(reference,
              removeGaps(insufficient.bestEffortReferenceAlignment));
    EXPECT_EQ(insufficient.score,
              scoreAlignment(insufficient.bestEffortQueryAlignment,
                             insufficient.bestEffortReferenceAlignment));
}

TEST(ExactAlignmentConsistency,
     ScoreCertifiedKsw2UsesExactAnchorChainAcrossLargeDiagonalDrift) {
    std::mt19937 generator(20260812);
    const char bases[] = {'A', 'C', 'G', 'T'};
    const auto randomSequence = [&](std::size_t length) {
        std::string sequence;
        sequence.reserve(length);
        for (std::size_t i = 0; i < length; ++i) {
            sequence.push_back(bases[generator() % 4]);
        }
        return sequence;
    };

    std::string query;
    std::string reference;
    anchorwave::AlignmentSelectionPlan plan;
    for (int block = 0; block <= 20; ++block) {
        const std::string shared = randomSequence(64);
        plan.certificationAnchors.push_back(
                anchorwave::AlignmentChainAnchor{
                        static_cast<uint32_t>(query.size()),
                        static_cast<uint32_t>(reference.size()), 15});
        query += shared;
        reference += shared;
        if (block < 10) {
            query += randomSequence(20);
        } else if (block < 20) {
            reference += randomSequence(20);
        }
    }

    std::string fullQuery;
    std::string fullReference;
    const int64_t fullScore = alignSlidingWindow_minimap2(
            query, reference,
            static_cast<int64_t>(query.size()),
            static_cast<int64_t>(reference.size()),
            fullQuery, fullReference, -1,
            -6, -8, -2, -75, -1);

    std::string fixedQuery;
    std::string fixedReference;
    const auto fixed = alignScoreCertifiedKsw2(
            query, reference, fixedQuery, fixedReference,
            32, 32, 0.0, -6, -8, -2, -75, -1);
    EXPECT_EQ(anchorwave::Ksw2CertifiedStatus::NotCertified, fixed.status);

    std::string chainedQuery;
    std::string chainedReference;
    const auto chained = alignScoreCertifiedKsw2(
            query, reference, chainedQuery, chainedReference,
            32, 32, 0.0, -6, -8, -2, -75, -1, &plan);
    ASSERT_EQ(anchorwave::Ksw2CertifiedStatus::Completed, chained.status);
    EXPECT_EQ(-2, chained.finalBandWidth);
    EXPECT_EQ(fullScore, chained.score);
    EXPECT_EQ(query, removeGaps(chainedQuery));
    EXPECT_EQ(reference, removeGaps(chainedReference));
    EXPECT_EQ(fullScore, scoreAlignment(chainedQuery, chainedReference));
}

TEST(ExactAlignmentConsistency,
     ScoreCertifiedKsw2RuntimeLimitReturnsNoPartialAlignment) {
    const std::string query(4096, 'A');
    const std::string reference(4096, 'T');
    std::string alignedQuery = "must be cleared";
    std::string alignedReference = "must be cleared";
    const auto result = alignScoreCertifiedKsw2(
            query, reference, alignedQuery, alignedReference,
            64, 4096, 1e-9, -6, -8, -2, -75, -1);
    EXPECT_EQ(anchorwave::Ksw2CertifiedStatus::TimeLimit, result.status);
    EXPECT_TRUE(alignedQuery.empty());
    EXPECT_TRUE(alignedReference.empty());
}

TEST(ExactAlignmentConsistency,
     SlidingWindowReportsTheScoreOfItsEmittedEndToEndAlignment) {
    std::string query =
            "ACGTACGTACGTGGGGTTTTACGTACGTACGTACGTACGT";
    std::string reference =
            "ACGTACGTACGTAAAACCCCTTTTACGTACGTACGT";
    std::string alignedQuery;
    std::string alignedReference;
    const int64_t score = alignSlidingWindowNW(
            query, reference, alignedQuery, alignedReference, 8,
            2, -6, -8, -2, -75, -1);

    EXPECT_EQ(query, removeGaps(alignedQuery));
    EXPECT_EQ(reference, removeGaps(alignedReference));
    EXPECT_EQ(scoreAlignment(alignedQuery, alignedReference), score);
}

TEST(ExactAlignmentConsistency, ExecutorChoosesFastestPredictedExactMode) {
    std::string query = "ACGTTGCAACCGGTTAACCGGTTACGT";
    std::string reference = "ACGTTGCAACCGGTTAACCGGTTACGT";
    std::string queryAlignment;
    std::string referenceAlignment;
    std::string method;
    const int64_t score = alignSlidingWindow(
            query, reference, queryAlignment, referenceAlignment, method,
            100000, 2, -6, -8, -2, -75, -1);

    EXPECT_EQ(0, score);
    // Tiny full matrices are cheaper in KSW2 than paying WFA setup cost.  The
    // quality tier remains exact even though the old fixed-order test expected
    // Singletrack for every interval.
    EXPECT_EQ("MINIMAP2", method);
    EXPECT_EQ(query, queryAlignment);
    EXPECT_EQ(reference, referenceAlignment);
}

TEST(ExactAlignmentConsistency,
     ReusingPrecomputedSelectionPlanPreservesTheAlignment) {
    const std::string originalReference =
            "ACGTTGCAACCGGTTAACCGGTTACGTACGTAACCGGTTACGTACGTA";
    const std::string originalQuery =
            "ACGTTGCAACCGGCTAATCCGGTTACGTACGTAACCGGTTTCGTACGTA";
    const auto plan = prepareAlignmentSelectionPlan(
            originalQuery, originalReference, 100000,
            -6, -8, -2, -75, -1);
    ASSERT_NE(nullptr, plan);

    std::string regularQuery = originalQuery;
    std::string regularReference = originalReference;
    std::string regularQueryAlignment;
    std::string regularReferenceAlignment;
    std::string regularMethod;
    const int64_t regularScore = alignSlidingWindow(
            regularQuery, regularReference,
            regularQueryAlignment, regularReferenceAlignment, regularMethod,
            100000, 2, -6, -8, -2, -75, -1);

    std::string reusedQuery = originalQuery;
    std::string reusedReference = originalReference;
    std::string reusedQueryAlignment;
    std::string reusedReferenceAlignment;
    std::string reusedMethod;
    const int64_t reusedScore = alignSlidingWindow(
            reusedQuery, reusedReference,
            reusedQueryAlignment, reusedReferenceAlignment, reusedMethod,
            100000, 2, -6, -8, -2, -75, -1, plan);

    EXPECT_EQ(regularScore, reusedScore);
    EXPECT_EQ(regularMethod, reusedMethod);
    EXPECT_EQ(regularQueryAlignment, reusedQueryAlignment);
    EXPECT_EQ(regularReferenceAlignment, reusedReferenceAlignment);
    EXPECT_EQ(originalQuery, removeGaps(reusedQueryAlignment));
    EXPECT_EQ(originalReference, removeGaps(reusedReferenceAlignment));
}

TEST(ExactAlignmentConsistency,
     ApproximateTierExecutesTheSelectorPrediction) {
    struct ExactLimitGuard {
        double previous =
                anchorwave::configuredExactAlignmentMaximumEstimatedMinutes();
        ~ExactLimitGuard() {
            anchorwave::configureExactAlignmentMaximumEstimatedMinutes(
                    previous);
        }
    } restoreExactLimit;
    constexpr double exactLimitMinutes = 1e-9;
    anchorwave::configureExactAlignmentMaximumEstimatedMinutes(
            exactLimitMinutes);

    std::mt19937 generator(20260812);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string reference;
    reference.reserve(4000);
    for (int position = 0; position < 4000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference.substr(0, 1600);
    for (int position = 0; position < 800; ++position) {
        query.push_back(bases[base(generator)]);
    }
    query.append(reference, 1600, 1600);
    ASSERT_EQ(reference.size(), query.size());

    int64_t windowWidth = 0;
    std::shared_ptr<const anchorwave::AlignmentSelectionPlan> plan;
    for (const int64_t candidateWidth :
         {250, 300, 400, 500, 650, 800, 1000, 1250, 1500,
          1750, 2000, 2500, 3000, 3500, 4000}) {
        auto candidatePlan =
                std::make_shared<const anchorwave::AlignmentSelectionPlan>(
                        anchorwave::makeAlignmentSelectionPlan(
                                query, reference, candidateWidth,
                                -6, -8, -2, -75, -1,
                                exactLimitMinutes));
        if (candidatePlan->exactCandidates.empty() &&
            candidatePlan->approximateCandidates.size() == 2) {
            windowWidth = candidateWidth;
            plan = std::move(candidatePlan);
            break;
        }
    }
    ASSERT_NE(nullptr, plan);
    ASSERT_EQ(2ULL, plan->approximateCandidates.size());
    ASSERT_GT(plan->ksw2BandWidth, 0);
    ASSERT_LT(plan->ksw2BandWidth,
              static_cast<int64_t>(std::max(query.size(), reference.size())));
    ASSERT_GT(plan->slidingWindowWidth, 0);

    std::string bandedQuery;
    std::string bandedReference;
    const int64_t bandedScore = alignSlidingWindow_minimap2(
            query, reference, query.size(), reference.size(),
            bandedQuery, bandedReference, plan->ksw2BandWidth,
            -6, -8, -2, -75, -1);

    std::string slidingQueryInput = query;
    std::string slidingReferenceInput = reference;
    std::string slidingQuery;
    std::string slidingReference;
    const int64_t slidingScore = alignSlidingWindowNW(
            slidingQueryInput, slidingReferenceInput,
            slidingQuery, slidingReference, plan->slidingWindowWidth,
            2, -6, -8, -2, -75, -1);
    ASSERT_NE(bandedScore, slidingScore);

    std::string selectedQueryInput = query;
    std::string selectedReferenceInput = reference;
    std::string selectedQuery;
    std::string selectedReference;
    std::string selectedMethod;
    const int64_t selectedScore = alignSlidingWindow(
            selectedQueryInput, selectedReferenceInput,
            selectedQuery, selectedReference, selectedMethod,
            windowWidth, 2, -6, -8, -2, -75, -1, plan);

    const bool predictedBanded = plan->approximateCandidates.front() ==
            anchorwave::AlignmentCandidate::Ksw2Banded;
    EXPECT_EQ(predictedBanded ? bandedScore : slidingScore, selectedScore);
    EXPECT_EQ(query, removeGaps(selectedQuery));
    EXPECT_EQ(reference, removeGaps(selectedReference));
    EXPECT_EQ(selectedScore, scoreAlignment(selectedQuery,
                                             selectedReference));
    EXPECT_EQ(predictedBanded ? "BANDED_MINIMAP2" : "SLIDING_WINDOW",
              selectedMethod);
}

TEST(ExactAlignmentConsistency,
     ProcessLimitDoesNotLetHighWfaBorrowBeyondWindowSquare) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const uint64_t baseline = gib / 4;
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            1, 3 * gib, 1000, baseline);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [baseline]() { return baseline; });
    anchorwave::ScopedAlignmentMemoryScheduler installed(scheduler);

    // Even a large -M cannot expand a single exact attempt beyond w^2.
    std::string query(20000, 'A');
    std::string reference = query;
    std::string queryAlignment;
    std::string referenceAlignment;
    std::string method;
    alignSlidingWindow(
            query, reference, queryAlignment, referenceAlignment, method,
            1000, 2, -6, -8, -2, -75, -1);

    EXPECT_NE("WAVEFRONT_SINGLETRACK", method);
    EXPECT_NE("WAVEFRONT", method);
    EXPECT_EQ(query, removeGaps(queryAlignment));
    EXPECT_EQ(reference, removeGaps(referenceAlignment));
    const auto stats = scheduler.stats();
    EXPECT_EQ(1ULL, stats.reservationCount);
    EXPECT_GT(stats.peakReservedBytes,
              anchorwave::wfaMemoryBudgetBytes(1000));
    EXPECT_LE(stats.peakProjectedProcessBytes + plan.safetyReserveBytes,
              plan.maxProcessMemoryBytes);
}

TEST(ExactAlignmentConsistency,
     TemporaryExactPressureDoesNotOpenTheBandedTier) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            1, 3 * gib, 4500, 0);
    // Leave less than one KSW/WFA reservation after the analytic resident
    // model; the exact tier must park/fail instead of silently opening a
    // lower-quality tier.
    const uint64_t available = 1ULL * 1024ULL * 1024ULL;
    const uint64_t usable = plan.maxProcessMemoryBytes -
                            plan.safetyReserveBytes;
    ASSERT_GT(usable, available);
    const uint64_t observed = usable - available;
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [observed]() { return observed; });
    anchorwave::ScopedAlignmentMemoryScheduler installed(scheduler);

    std::string query(1500, 'A');
    std::string reference(1600, 'T');
    std::string queryAlignment;
    std::string referenceAlignment;
    std::string method;
    EXPECT_THROW(
            alignSlidingWindow(
                    query, reference,
                    queryAlignment, referenceAlignment, method,
                    4500, 2, -6, -8, -2, -75, -1),
            std::runtime_error);
    EXPECT_TRUE(method.empty());
}

}  // namespace
