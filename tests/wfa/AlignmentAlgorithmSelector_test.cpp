#include "src/myImportandFunction/AlignmentAlgorithmSelector.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

const anchorwave::AlgorithmCostEstimate *findEstimate(
        const anchorwave::AlignmentSelectionPlan &plan,
        anchorwave::AlignmentCandidate candidate) {
    const auto found = std::find_if(
            plan.estimates.begin(), plan.estimates.end(),
            [candidate](const anchorwave::AlgorithmCostEstimate &estimate) {
                return estimate.candidate == candidate;
            });
    return found == plan.estimates.end() ? nullptr : &*found;
}

anchorwave::AlignmentSelectionPlan select(const std::string &query,
                                           const std::string &reference,
                                           int64_t windowWidth,
                                           double succinctWfaMaxMinutes = 0.0) {
    return anchorwave::makeAlignmentSelectionPlan(
            query, reference, windowWidth,
            -6, -8, -2, -75, -1, succinctWfaMaxMinutes);
}

TEST(AlignmentAlgorithmSelector, UsesSingletrackWfaForIdenticalSequence) {
    const std::string sequence =
            "ACGTTGCAACCGGTTAACCGGTTACGTACGTTGCAACCGGTTA";
    const auto plan = select(sequence, sequence, 100000);

    ASSERT_FALSE(plan.exactCandidates.empty());
    EXPECT_EQ(anchorwave::AlignmentCandidate::SingletrackWfa,
              plan.exactCandidates.front());
    ASSERT_EQ(1ULL, plan.lastResortCandidates.size());
    EXPECT_EQ(anchorwave::AlignmentCandidate::SlidingWindow,
              plan.lastResortCandidates.front());
    EXPECT_TRUE(plan.profile.identical);
    EXPECT_EQ(0ULL, plan.profile.conservativeScore);
    EXPECT_DOUBLE_EQ(1.0, plan.profile.confidence);
    EXPECT_NE(plan.exactCandidates.end(),
              std::find(plan.exactCandidates.begin(),
                        plan.exactCandidates.end(),
                        anchorwave::AlignmentCandidate::MediumWfa));
    EXPECT_NE(plan.exactCandidates.end(),
              std::find(plan.exactCandidates.begin(),
                        plan.exactCandidates.end(),
                        anchorwave::AlignmentCandidate::LowWfa));
}

TEST(AlignmentAlgorithmSelector,
     KeepsStableTieBreakOrderButExposesDynamicRuntimeEstimates) {
    std::mt19937 generator(20260809);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string reference;
    reference.reserve(20000);
    for (int position = 0; position < 20000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 50) {
        query[position] = query[position] == 'A' ? 'C' :
                          query[position] == 'C' ? 'G' :
                          query[position] == 'G' ? 'T' : 'A';
    }

    const auto plan = select(query, reference, 100000);
    const auto *singletrack = findEstimate(
            plan, anchorwave::AlignmentCandidate::SingletrackWfa);
    const auto *high = findEstimate(
            plan, anchorwave::AlignmentCandidate::StandardWfa);
    const auto *medium = findEstimate(
            plan, anchorwave::AlignmentCandidate::MediumWfa);
    const auto *low = findEstimate(
            plan, anchorwave::AlignmentCandidate::LowWfa);
    ASSERT_NE(nullptr, singletrack);
    ASSERT_NE(nullptr, high);
    ASSERT_NE(nullptr, medium);
    ASSERT_NE(nullptr, low);
    ASSERT_TRUE(singletrack->memoryFeasible);
    ASSERT_TRUE(high->memoryFeasible);
    ASSERT_TRUE(medium->memoryFeasible);
    ASSERT_TRUE(low->memoryFeasible);
    EXPECT_LT(singletrack->memoryBytes, high->memoryBytes);
    EXPECT_LT(singletrack->estimatedMinutes, high->estimatedMinutes);
    EXPECT_LT(low->memoryBytes, medium->memoryBytes);
    EXPECT_LT(medium->memoryBytes, high->memoryBytes);
    EXPECT_LT(high->estimatedMinutes, medium->estimatedMinutes);
    EXPECT_LT(medium->estimatedMinutes, low->estimatedMinutes);

    const auto singletrackPosition = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            anchorwave::AlignmentCandidate::SingletrackWfa);
    const auto highPosition = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            anchorwave::AlignmentCandidate::StandardWfa);
    const auto mediumPosition = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            anchorwave::AlignmentCandidate::MediumWfa);
    const auto lowPosition = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            anchorwave::AlignmentCandidate::LowWfa);
    ASSERT_NE(plan.exactCandidates.end(), singletrackPosition);
    ASSERT_NE(plan.exactCandidates.end(), highPosition);
    ASSERT_NE(plan.exactCandidates.end(), mediumPosition);
    ASSERT_NE(plan.exactCandidates.end(), lowPosition);
    EXPECT_LT(singletrackPosition, highPosition);
    EXPECT_LT(highPosition, mediumPosition);
    EXPECT_LT(mediumPosition, lowPosition);
}

TEST(AlignmentAlgorithmSelector,
     SuccinctModesBridgeHighWfaAndFullKswOnLargeExactTasks) {
    std::mt19937 generator(20260810);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string reference;
    reference.reserve(110000);
    for (int position = 0; position < 110000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 80) {
        query[position] = query[position] == 'A' ? 'C' :
                          query[position] == 'C' ? 'G' :
                          query[position] == 'G' ? 'T' : 'A';
    }

    const auto plan = select(query, reference, 100000);
    EXPECT_TRUE(findEstimate(
            plan,
            anchorwave::AlignmentCandidate::SingletrackWfa)->memoryFeasible);
    EXPECT_FALSE(findEstimate(
            plan, anchorwave::AlignmentCandidate::StandardWfa)->memoryFeasible);
    EXPECT_TRUE(findEstimate(
            plan, anchorwave::AlignmentCandidate::MediumWfa)->memoryFeasible);
    EXPECT_TRUE(findEstimate(
            plan, anchorwave::AlignmentCandidate::LowWfa)->memoryFeasible);
    EXPECT_FALSE(findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Full)->memoryFeasible);
    ASSERT_FALSE(plan.exactCandidates.empty());
    EXPECT_EQ(anchorwave::AlignmentCandidate::SingletrackWfa,
              plan.exactCandidates.front());
    const auto *singletrack = findEstimate(
            plan, anchorwave::AlignmentCandidate::SingletrackWfa);
    ASSERT_NE(nullptr, singletrack);
    EXPECT_GT(singletrack->estimatedMinutes, 0.0)
            << "even score-zero alignment has linear scan/CIGAR cost";
    EXPECT_GE(singletrack->estimatedMinutesP90,
              singletrack->estimatedMinutes);
}

TEST(AlignmentAlgorithmSelector, ModelsFullKsw2RuntimeInsideExactTier) {
    const std::string query(5000, 'A');
    const std::string reference(5000, 'T');
    const auto plan = select(query, reference, 100000);
    const auto *full = findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Full);
    ASSERT_NE(nullptr, full);
    ASSERT_TRUE(full->memoryFeasible);
    EXPECT_GT(full->estimatedMinutes, 0.0);
    EXPECT_GE(full->estimatedMinutesP90, full->estimatedMinutes);
    EXPECT_GT(full->successProbability, 0.99);
}

TEST(AlignmentAlgorithmSelector,
     ProcessMemoryCannotRaiseFullKsw2AboveWindowSquare) {
    const std::string query(20000, 'A');
    const std::string reference(20000, 'T');
    const auto legacy = select(query, reference, 10001);
    const auto *legacyFull = findEstimate(
            legacy, anchorwave::AlignmentCandidate::Ksw2Full);
    const auto *legacyBand = findEstimate(
            legacy, anchorwave::AlignmentCandidate::Ksw2Banded);
    ASSERT_NE(nullptr, legacyFull);
    ASSERT_NE(nullptr, legacyBand);
    ASSERT_FALSE(legacyFull->memoryFeasible);
    ASSERT_TRUE(legacyBand->memoryFeasible);
    EXPECT_LE(legacyBand->memoryBytes, 10001ULL * 10001ULL);
    EXPECT_LT(legacy.ksw2BandWidth,
              static_cast<int64_t>(std::max(query.size(), reference.size())));

    const auto processRich = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 10001,
            -6, -8, -2, -75, -1, 0.0, 0,
            legacyFull->memoryBytes);
    const auto *processRichFull = findEstimate(
            processRich, anchorwave::AlignmentCandidate::Ksw2Full);
    const auto *processRichBand = findEstimate(
            processRich, anchorwave::AlignmentCandidate::Ksw2Banded);
    ASSERT_NE(nullptr, processRichFull);
    ASSERT_NE(nullptr, processRichBand);
    EXPECT_FALSE(processRichFull->memoryFeasible);
    EXPECT_FALSE(processRichFull->feasible);
    EXPECT_EQ(processRich.exactCandidates.end(),
              std::find(processRich.exactCandidates.begin(),
                        processRich.exactCandidates.end(),
                        anchorwave::AlignmentCandidate::Ksw2Full));
    EXPECT_TRUE(processRichBand->memoryFeasible);
    EXPECT_EQ(legacyBand->memoryBytes, processRichBand->memoryBytes);
    EXPECT_EQ(legacy.ksw2BandWidth, processRich.ksw2BandWidth);
}

TEST(AlignmentAlgorithmSelector,
     ScoreCertifiedKsw2CannotBorrowAboveWindowSquare) {
    const std::string query(50000, 'A');
    const std::string reference(50000, 'T');
    const auto constrained = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 1000,
            -6, -8, -2, -75, -1, 0.0, 0, 0, false);
    const auto *constrainedCertified = findEstimate(
            constrained,
            anchorwave::AlignmentCandidate::Ksw2ScoreCertified);
    ASSERT_NE(nullptr, constrainedCertified);
    EXPECT_FALSE(constrainedCertified->memoryFeasible);

    constexpr uint64_t largeProcessHint = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    const auto processRich = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 1000,
            -6, -8, -2, -75, -1, 0.0, 0, largeProcessHint, true);
    const auto *processRichCertified = findEstimate(
            processRich,
            anchorwave::AlignmentCandidate::Ksw2ScoreCertified);
    ASSERT_NE(nullptr, processRichCertified);
    EXPECT_FALSE(processRichCertified->memoryFeasible);
    EXPECT_EQ(constrained.ksw2CertifiedMaximumBandWidth,
              processRich.ksw2CertifiedMaximumBandWidth);

    constexpr int64_t fittingWindow = 30000;
    const uint64_t fittingBudget = static_cast<uint64_t>(fittingWindow) *
                                   static_cast<uint64_t>(fittingWindow);
    const auto fitting = anchorwave::makeAlignmentSelectionPlan(
            query, reference, fittingWindow,
            -6, -8, -2, -75, -1, 0.0, 0, 0, true);
    const auto *fittingCertified = findEstimate(
            fitting, anchorwave::AlignmentCandidate::Ksw2ScoreCertified);
    const auto *fittingFull = findEstimate(
            fitting, anchorwave::AlignmentCandidate::Ksw2Full);
    ASSERT_NE(nullptr, fittingCertified);
    ASSERT_NE(nullptr, fittingFull);
    EXPECT_TRUE(fittingCertified->memoryFeasible);
    EXPECT_FALSE(fittingFull->memoryFeasible);
    EXPECT_LE(fittingCertified->memoryBytes, fittingBudget);

    const auto certifiedReservation = [&query, &reference](
            int64_t initialBand, int64_t maximumBand, bool cumulative) {
        const uint64_t totalLength = query.size() + reference.size();
        const uint64_t scoreOnly = 32ULL * 1024ULL * 1024ULL +
                                   96ULL * totalLength;
        uint64_t bytes = scoreOnly;
        uint64_t band = static_cast<uint64_t>(initialBand);
        const uint64_t maximum = static_cast<uint64_t>(maximumBand);
        while (true) {
            const uint64_t traceback =
                    anchorwave::ksw2TracebackResidentMemoryBytes(
                    query.size(), reference.size(), band);
            bytes = cumulative ? bytes + traceback
                               : std::max(bytes, traceback);
            if (band == maximum) {
                return bytes;
            }
            band = std::min<uint64_t>(maximum, band * 2);
        }
    };
    EXPECT_EQ(certifiedReservation(
                      fitting.ksw2CertifiedInitialBandWidth,
                      fitting.ksw2CertifiedMaximumBandWidth, false),
              fittingCertified->memoryBytes);
    ASSERT_GT(fitting.ksw2CertifiedMaximumBandWidth,
              fitting.ksw2CertifiedInitialBandWidth);
    EXPECT_LT(fittingCertified->memoryBytes,
              certifiedReservation(
                      fitting.ksw2CertifiedInitialBandWidth,
                      fitting.ksw2CertifiedMaximumBandWidth, true))
            << "sequential score/trace attempts reserve peak, not cumulative RSS";
}

TEST(AlignmentAlgorithmSelector,
     ScoreCertifiedInitialBandCoversUnavoidableLengthDifference) {
    const std::string query(50000, 'A');
    const std::string reference(48000, 'T');
    const auto plan = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 30000,
            -6, -8, -2, -75, -1, 0.0,
            0, 0, true);
    ASSERT_TRUE(findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2ScoreCertified)
                        ->memoryFeasible);
    EXPECT_GE(plan.ksw2CertifiedInitialBandWidth,
              static_cast<int64_t>(query.size() - reference.size()));
}

TEST(AlignmentAlgorithmSelector,
     HighWfaFastLaneRequiresImmediateFastestCompletion) {
    anchorwave::AlgorithmCostEstimate high;
    high.candidate = anchorwave::AlignmentCandidate::SingletrackWfa;
    high.estimatedMinutes = 1.0;
    high.estimatedMinutesP90 = 1.2;
    high.successProbability = 0.99;
    const double highRisk = anchorwave::alignmentRiskAdjustedMinutes(high);

    EXPECT_TRUE(anchorwave::highWfaFastLaneEligible(
            high.candidate, high, 0.0, highRisk));
    EXPECT_FALSE(anchorwave::highWfaFastLaneEligible(
            high.candidate, high, 0.01, highRisk));
    EXPECT_FALSE(anchorwave::highWfaFastLaneEligible(
            high.candidate, high, 0.0, highRisk / 1.05));
    EXPECT_FALSE(anchorwave::highWfaFastLaneEligible(
            anchorwave::AlignmentCandidate::Ksw2Full,
            high, 0.0, highRisk));
}

TEST(AlignmentAlgorithmSelector,
     ParksForShortFastExactWaitInsteadOfSlowLowWfa) {
    anchorwave::AlgorithmCostEstimate ksw;
    ksw.candidate = anchorwave::AlignmentCandidate::Ksw2Full;
    ksw.estimatedMinutes = 0.40;
    ksw.estimatedMinutesP90 = 0.50;
    anchorwave::AlgorithmCostEstimate lowWfa;
    lowWfa.candidate = anchorwave::AlignmentCandidate::LowWfa;
    lowWfa.estimatedMinutes = 6.0;
    lowWfa.estimatedMinutesP90 = 18.0;

    EXPECT_TRUE(anchorwave::fastExactDominatesSlowExact(
            ksw.candidate, ksw, 0.20,
            lowWfa.candidate, lowWfa));
    EXPECT_TRUE(anchorwave::fastExactDominatesSlowExact(
            anchorwave::AlignmentCandidate::SingletrackWfa, ksw, 0.20,
            lowWfa.candidate, lowWfa));
}

TEST(AlignmentAlgorithmSelector,
     LongFastExactWaitDoesNotBlockRunnableLowWfa) {
    anchorwave::AlgorithmCostEstimate ksw;
    ksw.candidate = anchorwave::AlignmentCandidate::Ksw2Full;
    ksw.estimatedMinutes = 0.500550589;
    ksw.estimatedMinutesP90 = 0.609016108;
    anchorwave::AlgorithmCostEstimate lowWfa;
    lowWfa.candidate = anchorwave::AlignmentCandidate::LowWfa;
    lowWfa.estimatedMinutes = 2.09085546;
    lowWfa.estimatedMinutesP90 = 6.14282379;

    EXPECT_FALSE(anchorwave::fastExactDominatesSlowExact(
            ksw.candidate, ksw, 19.0999198,
            lowWfa.candidate, lowWfa));
}

TEST(AlignmentAlgorithmSelector,
     FastExactDominanceRejectsWrongTiersAndUnboundedWaits) {
    anchorwave::AlgorithmCostEstimate fast;
    fast.estimatedMinutes = 0.1;
    fast.estimatedMinutesP90 = 0.2;
    anchorwave::AlgorithmCostEstimate slow;
    slow.estimatedMinutes = 5.0;
    slow.estimatedMinutesP90 = 10.0;

    EXPECT_FALSE(anchorwave::fastExactDominatesSlowExact(
            anchorwave::AlignmentCandidate::MediumWfa, fast, 0.1,
            anchorwave::AlignmentCandidate::Ksw2Banded, slow));
    EXPECT_FALSE(anchorwave::fastExactDominatesSlowExact(
            anchorwave::AlignmentCandidate::Ksw2Full, fast, 0.1,
            anchorwave::AlignmentCandidate::Ksw2Banded, slow));
    EXPECT_FALSE(anchorwave::fastExactDominatesSlowExact(
            anchorwave::AlignmentCandidate::Ksw2Full, fast,
            std::numeric_limits<double>::infinity(),
            anchorwave::AlignmentCandidate::LowWfa, slow));
    EXPECT_FALSE(anchorwave::fastExactDominatesSlowExact(
            anchorwave::AlignmentCandidate::Ksw2Full, fast, 2.1,
            anchorwave::AlignmentCandidate::LowWfa, slow));
}

TEST(AlignmentAlgorithmSelector,
     ProductionSelectorDoesNotExposeBiWfa) {
    const std::string query(5000, 'A');
    const std::string reference(5000, 'T');
    const auto plan = select(query, reference, 100000);
    const auto *singletrack = findEstimate(
            plan, anchorwave::AlignmentCandidate::SingletrackWfa);
    ASSERT_NE(nullptr, singletrack);
    EXPECT_GE(singletrack->estimatedMinutesP90,
              singletrack->estimatedMinutes);
}

TEST(AlignmentAlgorithmSelector,
     PricesGenericKswScoringForAmbiguousGenomeBases) {
    std::string reference;
    reference.reserve(5000);
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int position = 0; position < 5000; ++position) {
        reference.push_back(bases[position % 4]);
    }
    std::string cleanQuery = reference;
    cleanQuery[2500] = cleanQuery[2500] == 'A' ? 'C' : 'A';
    std::string ambiguousQuery = cleanQuery;
    ambiguousQuery[2000] = 'N';

    const auto clean = select(cleanQuery, reference, 100000);
    const auto ambiguous = select(ambiguousQuery, reference, 100000);
    const auto *cleanKsw = findEstimate(
            clean, anchorwave::AlignmentCandidate::Ksw2Full);
    const auto *ambiguousKsw = findEstimate(
            ambiguous, anchorwave::AlignmentCandidate::Ksw2Full);
    ASSERT_NE(nullptr, cleanKsw);
    ASSERT_NE(nullptr, ambiguousKsw);
    EXPECT_DOUBLE_EQ(0.0, clean.profile.ambiguousBaseFraction);
    EXPECT_GT(ambiguous.profile.ambiguousBaseFraction, 0.0);
    EXPECT_GT(ambiguousKsw->workUnits, cleanKsw->workUnits * 1.8L);
    EXPECT_GT(ambiguousKsw->estimatedMinutes,
              cleanKsw->estimatedMinutes * 1.7);
}

TEST(AlignmentAlgorithmSelector,
     IdenticalNContainingSequenceStillPricesGenericKswScoring) {
    std::string sequence;
    sequence.reserve(5000);
    for (int position = 0; position < 5000; ++position) {
        sequence.push_back(position % 101 == 0 ? 'N' :
                           "ACGT"[position % 4]);
    }
    const auto plan = select(sequence, sequence, 100000);
    const auto *full = findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Full);
    ASSERT_NE(nullptr, full);
    EXPECT_TRUE(plan.profile.identical);
    EXPECT_GT(plan.profile.ambiguousBaseFraction, 0.0);
    EXPECT_GT(full->workUnits,
              static_cast<long double>(sequence.size()) * sequence.size() *
                      1.8L);
}

TEST(AlignmentAlgorithmSelector, DirectlyBypassesAllLargeMatrixAlgorithms) {
    const std::string query(20000, 'A');
    const std::string reference(20000, 'T');
    const auto plan = select(query, reference, 100);

    EXPECT_TRUE(plan.exactCandidates.empty());
    EXPECT_TRUE(plan.bandedCandidates.empty());
    ASSERT_EQ(1ULL, plan.lastResortCandidates.size());
    EXPECT_EQ(anchorwave::AlignmentCandidate::SlidingWindow,
              plan.lastResortCandidates.front());
    ASSERT_NE(nullptr, findEstimate(
            plan, anchorwave::AlignmentCandidate::StandardWfa));
    EXPECT_FALSE(findEstimate(
            plan, anchorwave::AlignmentCandidate::StandardWfa)->feasible);
    EXPECT_FALSE(findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Full)->feasible);
    EXPECT_FALSE(findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Banded)->feasible);
}

TEST(AlignmentAlgorithmSelector,
     FullKswSuppressesScoreCertifiedRescue) {
    const std::string query(5000, 'A');
    const std::string reference(5000, 'T');
    const auto plan = select(query, reference, 100000);

    ASSERT_FALSE(plan.exactCandidates.empty());
    const auto certifiedKsw2 = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            anchorwave::AlignmentCandidate::Ksw2ScoreCertified);
    const auto fullKsw2 = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            anchorwave::AlignmentCandidate::Ksw2Full);
    EXPECT_EQ(plan.exactCandidates.end(), certifiedKsw2);
    ASSERT_NE(plan.exactCandidates.end(), fullKsw2);
    ASSERT_EQ(1ULL, plan.lastResortCandidates.size());
    EXPECT_EQ(anchorwave::AlignmentCandidate::SlidingWindow,
              plan.lastResortCandidates.front());
    EXPECT_TRUE(findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Full)->feasible);
}

TEST(AlignmentAlgorithmSelector,
     BandedKswUsesTheWidestReachableBandInsideWindowSquare) {
    const std::string query(20000, 'A');
    const std::string reference(20000, 'T');
    const auto constrained = select(query, reference, 10001);

    const auto *constrainedBand = findEstimate(
            constrained, anchorwave::AlignmentCandidate::Ksw2Banded);
    ASSERT_NE(nullptr, constrainedBand);
    EXPECT_TRUE(constrainedBand->feasible);
    EXPECT_GT(constrained.ksw2BandWidth, 0);
    EXPECT_LT(constrained.ksw2BandWidth,
              static_cast<int64_t>(std::max(query.size(), reference.size())));
    EXPECT_LE(constrainedBand->memoryBytes, 10001ULL * 10001ULL);

    const auto fitting = select(query, reference, 30000);
    EXPECT_TRUE(findEstimate(
            fitting, anchorwave::AlignmentCandidate::Ksw2Banded)->feasible);
    EXPECT_GE(fitting.ksw2BandWidth,
              static_cast<int64_t>(std::max(query.size(), reference.size())));
    EXPECT_NE(fitting.bandedCandidates.end(),
              std::find(fitting.bandedCandidates.begin(),
                        fitting.bandedCandidates.end(),
                        anchorwave::AlignmentCandidate::Ksw2Banded));
}

TEST(AlignmentAlgorithmSelector,
     LowSimilarityBandHasNonzeroUnresolvedPathLoss) {
    const std::string query(20000, 'A');
    const std::string reference(20000, 'T');
    const auto plan = select(query, reference, 10001);
    const auto *banded = findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Banded);
    const auto *sliding = findEstimate(
            plan, anchorwave::AlignmentCandidate::SlidingWindow);
    ASSERT_NE(nullptr, banded);
    ASSERT_NE(nullptr, sliding);
    ASSERT_TRUE(banded->feasible);
    ASSERT_EQ(2ULL, plan.approximateCandidates.size());
    EXPECT_LT(plan.profile.sketchJaccard, 0.01);
    EXPECT_GT(banded->predictedScoreLoss, 0.0);
    EXPECT_GT(banded->predictedScoreLossP90,
              banded->predictedScoreLoss);
}

TEST(AlignmentAlgorithmSelector,
     ApproximateFallbackPrefersABandThatCoversThePredictedPath) {
    std::mt19937 generator(20260810);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string reference;
    reference.reserve(200000);
    for (int position = 0; position < 200000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 1000) {
        query[position] = query[position] == 'A' ? 'C' : 'A';
    }

    const auto plan = select(query, reference, 10000);
    const auto *banded = findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Banded);
    const auto *sliding = findEstimate(
            plan, anchorwave::AlignmentCandidate::SlidingWindow);
    ASSERT_NE(nullptr, banded);
    ASSERT_NE(nullptr, sliding);
    ASSERT_TRUE(banded->feasible);
    ASSERT_EQ(2ULL, plan.approximateCandidates.size());
    EXPECT_EQ(anchorwave::AlignmentCandidate::Ksw2Banded,
              plan.approximateCandidates.front());
    EXPECT_LE(banded->predictedScoreLoss,
              sliding->predictedScoreLoss);
}

TEST(AlignmentAlgorithmSelector,
     ApproximateFallbackPrefersSlidingWhenTheSafeBandMissesAStructuralJump) {
    std::mt19937 generator(20260811);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string reference;
    reference.reserve(200000);
    for (int position = 0; position < 200000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference.substr(0, 80000);
    for (int position = 0; position < 40000; ++position) {
        query.push_back(bases[base(generator)]);
    }
    // A 40 kb insertion followed by a compensating 40 kb terminal deletion
    // preserves total length while shifting the second collinear segment far
    // outside the memory-safe narrow band.
    query.append(reference, 80000, 80000);

    const auto plan = select(query, reference, 10000);
    const auto *banded = findEstimate(
            plan, anchorwave::AlignmentCandidate::Ksw2Banded);
    const auto *sliding = findEstimate(
            plan, anchorwave::AlignmentCandidate::SlidingWindow);
    ASSERT_NE(nullptr, banded);
    ASSERT_NE(nullptr, sliding);
    ASSERT_TRUE(banded->feasible);
    ASSERT_EQ(2ULL, plan.approximateCandidates.size());
    EXPECT_GT(plan.profile.maximumDiagonalJump,
              static_cast<double>(plan.ksw2BandWidth));
    EXPECT_EQ(anchorwave::AlignmentCandidate::SlidingWindow,
              plan.approximateCandidates.front());
    EXPECT_LT(sliding->predictedScoreLossP90,
              banded->predictedScoreLossP90);
}

TEST(AlignmentAlgorithmSelector,
     SlidingWindowIsReducedUntilItsKswTracebackFitsWindowSquare) {
    const std::string query(110000, 'A');
    const std::string reference(110000, 'T');
    constexpr int64_t requestedWindow = 100000;
    const auto plan = select(query, reference, requestedWindow);
    const uint64_t budget = static_cast<uint64_t>(requestedWindow) *
                            static_cast<uint64_t>(requestedWindow);

    ASSERT_GT(plan.slidingWindowWidth, 0);
    EXPECT_LE(plan.slidingWindowWidth, requestedWindow);
    EXPECT_LE(anchorwave::ksw2TracebackResidentMemoryBytes(
                      static_cast<uint64_t>(plan.slidingWindowWidth),
                      static_cast<uint64_t>(plan.slidingWindowWidth),
                      static_cast<uint64_t>(plan.slidingWindowWidth)),
              budget);
    if (plan.slidingWindowWidth < requestedWindow) {
        const uint64_t next = static_cast<uint64_t>(
                plan.slidingWindowWidth + 1);
        EXPECT_GT(anchorwave::ksw2TracebackResidentMemoryBytes(
                          next, next, next),
                  budget);
    }
}

TEST(AlignmentAlgorithmSelector,
     Ksw2ResidentModelCountsTouchedPagesInsteadOfVirtualRectangle) {
    constexpr uint64_t queryLength = 10000;
    constexpr uint64_t referenceLength = 12000;
    const uint64_t full = anchorwave::ksw2TracebackResidentMemoryBytes(
            queryLength, referenceLength,
            std::max(queryLength, referenceLength));
    const uint64_t narrow = anchorwave::ksw2TracebackResidentMemoryBytes(
            queryLength, referenceLength, 128);
    // The virtual row rectangle is an upper bound on the traceback region
    // alone. A full resident prediction must be materially below it because
    // rising/falling anti-diagonals touch only their live SIMD prefixes.
    const uint64_t conservativeVirtualRectangle =
            (queryLength + referenceLength - 1) *
            (queryLength + 64);
    EXPECT_GT(full, narrow);
    EXPECT_LT(full, conservativeVirtualRectangle);
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::ksw2TracebackResidentMemoryBytes(
                      static_cast<uint64_t>(INT_MAX) + 1, 10, 10));
}

TEST(AlignmentAlgorithmSelector,
     RemovedBiWfaNeverAppearsInProductionPlanOrTelemetry) {
    std::mt19937 generator(20260808);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string reference;
    reference.reserve(110000);
    for (int position = 0; position < 110000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 20) {
        query[position] = query[position] == 'A' ? 'C' :
                          query[position] == 'C' ? 'G' :
                          query[position] == 'G' ? 'T' : 'A';
    }

    // BiWFA remains independently benchmarkable through WfaAlignment, but is
    // absent from the selector's type system, plan and telemetry.
    const auto plan = select(query, reference, 55001);
    for (const auto &estimate : plan.estimates) {
        EXPECT_STRNE("BIWFA",
                     anchorwave::alignmentCandidateName(estimate.candidate));
    }
    for (const auto candidate : plan.exactCandidates) {
        EXPECT_STRNE("BIWFA", anchorwave::alignmentCandidateName(candidate));
    }
    for (const auto candidate : plan.approximateCandidates) {
        EXPECT_STRNE("BIWFA", anchorwave::alignmentCandidateName(candidate));
    }

    anchorwave::resetAlignmentSelectionTelemetry();
    anchorwave::recordAlignmentSelectionPlan(plan);
    const auto telemetry =
            anchorwave::alignmentSelectionTelemetrySnapshot();
    EXPECT_EQ(1ULL, telemetry.exactTierIntervals);
    EXPECT_EQ(0ULL, telemetry.bandedOnlyIntervals);
    EXPECT_EQ(0ULL, telemetry.slidingOnlyIntervals);
}

TEST(AlignmentAlgorithmSelector,
     ExactFirstKeepsEveryMemoryFeasibleExactCandidate) {
    std::string reference;
    reference.reserve(110000);
    std::mt19937 generator(20260811);
    std::uniform_int_distribution<int> base(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int position = 0; position < 110000; ++position) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 80) {
        query[position] = query[position] == 'A' ? 'C' :
                          query[position] == 'C' ? 'G' :
                          query[position] == 'G' ? 'T' : 'A';
    }

    const auto unlimited = select(query, reference, 100000, 0.0);
    const auto *unlimitedMediumWfa = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::MediumWfa);
    const auto *unlimitedLowWfa = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::LowWfa);
    const auto *unlimitedSingletrack = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::SingletrackWfa);
    const auto *unlimitedHigh = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::StandardWfa);
    ASSERT_NE(nullptr, unlimitedMediumWfa);
    ASSERT_NE(nullptr, unlimitedLowWfa);
    ASSERT_NE(nullptr, unlimitedSingletrack);
    ASSERT_NE(nullptr, unlimitedHigh);
    ASSERT_TRUE(unlimitedSingletrack->memoryFeasible ||
                unlimitedHigh->memoryFeasible);
    ASSERT_TRUE(unlimitedMediumWfa->memoryFeasible);
    ASSERT_TRUE(unlimitedLowWfa->memoryFeasible);
    EXPECT_GT(unlimitedMediumWfa->estimatedMinutes, 0.0);
    EXPECT_GT(unlimitedLowWfa->estimatedMinutes, 0.0);
    EXPECT_TRUE(unlimitedMediumWfa->timeFeasible);
    EXPECT_TRUE(unlimitedLowWfa->timeFeasible);
    EXPECT_TRUE(unlimitedMediumWfa->feasible);
    EXPECT_TRUE(unlimitedLowWfa->feasible);

    const double cutoff = 0.5 * std::min({
            unlimitedMediumWfa->estimatedMinutes,
            unlimitedLowWfa->estimatedMinutes});
    const auto limited = select(query, reference, 100000, cutoff);
    const auto *limitedMediumWfa = findEstimate(
            limited, anchorwave::AlignmentCandidate::MediumWfa);
    const auto *limitedLowWfa = findEstimate(
            limited, anchorwave::AlignmentCandidate::LowWfa);
    ASSERT_NE(nullptr, limitedMediumWfa);
    ASSERT_NE(nullptr, limitedLowWfa);
    EXPECT_TRUE(limitedMediumWfa->memoryFeasible);
    EXPECT_TRUE(limitedLowWfa->memoryFeasible);
    EXPECT_GT(limitedMediumWfa->estimatedMinutes, cutoff);
    EXPECT_GT(limitedLowWfa->estimatedMinutes, cutoff);
    EXPECT_FALSE(limitedMediumWfa->timeFeasible);
    EXPECT_FALSE(limitedLowWfa->timeFeasible);
    EXPECT_FALSE(limitedMediumWfa->withinConfiguredTimeLimit);
    EXPECT_FALSE(limitedLowWfa->withinConfiguredTimeLimit);
    EXPECT_FALSE(limitedMediumWfa->feasible);
    EXPECT_FALSE(limitedLowWfa->feasible);
    for (const auto candidate : {
            anchorwave::AlignmentCandidate::MediumWfa,
            anchorwave::AlignmentCandidate::LowWfa}) {
        EXPECT_EQ(limited.exactCandidates.end(),
                  std::find(limited.exactCandidates.begin(),
                            limited.exactCandidates.end(), candidate));
    }
    EXPECT_FALSE(limited.lastResortCandidates.empty());

    const auto defaultLimited = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 100000,
            -6, -8, -2, -75, -1);
    for (const auto candidate : {
            anchorwave::AlignmentCandidate::MediumWfa,
            anchorwave::AlignmentCandidate::LowWfa}) {
        const auto *cost = findEstimate(defaultLimited, candidate);
        ASSERT_NE(nullptr, cost);
        EXPECT_TRUE(cost->timeFeasible);
        EXPECT_EQ(cost->memoryFeasible && cost->timeFeasible,
                  cost->feasible);
    }

    anchorwave::resetAlignmentSelectionTelemetry();
    anchorwave::recordAlignmentSelectionPlan(limited);
    const auto telemetry =
            anchorwave::alignmentSelectionTelemetrySnapshot();
    EXPECT_EQ(1ULL, telemetry.mediumWfaTimeRejects);
    EXPECT_EQ(1ULL, telemetry.lowWfaTimeRejects);
}

TEST(AlignmentAlgorithmSelector,
     ExactTimeLimitAppliesAfterBothHighWfasExceedMemory) {
    const std::string query(5000, 'A');
    const std::string reference(5000, 'T');
    const auto unlimited = select(query, reference, 100000, 0.0);
    const auto *singletrack = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::SingletrackWfa);
    const auto *high = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::StandardWfa);
    const auto *medium = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::MediumWfa);
    const auto *low = findEstimate(
            unlimited, anchorwave::AlignmentCandidate::LowWfa);
    ASSERT_NE(nullptr, singletrack);
    ASSERT_NE(nullptr, high);
    ASSERT_NE(nullptr, medium);
    ASSERT_NE(nullptr, low);
    ASSERT_FALSE(singletrack->memoryFeasible);
    ASSERT_FALSE(high->memoryFeasible);

    const double cutoff = 0.5 * std::min(
            medium->estimatedMinutes, low->estimatedMinutes);
    const auto limited = select(query, reference, 100000, cutoff);
    for (const auto candidate : {
            anchorwave::AlignmentCandidate::MediumWfa,
            anchorwave::AlignmentCandidate::LowWfa}) {
        const auto *cost = findEstimate(limited, candidate);
        ASSERT_NE(nullptr, cost);
        EXPECT_GT(cost->estimatedMinutes, cutoff);
        EXPECT_FALSE(cost->withinConfiguredTimeLimit);
        EXPECT_FALSE(cost->timeFeasible);
        EXPECT_EQ(limited.exactCandidates.end(),
                  std::find(limited.exactCandidates.begin(),
                            limited.exactCandidates.end(), candidate));
    }
    anchorwave::resetAlignmentSelectionTelemetry();
    anchorwave::recordAlignmentSelectionPlan(limited);
    const auto telemetry =
            anchorwave::alignmentSelectionTelemetrySnapshot();
    EXPECT_EQ(0ULL, telemetry.mediumWfaTimeRejects);
    EXPECT_EQ(0ULL, telemetry.lowWfaTimeRejects);
}

TEST(AlignmentAlgorithmSelector,
     BalancedModeAppliesTimeGatesWhenHighWfaFitsWindowSquare) {
    const std::string query(5000, 'A');
    const std::string reference(5000, 'T');
    constexpr int64_t fittingWindow = 1000000;
    const auto legacy = select(query, reference, fittingWindow, 0.0);
    const auto *legacySingletrack = findEstimate(
            legacy, anchorwave::AlignmentCandidate::SingletrackWfa);
    const auto *legacyHigh = findEstimate(
            legacy, anchorwave::AlignmentCandidate::StandardWfa);
    ASSERT_NE(nullptr, legacySingletrack);
    ASSERT_NE(nullptr, legacyHigh);
    ASSERT_TRUE(legacySingletrack->memoryFeasible);
    ASSERT_FALSE(legacyHigh->memoryFeasible)
            << "standard high retains its smaller historical trial ceiling";

    const auto *legacyMedium = findEstimate(
            legacy, anchorwave::AlignmentCandidate::MediumWfa);
    ASSERT_NE(nullptr, legacyMedium);
    const double cutoff = 0.5 * legacyMedium->estimatedMinutes;
    const auto elastic = select(query, reference, fittingWindow, cutoff);
    EXPECT_TRUE(findEstimate(
            elastic,
            anchorwave::AlignmentCandidate::SingletrackWfa)->memoryFeasible);
    EXPECT_FALSE(findEstimate(
            elastic,
            anchorwave::AlignmentCandidate::StandardWfa)->memoryFeasible);
    for (const auto candidate : {
            anchorwave::AlignmentCandidate::MediumWfa,
            anchorwave::AlignmentCandidate::LowWfa}) {
        const auto *estimate = findEstimate(elastic, candidate);
        ASSERT_NE(nullptr, estimate);
        EXPECT_GT(estimate->estimatedMinutes, cutoff);
        EXPECT_FALSE(estimate->withinConfiguredTimeLimit);
        EXPECT_FALSE(estimate->timeFeasible);
        EXPECT_EQ(elastic.exactCandidates.end(),
                  std::find(elastic.exactCandidates.begin(),
                            elastic.exactCandidates.end(), candidate));
    }
}

TEST(AlignmentAlgorithmSelector,
     DefaultBalancedModeUsesStrictCalibratedP50AndP90) {
    const std::string query(70000, 'A');
    const std::string reference(70000, 'T');
    const auto defaultPlan = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 100000,
            -6, -8, -2, -75, -1);
    const auto *defaultFull = findEstimate(
            defaultPlan, anchorwave::AlignmentCandidate::Ksw2Full);
    ASSERT_NE(nullptr, defaultFull);
    ASSERT_TRUE(defaultFull->memoryFeasible);
    ASSERT_GT(defaultFull->estimatedMinutes, 0.0);
    ASSERT_GT(defaultFull->estimatedMinutesP90,
              defaultFull->estimatedMinutes);
    EXPECT_EQ(anchorwave::exactCandidateWithinTimeLimit(
                      defaultFull->estimatedMinutes,
                      defaultFull->estimatedMinutesP90, 30.0),
              defaultFull->timeFeasible);
    const double strictCutoff = std::min(
            defaultFull->estimatedMinutes * 0.5,
            defaultFull->estimatedMinutesP90 / 3.0);
    const auto strictPlan = select(
            query, reference, 100000, strictCutoff);
    const auto *strictFull = findEstimate(
            strictPlan, anchorwave::AlignmentCandidate::Ksw2Full);
    ASSERT_NE(nullptr, strictFull);
    EXPECT_FALSE(strictFull->timeFeasible);
}

TEST(AlignmentAlgorithmSelector,
     ExactTimePolicyUsesZeroForExactFirstAndPositiveForBalanced) {
    EXPECT_TRUE(anchorwave::exactCandidateWithinTimeLimit(
            1000000.0, 2000000.0, 0.0));
    EXPECT_FALSE(anchorwave::exactCandidateWithinTimeLimit(
            29.0, 80.0, 30.0));
    EXPECT_FALSE(anchorwave::exactCandidateWithinTimeLimit(
            31.0, 44.0, 30.0));
    EXPECT_TRUE(anchorwave::exactCandidateWithinTimeLimit(
            29.0, 30.0, 30.0));
    EXPECT_FALSE(anchorwave::exactCandidateWithinTimeLimit(
            31.0, 46.0, 30.0));
    EXPECT_THROW(anchorwave::exactCandidateWithinTimeLimit(
                         1.0, 1.0, -1.0), std::invalid_argument);

    anchorwave::configureExactAlignmentMaximumEstimatedMinutes(0.0);
    EXPECT_DOUBLE_EQ(0.0,
                     anchorwave::configuredExactAlignmentMaximumEstimatedMinutes());
    EXPECT_THROW(anchorwave::configureExactAlignmentMaximumEstimatedMinutes(
                         -1.0),
                 std::invalid_argument);
    EXPECT_THROW(anchorwave::configureExactAlignmentMaximumEstimatedMinutes(
                         std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    anchorwave::configureExactAlignmentMaximumEstimatedMinutes(30.0);
    EXPECT_DOUBLE_EQ(30.0,
                     anchorwave::configuredExactAlignmentMaximumEstimatedMinutes());
}

TEST(AlignmentAlgorithmSelector,
     LongTaskPriorityUsesTheFastestRiskAdjustedCandidateInTheBestTier) {
    const std::string query(20000, 'A');
    const std::string reference(20000, 'T');
    const auto plan = select(query, reference, 100000, 0.0);
    ASSERT_FALSE(plan.exactCandidates.empty());

    double expected = std::numeric_limits<double>::infinity();
    for (const auto candidate : plan.exactCandidates) {
        const auto *estimate = findEstimate(plan, candidate);
        ASSERT_NE(nullptr, estimate);
        expected = std::min(
                expected,
                anchorwave::alignmentRiskAdjustedMinutes(*estimate));
    }
    EXPECT_DOUBLE_EQ(expected,
                     anchorwave::alignmentSelectionPriorityMinutes(plan));
}

TEST(AlignmentAlgorithmSelector,
     ReusablePlanProvenanceRejectsContentPolicyAndResourceChanges) {
    const std::string query = "ACGTTGCAACCGGTTA";
    const std::string reference = "ACGTCGCAACCGGTTA";
    const auto plan = anchorwave::makeAlignmentSelectionPlan(
            query, reference, 100000,
            -6, -8, -2, -75, -1, 30.0, 1234, 5678);
    EXPECT_TRUE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 100000,
            -6, -8, -2, -75, -1, 30.0, 1234, 5678));

    std::string changedQuery = query;
    changedQuery[0] = changedQuery[0] == 'A' ? 'C' : 'A';
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, changedQuery, reference, 100000,
            -6, -8, -2, -75, -1, 30.0, 1234, 5678));
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 99999,
            -6, -8, -2, -75, -1, 30.0, 1234, 5678));
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 100000,
            -7, -8, -2, -75, -1, 30.0, 1234, 5678));
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 100000,
            -6, -8, -2, -75, -1, 0.0, 1234, 5678));
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 100000,
            -6, -8, -2, -75, -1, 30.0, 1235, 5678));
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 100000,
            -6, -8, -2, -75, -1, 30.0, 1234, 5679));
    EXPECT_FALSE(anchorwave::alignmentSelectionPlanMatches(
            plan, query, reference, 100000,
            -6, -8, -2, -75, -1, 30.0, 1234, 5678, true));
}

TEST(AlignmentAlgorithmSelector,
     BalancedModeAppliesTheSameTimeGateToEveryExactEngine) {
    std::string reference;
    reference.reserve(20000);
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int position = 0; position < 20000; ++position) {
        reference.push_back(bases[position % 4]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 31) {
        query[position] = query[position] == 'A' ? 'T' : 'A';
    }

    const auto exactFirst = select(query, reference, 100000, 0.0);
    double strictLimit = std::numeric_limits<double>::infinity();
    for (const auto candidate : {
            anchorwave::AlignmentCandidate::SingletrackWfa,
            anchorwave::AlignmentCandidate::StandardWfa,
            anchorwave::AlignmentCandidate::MediumWfa,
            anchorwave::AlignmentCandidate::LowWfa,
            anchorwave::AlignmentCandidate::Ksw2ScoreCertified,
            anchorwave::AlignmentCandidate::Ksw2Full}) {
        const auto *estimate = findEstimate(exactFirst, candidate);
        ASSERT_NE(nullptr, estimate);
        if (estimate->memoryFeasible) {
            ASSERT_GT(estimate->estimatedMinutes, 0.0);
            strictLimit = std::min(
                    strictLimit,
                    std::min(estimate->estimatedMinutes,
                             estimate->estimatedMinutesP90));
        }
    }
    ASSERT_TRUE(std::isfinite(strictLimit));
    strictLimit *= 0.5;

    const auto balanced = select(
            query, reference, 100000, strictLimit);
    for (const auto candidate : {
            anchorwave::AlignmentCandidate::SingletrackWfa,
            anchorwave::AlignmentCandidate::StandardWfa,
            anchorwave::AlignmentCandidate::MediumWfa,
            anchorwave::AlignmentCandidate::LowWfa,
            anchorwave::AlignmentCandidate::Ksw2ScoreCertified,
            anchorwave::AlignmentCandidate::Ksw2Full}) {
        const auto *estimate = findEstimate(balanced, candidate);
        ASSERT_NE(nullptr, estimate);
        if (!estimate->memoryFeasible) {
            continue;
        }
        EXPECT_FALSE(estimate->withinConfiguredTimeLimit)
                << anchorwave::alignmentCandidateName(candidate);
        EXPECT_FALSE(estimate->timeFeasible)
                << anchorwave::alignmentCandidateName(candidate);
        EXPECT_EQ(balanced.exactCandidates.end(),
                  std::find(balanced.exactCandidates.begin(),
                            balanced.exactCandidates.end(), candidate));
    }
    EXPECT_TRUE(balanced.exactCandidates.empty());
}

TEST(AlignmentAlgorithmSelector, WindowChangesResourcesNotSequenceEstimate) {
    std::string reference;
    reference.reserve(10000);
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int position = 0; position < 10000; ++position) {
        reference.push_back(bases[position % 4]);
    }
    std::string query = reference;
    for (std::size_t position = 0; position < query.size(); position += 17) {
        query[position] = query[position] == 'A' ? 'T' : 'A';
    }

    const auto small = select(query, reference, 100);
    const auto large = select(query, reference, 100000);
    EXPECT_EQ(small.profile.estimatedScore, large.profile.estimatedScore);
    EXPECT_EQ(small.profile.conservativeScore,
              large.profile.conservativeScore);
    EXPECT_EQ(small.profile.estimatedMismatchRate,
              large.profile.estimatedMismatchRate);
    EXPECT_LE(small.exactCandidates.size(), large.exactCandidates.size());
}

TEST(AlignmentAlgorithmSelector, CandidateLabelsRemainOutputCompatible) {
    EXPECT_STREQ("WAVEFRONT_SINGLETRACK", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::SingletrackWfa));
    EXPECT_STREQ("WAVEFRONT", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::StandardWfa));
    EXPECT_STREQ("WAVEFRONT_MEDIUM", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::MediumWfa));
    EXPECT_STREQ("WAVEFRONT_LOW", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::LowWfa));
    EXPECT_STREQ("KSW2_SCORE_CERTIFIED", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::Ksw2ScoreCertified));
    EXPECT_STREQ("MINIMAP2", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::Ksw2Full));
    EXPECT_STREQ("BANDED_MINIMAP2", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::Ksw2Banded));
    EXPECT_STREQ("SLIDING_WINDOW", anchorwave::alignmentCandidateName(
            anchorwave::AlignmentCandidate::SlidingWindow));
}

TEST(AlignmentAlgorithmSelector, BedLabelHidesInternalExactEngineModes) {
    EXPECT_EQ("WAVEFRONT", anchorwave::alignmentMethodBedLabel(
            "WAVEFRONT_SINGLETRACK"));
    EXPECT_EQ("WAVEFRONT", anchorwave::alignmentMethodBedLabel(
            "WAVEFRONT"));
    EXPECT_EQ("WAVEFRONT", anchorwave::alignmentMethodBedLabel(
            "WAVEFRONT_MEDIUM"));
    EXPECT_EQ("WAVEFRONT", anchorwave::alignmentMethodBedLabel(
            "WAVEFRONT_LOW"));
    EXPECT_EQ("MINIMAP2", anchorwave::alignmentMethodBedLabel(
            "KSW2_SCORE_CERTIFIED"));
    EXPECT_EQ("MINIMAP2", anchorwave::alignmentMethodBedLabel("MINIMAP2"));
    EXPECT_EQ("BANDED_MINIMAP2", anchorwave::alignmentMethodBedLabel(
            "BANDED_MINIMAP2"));
    EXPECT_EQ("SLIDING_WINDOW", anchorwave::alignmentMethodBedLabel(
            "SLIDING_WINDOW"));
    EXPECT_EQ("FILLING", anchorwave::alignmentMethodBedLabel("FILLING"));
    EXPECT_EQ("UNKNOWN", anchorwave::alignmentMethodBedLabel("UNKNOWN"));
}

TEST(AlignmentAlgorithmSelector,
     TraceSeparatesAttemptedCandidateFromRetainedBedMethod) {
    char path[] = "/tmp/anchorwave-alignment-trace-XXXXXX";
    const int descriptor = mkstemp(path);
    ASSERT_GE(descriptor, 0);
    close(descriptor);

    anchorwave::configureAlignmentTraceFile(path);
    anchorwave::AlignmentAttemptTrace record;
    record.intervalId = 1;
    record.attempt = 1;
    record.candidate = anchorwave::AlignmentCandidate::SlidingWindow;
    record.resultMethod = "BANDED_MINIMAP2";
    record.actualScore = -123;
    record.status = "completed";
    anchorwave::recordAlignmentAttemptTrace(record);
    anchorwave::configureAlignmentTraceFile("");

    std::ifstream input(path);
    std::string headerLine;
    std::string recordLine;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, headerLine)));
    ASSERT_TRUE(static_cast<bool>(std::getline(input, recordLine)));
    unlink(path);

    auto split = [](const std::string &line) {
        std::vector<std::string> fields;
        std::istringstream stream(line);
        std::string field;
        while (std::getline(stream, field, '\t')) {
            fields.push_back(field);
        }
        return fields;
    };
    const std::vector<std::string> header = split(headerLine);
    const std::vector<std::string> values = split(recordLine);
    ASSERT_EQ(header.size(), values.size());
    const auto candidate = std::find(
            header.begin(), header.end(), "candidate");
    const auto resultMethod = std::find(
            header.begin(), header.end(), "result_method");
    ASSERT_NE(header.end(), candidate);
    ASSERT_NE(header.end(), resultMethod);
    EXPECT_EQ("SLIDING_WINDOW",
              values[static_cast<std::size_t>(candidate - header.begin())]);
    EXPECT_EQ("BANDED_MINIMAP2",
              values[static_cast<std::size_t>(
                      resultMethod - header.begin())]);
}

TEST(AlignmentAlgorithmSelector, ReportsStrictTierTelemetry) {
    anchorwave::resetAlignmentSelectionTelemetry();
    const auto slidingOnly = select(std::string(20000, 'A'),
                                    std::string(20000, 'T'), 100);
    const auto bandedOnly = select(std::string(20000, 'A'),
                                   std::string(20000, 'T'), 30000,
                                   1.0e-12);
    const std::string exactSequence =
            "ACGTTGCAACCGGTTAACCGGTTACGTACGTTGCAACCGGTTA";
    const auto exact = select(exactSequence, exactSequence, 100000);
    anchorwave::recordAlignmentSelectionPlan(slidingOnly);
    anchorwave::recordAlignmentSelectionPlan(bandedOnly);
    anchorwave::recordAlignmentSelectionPlan(exact);
    anchorwave::recordExactAlignmentRuntimeFailure(true);
    anchorwave::recordExactAlignmentRuntimeFailure(false);
    anchorwave::recordBandedFallbackExecution();
    anchorwave::recordSlidingFallbackExecution();
    const auto telemetry =
            anchorwave::alignmentSelectionTelemetrySnapshot();
    EXPECT_EQ(3ULL, telemetry.evaluatedIntervals);
    EXPECT_EQ(1ULL, telemetry.exactTierIntervals);
    EXPECT_EQ(1ULL, telemetry.bandedOnlyIntervals);
    EXPECT_EQ(1ULL, telemetry.slidingOnlyIntervals);
    EXPECT_EQ(1ULL, telemetry.exactRuntimeMemoryFailures);
    EXPECT_EQ(1ULL, telemetry.exactRuntimeOtherFailures);
    EXPECT_EQ(1ULL, telemetry.bandedFallbackExecutions);
    EXPECT_EQ(1ULL, telemetry.slidingFallbackExecutions);
    EXPECT_EQ(2ULL, telemetry.mediumWfaMemoryRejects);
    EXPECT_EQ(2ULL, telemetry.lowWfaMemoryRejects);
    EXPECT_EQ(2ULL, telemetry.singletrackWfaMemoryRejects);
}

}  // namespace
