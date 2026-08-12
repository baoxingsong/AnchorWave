#include "src/myImportandFunction/WfaAlignment.h"
#include "src/myImportandFunction/AlignmentResourceScheduler.h"

extern "C" {
#include "WFA2-lib/wavefront/wfa.h"
}

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr uint64_t kDefaultWorkerMemoryBudget = 10000000000ULL;

struct MemoryProbeState {
    uint64_t calls = 0;
    uint64_t peakBytes = 0;
    uint64_t observedLimitBytes = 0;
    uint64_t rejectAboveBytes = std::numeric_limits<uint64_t>::max();
};

bool recordMemoryProbe(void *arguments,
                       uint64_t memoryUsedBytes,
                       uint64_t memoryLimitBytes,
                       int) {
    auto *const state = static_cast<MemoryProbeState *>(arguments);
    ++state->calls;
    state->peakBytes = std::max(state->peakBytes, memoryUsedBytes);
    state->observedLimitBytes = memoryLimitBytes;
    return memoryUsedBytes <= state->rejectAboveBytes;
}

std::string withoutGaps(std::string sequence) {
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '-'),
                   sequence.end());
    return sequence;
}

int64_t scoreAlignment(const std::string &query,
                       const std::string &reference,
                       int mismatch,
                       int gapOpen1,
                       int gapExtend1,
                       int gapOpen2,
                       int gapExtend2) {
    int64_t score = 0;
    std::size_t i = 0;
    while (i < query.size()) {
        if (query[i] == '-' || reference[i] == '-') {
            const bool queryGap = query[i] == '-';
            std::size_t end = i + 1;
            while (end < query.size() &&
                   (query[end] == '-') == queryGap &&
                   (reference[end] == '-') != queryGap) {
                ++end;
            }
            const int64_t length = static_cast<int64_t>(end - i);
            score -= std::min<int64_t>(gapOpen1 + gapExtend1 * length,
                                       gapOpen2 + gapExtend2 * length);
            i = end;
            continue;
        }
        if (query[i] != reference[i]) {
            score -= mismatch;
        }
        ++i;
    }
    return score;
}

int64_t standardWfaScore(const std::string &query,
                         const std::string &reference) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = 6;
    attributes.affine2p_penalties.gap_opening1 = 8;
    attributes.affine2p_penalties.gap_extension1 = 2;
    attributes.affine2p_penalties.gap_opening2 = 75;
    attributes.affine2p_penalties.gap_extension2 = 1;
    attributes.alignment_form.span = alignment_end2end;
    attributes.alignment_scope = compute_alignment;
    attributes.memory_mode = wavefront_memory_high;
    attributes.heuristic.strategy = wf_heuristic_none;
    attributes.system.verbose = 0;
    attributes.system.max_num_threads = 1;
    wavefront_aligner_t *aligner = wavefront_aligner_new(&attributes);
    if (aligner == nullptr) {
        return std::numeric_limits<int64_t>::min();
    }
    const int status = wavefront_align(
            aligner,
            reference.data(), static_cast<int>(reference.size()),
            query.data(), static_cast<int>(query.size()));
    const int64_t score = status == WF_STATUS_ALG_COMPLETED
                          ? aligner->cigar->score
                          : std::numeric_limits<int64_t>::min();
    wavefront_aligner_delete(aligner);
    return score;
}

TEST(WfaMemoryBudget, PreservesHistoricalWindowSquareSemantics) {
    EXPECT_EQ(10000000000ULL, anchorwave::wfaMemoryBudgetBytes(100000));
    EXPECT_EQ(3333333333ULL,
              anchorwave::biWfaComponentMemoryLimitBytes(
                      kDefaultWorkerMemoryBudget));
    EXPECT_EQ(2147483648ULL,
              anchorwave::standardWfaTrialMemoryBudgetBytes(
                      kDefaultWorkerMemoryBudget));
    EXPECT_EQ(1048576ULL,
              anchorwave::standardWfaTrialMemoryBudgetBytes(1048576ULL));
    EXPECT_EQ(0ULL, anchorwave::wfaMemoryBudgetBytes(0));
    EXPECT_EQ(0ULL, anchorwave::wfaMemoryBudgetBytes(-1));
    EXPECT_EQ(0ULL, anchorwave::biWfaComponentMemoryLimitBytes(0));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::wfaMemoryBudgetBytes(
                      std::numeric_limits<int64_t>::max()));
}

TEST(WfaMemoryBudget, DerivesAConservativeMonotonicBiWfaLeafScore) {
    EXPECT_EQ(250, wavefront_aligner_attr_default.system.biwfa_leaf_score);
    EXPECT_EQ(250, anchorwave::biWfaLeafScoreFromMemoryBudgetBytes(0));
    const int small = anchorwave::biWfaLeafScoreFromMemoryBudgetBytes(
            16ULL * 1024ULL * 1024ULL);
    const int medium = anchorwave::biWfaLeafScoreFromMemoryBudgetBytes(
            1024ULL * 1024ULL * 1024ULL);
    const int large = anchorwave::biWfaLeafScoreFromMemoryBudgetBytes(
            80ULL * 1024ULL * 1024ULL * 1024ULL);
    EXPECT_GE(small, 250);
    EXPECT_GE(medium, small);
    EXPECT_GE(large, medium);
    EXPECT_GT(large, 250);
    EXPECT_LE(large, 16384);
}

TEST(AlignmentResourceScheduler, KeepsHeadroomInsideTheProcessLimit) {
    const uint64_t eightyGiB = anchorwave::memoryLimitBytesFromGiB(80.0);
    EXPECT_EQ(85899345920ULL, eightyGiB);

    const auto plan = anchorwave::makeAlignmentResourcePlan(
            10, eightyGiB, 100000,
            anchorwave::memoryLimitBytesFromGiB(8.0));
    EXPECT_EQ(10, plan.requestedMaxThreads);
    EXPECT_EQ(10, plan.effectiveMaxThreads);
    EXPECT_EQ(10000000000ULL, plan.perWorkerMemoryBytes);
    EXPECT_EQ(10536870912ULL, plan.perWorkerReservationBytes);
    EXPECT_GE(plan.safetyReserveBytes,
              anchorwave::memoryLimitBytesFromGiB(1.0));
    EXPECT_LE(plan.safetyReserveBytes,
              anchorwave::memoryLimitBytesFromGiB(32.0));
    EXPECT_EQ(plan.maxProcessMemoryBytes,
              plan.baselineResidentBytes + plan.safetyReserveBytes +
                      plan.taskMemoryCapacityBytes);

    const auto sixWorkerPlan = anchorwave::makeAlignmentResourcePlan(
            6, eightyGiB, 100000,
            anchorwave::memoryLimitBytesFromGiB(8.0));
    EXPECT_EQ(6, sixWorkerPlan.effectiveMaxThreads);
    EXPECT_EQ(sixWorkerPlan.maxProcessMemoryBytes,
              sixWorkerPlan.baselineResidentBytes +
                      sixWorkerPlan.safetyReserveBytes +
                      sixWorkerPlan.taskMemoryCapacityBytes);

    const auto legacyPlan = anchorwave::makeAlignmentResourcePlan(
            10, 0, 100000,
            anchorwave::memoryLimitBytesFromGiB(8.0));
    EXPECT_EQ(10, legacyPlan.effectiveMaxThreads);

    EXPECT_THROW(anchorwave::makeAlignmentResourcePlan(
                         10, anchorwave::memoryLimitBytesFromGiB(9.0),
                         100000, 0), std::invalid_argument);
    EXPECT_THROW(anchorwave::makeAlignmentResourcePlan(
                         10, eightyGiB, 100000, eightyGiB),
                 std::runtime_error);
}

TEST(AlignmentResourceScheduler, AdmitsTasksByPredictedMemoryAndReplenishes) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; });

    std::vector<anchorwave::AlignmentMemoryReservation> reservations;
    for (int task = 0; task < 4; ++task) {
        reservations.push_back(scheduler.acquire(4 * gib));
        ASSERT_TRUE(reservations.back());
    }
    EXPECT_FALSE(scheduler.acquire(18 * gib));

    std::atomic<bool> fifthEntered(false);
    std::thread fifth([&]() {
        auto reservation = scheduler.acquire(4 * gib);
        fifthEntered.store(static_cast<bool>(reservation));
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    EXPECT_FALSE(fifthEntered.load());
    reservations.pop_back();
    fifth.join();
    EXPECT_TRUE(fifthEntered.load());

    const auto stats = scheduler.stats();
    EXPECT_EQ(5ULL, stats.reservationCount);
    EXPECT_EQ(1ULL, stats.waitedReservationCount);
    EXPECT_EQ(1ULL, stats.impossibleReservationCount);
    EXPECT_EQ(4, stats.peakConcurrentReservations);
    EXPECT_EQ(16 * gib, stats.peakReservedBytes);
    EXPECT_LE(stats.peakProjectedProcessBytes + plan.safetyReserveBytes,
              plan.maxProcessMemoryBytes);
}

TEST(AlignmentResourceScheduler, RejectsResidentDriftWithoutDeadlock) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            4, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 19 * gib; });
    EXPECT_FALSE(scheduler.acquire(gib));
    EXPECT_EQ(1ULL, scheduler.stats().impossibleReservationCount);
}

TEST(AlignmentResourceScheduler, ScopedInstallationRoutesWorkerReservations) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            4, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; });
    EXPECT_FALSE(anchorwave::alignmentMemorySchedulingEnabled());
    {
        anchorwave::ScopedAlignmentMemoryScheduler installed(scheduler);
        EXPECT_TRUE(anchorwave::alignmentMemorySchedulingEnabled());
        auto reservation = anchorwave::acquireAlignmentMemory(3 * gib);
        ASSERT_TRUE(reservation);
        EXPECT_EQ(3 * gib, reservation.reservedBytes());
    }
    EXPECT_FALSE(anchorwave::alignmentMemorySchedulingEnabled());
    EXPECT_EQ(1ULL, scheduler.stats().reservationCount);
}

TEST(AlignmentResourceScheduler,
     DistinguishesTemporaryPressureFromPermanentInfeasibility) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; });

    auto active = scheduler.acquire(12 * gib, 0.02);
    ASSERT_TRUE(active);
    const auto temporary = scheduler.tryAcquire(6 * gib, 0.01);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::TemporarilyUnavailable,
              temporary.admission);
    EXPECT_GT(scheduler.estimatedWaitMinutes(6 * gib), 0.0);
    EXPECT_LE(scheduler.estimatedWaitMinutes(6 * gib), 0.02);

    const auto impossible = scheduler.tryAcquire(18 * gib);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::PermanentlyInfeasible,
              impossible.admission);
    EXPECT_EQ(plan.taskMemoryCapacityBytes,
              scheduler.maximumSingleTaskReservationBytes());

    active = anchorwave::AlignmentMemoryReservation();
    auto admitted = scheduler.tryAcquire(6 * gib);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::Acquired,
              admitted.admission);
    EXPECT_TRUE(admitted);
    EXPECT_EQ(1ULL, scheduler.stats().temporaryReservationDeferrals);
}

TEST(AlignmentResourceScheduler,
     TreatsUntrackedResidentPressureAsTemporaryForNonblockingAdmission) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 19 * gib; });

    const auto result = scheduler.tryAcquire(gib / 2, 0.01);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::TemporarilyUnavailable,
              result.admission);
    EXPECT_GT(scheduler.estimatedWaitMinutes(gib / 2), 0.0);
    EXPECT_EQ(0ULL, scheduler.stats().impossibleReservationCount);
}

TEST(AlignmentResourceScheduler,
     ClassifiesRoundOnePersistentRssAsAStableRuntimeFloor) {
    // Exact values from B73/Mo17 100-Mb Round1 interval 2541.  The request is
    // below the startup task capacity, but 21.3 GB of retained RSS leaves it
    // 16.6 GB beyond the live -M envelope with no tracked releaser.
    const uint64_t processLimit = 91268055040ULL;
    const uint64_t baseline = 381091840ULL;
    const uint64_t observedRss = 21301166080ULL;
    const uint64_t request = 83455277678ULL;
    std::atomic<uint64_t> observed(observedRss);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, processLimit, 100000, baseline);
    ASSERT_LT(request, plan.taskMemoryCapacityBytes);

    anchorwave::AlignmentMemoryScheduler::Clock::time_point now;
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [&observed]() { return observed.load(); },
            [&now]() { return now; });
    anchorwave::AlignmentMemoryPressureState pressure;

    auto availability = scheduler.memoryAvailability(request, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              availability.availability);
    EXPECT_EQ(16575396462ULL, availability.shortfallBytes);
    EXPECT_NEAR(0.5, availability.waitMinutes, 1e-9);

    now += std::chrono::seconds(29);
    availability = scheduler.memoryAvailability(request, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              availability.availability);

    now += std::chrono::seconds(1);
    availability = scheduler.memoryAvailability(request, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::StableRuntimeFloor,
              availability.availability);
    EXPECT_TRUE(std::isinf(availability.waitMinutes));

    // The tier-boundary probe is instantaneous: a late RSS release revives
    // exact DP even after the bounded stable-floor conclusion was reached.
    observed.store(4ULL * 1024ULL * 1024ULL * 1024ULL);
    availability = scheduler.memoryAvailability(request, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::Ready,
              availability.availability);
    EXPECT_FALSE(pressure.initialized);
}

TEST(AlignmentResourceScheduler,
     FallingRssBecomesReadyInsteadOfAStableRuntimeFloor) {
    const uint64_t processLimit = 91268055040ULL;
    const uint64_t baseline = 381091840ULL;
    const uint64_t request = 83455277678ULL;
    std::atomic<uint64_t> observed(21301166080ULL);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, processLimit, 100000, baseline);
    anchorwave::AlignmentMemoryScheduler::Clock::time_point now;
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [&observed]() { return observed.load(); },
            [&now]() { return now; });
    anchorwave::AlignmentMemoryPressureState pressure;

    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              scheduler.memoryAvailability(request, &pressure).availability);
    now += std::chrono::seconds(25);
    observed.store(4ULL * 1024ULL * 1024ULL * 1024ULL);
    const auto ready = scheduler.memoryAvailability(request, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::Ready, ready.availability);
    EXPECT_EQ(0ULL, ready.shortfallBytes);
    EXPECT_FALSE(pressure.initialized);
}

TEST(AlignmentResourceScheduler,
     MeaningfulBlockedRssDropRestartsTheNoProgressClock) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const uint64_t processLimit = 91268055040ULL;
    const uint64_t baseline = 381091840ULL;
    const uint64_t request = 83455277678ULL;
    std::atomic<uint64_t> observed(21301166080ULL);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, processLimit, 100000, baseline);
    anchorwave::AlignmentMemoryScheduler::Clock::time_point now;
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [&observed]() { return observed.load(); },
            [&now]() { return now; });
    anchorwave::AlignmentMemoryPressureState pressure;

    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              scheduler.memoryAvailability(request, &pressure).availability);
    now += std::chrono::seconds(29);
    observed.fetch_sub(gib);
    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              scheduler.memoryAvailability(request, &pressure).availability);

    // The request is still blocked at the 30-second total-age boundary, but
    // the meaningful decline happened only one second ago.
    now += std::chrono::seconds(1);
    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              scheduler.memoryAvailability(request, &pressure).availability);
    now += std::chrono::seconds(4);
    EXPECT_EQ(anchorwave::MemoryAvailability::StableRuntimeFloor,
              scheduler.memoryAvailability(request, &pressure).availability);
}

TEST(AlignmentResourceScheduler,
     TrackedReleaseHasAFiniteWaitAndCannotBecomeAStableFloor) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler::Clock::time_point now;
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; },
            [&now]() { return now; });
    auto active = scheduler.acquire(12 * gib, 0.02);
    ASSERT_TRUE(active);
    anchorwave::AlignmentMemoryPressureState pressure;
    const auto availability = scheduler.memoryAvailability(6 * gib,
                                                            &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::TrackedWait,
              availability.availability);
    EXPECT_GT(availability.waitMinutes, 0.0);
    EXPECT_LE(availability.waitMinutes, 0.02);
    EXPECT_FALSE(pressure.initialized);
}

TEST(AlignmentResourceScheduler,
     UnquantifiedActiveCoolingHasBoundedFastExactGrace) {
    // Exact values from the failed Round1 interval 1951 admission.  Releasing
    // the tracked reservation alone does not satisfy the conservative KSW2
    // request, but the completing aligner can also return private allocator
    // pages.  The selector must get a bounded opportunity to wait for that
    // cooling instead of immediately starting the 15-minute BiWFA.
    const uint64_t processLimit = 91268055040ULL;
    const uint64_t baseline = 383631360ULL;
    const uint64_t observedRss = 22835740672ULL;
    const uint64_t activeBytes = 2846840173ULL;
    const uint64_t kswReservation = 77553317402ULL;
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, processLimit, 100000, baseline);
    anchorwave::AlignmentMemoryScheduler::Clock::time_point now;
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, []() { return observedRss; },
            [&now]() { return now; });
    auto active = scheduler.acquire(activeBytes, 1.0);
    ASSERT_TRUE(active);
    anchorwave::AlignmentMemoryPressureState pressure;
    const auto availability = scheduler.memoryAvailability(
            kswReservation, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::CoolingRss,
              availability.availability);
    EXPECT_DOUBLE_EQ(0.5, availability.waitMinutes);
    EXPECT_FALSE(pressure.initialized);
}

TEST(AlignmentResourceScheduler,
     RequestAboveStartupCapacityIsStaticallyInfeasible) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; });
    anchorwave::AlignmentMemoryPressureState pressure;
    const auto availability = scheduler.memoryAvailability(
            plan.taskMemoryCapacityBytes + 1, &pressure);
    EXPECT_EQ(anchorwave::MemoryAvailability::StaticInfeasible,
              availability.availability);
    EXPECT_EQ(1ULL, availability.shortfallBytes);
}

TEST(AlignmentResourceScheduler,
     AdaptiveReserveScalesAcrossSupportedThreadMemoryScenarios) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const std::pair<int, uint64_t> scenarios[] = {
            {12, 80}, {36, 120}, {96, 250},
            {128, 250}, {128, 500}, {128, 1000}};
    uint64_t previousLargeNodeReserve = 0;
    for (const auto &scenario : scenarios) {
        const auto plan = anchorwave::makeAlignmentResourcePlan(
                scenario.first, scenario.second * gib, 100000, gib);
        EXPECT_EQ(scenario.first, plan.effectiveMaxThreads);
        EXPECT_GE(plan.safetyReserveBytes, gib);
        EXPECT_LE(plan.safetyReserveBytes, 32 * gib);
        EXPECT_LT(plan.safetyReserveBytes,
                  plan.maxProcessMemoryBytes / 10 + gib);
        EXPECT_EQ(plan.maxProcessMemoryBytes,
                  plan.baselineResidentBytes + plan.safetyReserveBytes +
                          plan.taskMemoryCapacityBytes);
        if (scenario.second >= 500) {
            EXPECT_GE(plan.safetyReserveBytes, previousLargeNodeReserve);
            previousLargeNodeReserve = plan.safetyReserveBytes;
        }
    }
}

TEST(AlignmentResourceScheduler, SnapshotReportsDynamicCapacity) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; });
    auto reservation = scheduler.acquire(4 * gib, 0.1);
    ASSERT_TRUE(reservation);
    const auto snapshot = scheduler.snapshot();
    EXPECT_TRUE(snapshot.enabled);
    EXPECT_EQ(4 * gib, snapshot.activeReservedBytes);
    EXPECT_EQ(1, snapshot.activeReservations);
    EXPECT_EQ(plan.maxProcessMemoryBytes, snapshot.maxProcessMemoryBytes);
    EXPECT_LE(snapshot.immediatelyAvailableBytes,
              plan.taskMemoryCapacityBytes);
}

TEST(AlignmentResourceScheduler,
     RetainedResidentDriftRemainsInsideLaterAdmissions) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    std::atomic<uint64_t> observed(2 * gib);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [&observed]() { return observed.load(); });

    auto first = scheduler.acquire(4 * gib, 0.01);
    ASSERT_TRUE(first);
    // The completed task releases its token, but output strings and allocator
    // caches can remain resident. Preserve that drift in the next projection.
    observed.store(8 * gib);
    first = anchorwave::AlignmentMemoryReservation();
    EXPECT_EQ(6 * gib, scheduler.snapshot().untrackedResidentBytes);

    const auto blocked = scheduler.tryAcquire(11 * gib, 0.01);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::TemporarilyUnavailable,
              blocked.admission);

    // Once the consumer/caches release the resident state, a no-active-task
    // sample is allowed to lower the untracked floor and admission resumes.
    observed.store(2 * gib);
    auto admitted = scheduler.tryAcquire(11 * gib, 0.01);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::Acquired,
              admitted.admission);
    EXPECT_TRUE(admitted);
}

TEST(AlignmentResourceScheduler,
     EstimatedWaitKeepsUntrackedRssAfterTrackedCompletions) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    std::atomic<uint64_t> observed(2 * gib);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [&observed]() { return observed.load(); });

    auto early = scheduler.acquire(4 * gib, 0.01);
    auto late = scheduler.acquire(4 * gib, 0.04);
    ASSERT_TRUE(early);
    ASSERT_TRUE(late);
    // Six GiB are not represented by either active token.  Releasing only
    // the early token still cannot admit an eight-GiB request.
    observed.store(16 * gib);
    const double wait = scheduler.estimatedWaitMinutes(8 * gib);
    EXPECT_GT(wait, 0.02);
    EXPECT_LE(wait, 0.04);
}

TEST(AlignmentResourceScheduler,
     PreferredRequestDrainsMemoryWithoutYoungerTaskStarvation) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 2 * gib; });
    auto active = scheduler.acquire(12 * gib, 0.02);
    ASSERT_TRUE(active);

    std::atomic<bool> preferredAcquired(false);
    std::thread preferred([&]() {
        auto admission = scheduler.acquirePreferred(8 * gib, 0.01);
        preferredAcquired.store(static_cast<bool>(admission));
    });
    while (scheduler.stats().preferredReservationCount == 0) {
        std::this_thread::yield();
    }

    const auto younger = scheduler.tryAcquire(gib, 0.01);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::TemporarilyUnavailable,
              younger.admission);
    EXPECT_GT(scheduler.estimatedWaitMinutes(gib), 0.0);
    EXPECT_LE(scheduler.estimatedWaitMinutes(gib), 0.02);
    EXPECT_FALSE(preferredAcquired.load());
    active = anchorwave::AlignmentMemoryReservation();
    preferred.join();
    EXPECT_TRUE(preferredAcquired.load());

    const auto stats = scheduler.stats();
    EXPECT_EQ(1ULL, stats.preferredReservationCount);
    EXPECT_EQ(1ULL, stats.preferredWaitCount);
}

TEST(AlignmentResourceScheduler,
     PreferredAdmissionKeepsPersistentRssPressureTemporary) {
    const uint64_t gib = anchorwave::memoryLimitBytesFromGiB(1.0);
    const auto plan = anchorwave::makeAlignmentResourcePlan(
            12, 20 * gib, 100000, 2 * gib);
    anchorwave::AlignmentMemoryScheduler scheduler(
            plan, [gib]() { return 19 * gib; });

    const auto result = scheduler.acquirePreferred(gib / 2, 0.01);
    EXPECT_EQ(anchorwave::AlignmentMemoryAdmission::TemporarilyUnavailable,
              result.admission);
    EXPECT_EQ(0ULL, scheduler.stats().impossibleReservationCount);
    EXPECT_EQ(1ULL, scheduler.stats().temporaryReservationDeferrals);
}

TEST(AlignmentResourceScheduler, IncludesTransientStringsInsideReservation) {
    EXPECT_EQ(16ULL * 1024ULL * 1024ULL + 6ULL * 3000ULL,
              anchorwave::alignmentTaskTransientMemoryBytes(1000, 2000));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::alignmentTaskTransientMemoryBytes(
                      std::numeric_limits<uint64_t>::max(), 1));
}

TEST(AlignmentResourceScheduler, ReplenishesWorkerSlotsWithoutBusyWaiting) {
    anchorwave::DynamicThreadLimiter limiter(4);
    std::atomic_int active(0);
    std::vector<std::thread> workers;
    for (int task = 0; task < 12; ++task) {
        limiter.acquire();
        workers.emplace_back([&limiter, &active]() {
            ++active;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            --active;
            limiter.release();
        });
    }
    for (std::thread &worker : workers) {
        worker.join();
    }
    limiter.waitUntilIdle();
    EXPECT_EQ(0, active.load());
    EXPECT_EQ(4, limiter.maximumActiveThreads());
    EXPECT_EQ(4, limiter.peakActiveThreads());
}

TEST(BiWfaAlignment, AlignsIdenticalSequencesExactly) {
    const std::string sequence = "ACGTTGCAACGTTGCA";
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithBiWfa(
            sequence, sequence, kDefaultWorkerMemoryBudget,
            -6, -8, -2, -75, -1, result);

    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed, status);
    EXPECT_EQ(0, result.score);
    EXPECT_EQ(sequence, result.queryAlignment);
    EXPECT_EQ(sequence, result.referenceAlignment);
    EXPECT_EQ(10000000000ULL, result.memoryBudgetBytes);
    EXPECT_EQ(anchorwave::WfaAlgorithm::Bidirectional, result.algorithm);
}

TEST(StandardWfaAlignment, ProducesTheSameExactTwoPieceAffineScore) {
    const std::string query =
            "ACGTTGCAAAAAAACCGGTTAACCGGTTACGT";
    const std::string reference =
            "ACGTTGCACCGGTTCACCGGTTACGT";
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithStandardWfa(
            query, reference,
            anchorwave::standardWfaTrialMemoryBudgetBytes(
                    kDefaultWorkerMemoryBudget),
            -6, -8, -2, -30, -1, result);

    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed, status);
    EXPECT_EQ(anchorwave::WfaAlgorithm::Standard, result.algorithm);
    EXPECT_EQ(query, withoutGaps(result.queryAlignment));
    EXPECT_EQ(reference, withoutGaps(result.referenceAlignment));
    EXPECT_EQ(result.score,
              scoreAlignment(result.queryAlignment, result.referenceAlignment,
                             6, 8, 2, 30, 1));
}

TEST(SingletrackWfaAlignment, ProducesAValidExactTwoPieceAffineAlignment) {
    const std::string query =
            "ACGTTGCAAAAAAACCGGTTAACCGGTTACGT";
    const std::string reference =
            "ACGTTGCACCGGTTCACCGGTTACGT";
    anchorwave::WfaAlignmentResult standard;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      query, reference, 67108864ULL,
                      -6, -8, -2, -30, -1, standard));

    anchorwave::WfaAlignmentResult singletrack;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithSingletrackWfa(
                      query, reference, 67108864ULL,
                      -6, -8, -2, -30, -1, singletrack));
    EXPECT_EQ(anchorwave::WfaAlgorithm::Singletrack,
              singletrack.algorithm);
    EXPECT_EQ(standard.score, singletrack.score);
    EXPECT_EQ(query, withoutGaps(singletrack.queryAlignment));
    EXPECT_EQ(reference, withoutGaps(singletrack.referenceAlignment));
    EXPECT_EQ(singletrack.score,
              scoreAlignment(singletrack.queryAlignment,
                             singletrack.referenceAlignment,
                             6, 8, 2, 30, 1));
}

TEST(SingletrackWfaAlignment, ReapsAFailedTrialAndRemainsReusable) {
    anchorwave::WfaAlignmentResult result;
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::MemoryLimit,
              anchorwave::alignWithSingletrackWfa(
                      std::string(4096, 'A'), std::string(4096, 'T'),
                      1, -6, -8, -2, -75, -1, result));
    EXPECT_TRUE(result.queryAlignment.empty());
    EXPECT_TRUE(result.referenceAlignment.empty());

    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithSingletrackWfa(
                      "ACGT", "ACGT", 1,
                      -6, -8, -2, -75, -1, result));
    EXPECT_EQ(anchorwave::WfaAlgorithm::Singletrack, result.algorithm);
    EXPECT_EQ("ACGT", result.queryAlignment);
    EXPECT_EQ("ACGT", result.referenceAlignment);
}

TEST(SuccinctWfaAlignment,
     MediumAndLowMatchTheStandardWfaScoreAndCigar) {
    const std::string query =
            "ACGTTGCAAAAAAACCGGTTAACCGGTTACGT";
    const std::string reference =
            "ACGTTGCACCGGTTCACCGGTTACGT";
    constexpr uint64_t budget = 67108864ULL;

    anchorwave::WfaAlignmentResult high;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      query, reference, budget,
                      -6, -8, -2, -30, -1, high));

    anchorwave::WfaAlignmentResult medium;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithMediumWfa(
                      query, reference, budget,
                      -6, -8, -2, -30, -1, medium));
    EXPECT_EQ(anchorwave::WfaAlgorithm::Medium, medium.algorithm);

    anchorwave::WfaAlignmentResult low;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithLowWfa(
                      query, reference, budget,
                      -6, -8, -2, -30, -1, low));
    EXPECT_EQ(anchorwave::WfaAlgorithm::Low, low.algorithm);

    EXPECT_EQ(high.score, medium.score);
    EXPECT_EQ(high.score, low.score);
    EXPECT_EQ(high.queryAlignment, medium.queryAlignment);
    EXPECT_EQ(high.referenceAlignment, medium.referenceAlignment);
    EXPECT_EQ(high.queryAlignment, low.queryAlignment);
    EXPECT_EQ(high.referenceAlignment, low.referenceAlignment);
}

TEST(SuccinctWfaAlignment, ReapsMemoryLimitedModesAndRemainsReusable) {
    const std::string query(4096, 'A');
    const std::string reference(4096, 'T');
    anchorwave::WfaAlignmentResult medium;
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::MemoryLimit,
              anchorwave::alignWithMediumWfa(
                      query, reference, 1,
                      -6, -8, -2, -75, -1, medium));
    EXPECT_TRUE(medium.queryAlignment.empty());
    EXPECT_TRUE(medium.referenceAlignment.empty());

    anchorwave::WfaAlignmentResult low;
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::MemoryLimit,
              anchorwave::alignWithLowWfa(
                      query, reference, 1,
                      -6, -8, -2, -75, -1, low));
    EXPECT_TRUE(low.queryAlignment.empty());
    EXPECT_TRUE(low.referenceAlignment.empty());

    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithMediumWfa(
                      "ACGT", "ACGT", 1,
                      -6, -8, -2, -75, -1, medium));
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithLowWfa(
                      "ACGT", "ACGT", 1,
                      -6, -8, -2, -75, -1, low));
}

TEST(StandardWfaAlignment, ReapsAFailedTrialAndRemainsReusable) {
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithStandardWfa(
            std::string(4096, 'A'), std::string(4096, 'T'),
            1, -6, -8, -2, -75, -1, result);

    EXPECT_EQ(anchorwave::WfaAlignmentStatus::MemoryLimit, status);
    EXPECT_TRUE(result.queryAlignment.empty());
    EXPECT_TRUE(result.referenceAlignment.empty());

    const auto recoveryStatus = anchorwave::alignWithStandardWfa(
            "ACGT", "ACGT", 1, -6, -8, -2, -75, -1, result);
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed, recoveryStatus);
    EXPECT_EQ(anchorwave::WfaAlgorithm::Standard, result.algorithm);
    EXPECT_EQ("ACGT", result.queryAlignment);
    EXPECT_EQ("ACGT", result.referenceAlignment);
}

TEST(AdaptiveWfaAlignment, UsesSingletrackWfaAsTheFirstFastPath) {
    const std::string reference = "ACGTTGCACCGGTTAACCGGTTACGT";
    const std::string query = "ACGTTGCAACCGGTTAACCGGTTACGT";
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithAdaptiveWfa(
            query, reference, kDefaultWorkerMemoryBudget,
            -6, -8, -2, -75, -1, result);

    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed, status);
    EXPECT_EQ(anchorwave::WfaAlgorithm::Singletrack, result.algorithm);
    EXPECT_EQ(kDefaultWorkerMemoryBudget, result.memoryBudgetBytes);
    EXPECT_EQ(query, withoutGaps(result.queryAlignment));
    EXPECT_EQ(reference, withoutGaps(result.referenceAlignment));
}

TEST(WfaAlignment, KeepsAllMemoryModeThreadCachesIndependent) {
    const std::string reference = "ACGTTGCACCGGTTAACCGGTTACGT";
    const std::string query = "ACGTTGCAACCGGTTAACCGGTTACGT";
    for (int iteration = 0; iteration < 10; ++iteration) {
        anchorwave::WfaAlignmentResult singletrack;
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
                  anchorwave::alignWithSingletrackWfa(
                          query, reference, 67108864ULL,
                          -6, -8, -2, -75, -1, singletrack));
        EXPECT_EQ(anchorwave::WfaAlgorithm::Singletrack,
                  singletrack.algorithm);

        anchorwave::WfaAlignmentResult standard;
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
                  anchorwave::alignWithStandardWfa(
                          query, reference, 67108864ULL,
                          -6, -8, -2, -75, -1, standard));
        EXPECT_EQ(anchorwave::WfaAlgorithm::Standard, standard.algorithm);

        anchorwave::WfaAlignmentResult medium;
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
                  anchorwave::alignWithMediumWfa(
                          query, reference, 67108864ULL,
                          -6, -8, -2, -75, -1, medium));
        EXPECT_EQ(anchorwave::WfaAlgorithm::Medium, medium.algorithm);

        anchorwave::WfaAlignmentResult low;
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
                  anchorwave::alignWithLowWfa(
                          query, reference, 67108864ULL,
                          -6, -8, -2, -75, -1, low));
        EXPECT_EQ(anchorwave::WfaAlgorithm::Low, low.algorithm);

        anchorwave::WfaAlignmentResult bidirectional;
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
                  anchorwave::alignWithBiWfa(
                          query, reference, 67108864ULL,
                          -6, -8, -2, -75, -1, bidirectional));
        EXPECT_EQ(anchorwave::WfaAlgorithm::Bidirectional,
                  bidirectional.algorithm);
        EXPECT_EQ(standard.score, singletrack.score);
        EXPECT_EQ(standard.score, medium.score);
        EXPECT_EQ(standard.score, low.score);
        EXPECT_EQ(standard.score, bidirectional.score);
    }
}

TEST(BiWfaAlignment, ProducesAValidTwoPieceAffineAlignment) {
    const std::string query =
            "ACGTTGCAAAAAAACCGGTTAACCGGTTACGT";
    const std::string reference =
            "ACGTTGCACCGGTTCACCGGTTACGT";
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithBiWfa(
            query, reference, kDefaultWorkerMemoryBudget,
            -6, -8, -2, -30, -1, result);

    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed, status);
    ASSERT_EQ(result.queryAlignment.size(), result.referenceAlignment.size());
    EXPECT_EQ(query, withoutGaps(result.queryAlignment));
    EXPECT_EQ(reference, withoutGaps(result.referenceAlignment));
    EXPECT_EQ(result.score,
              scoreAlignment(result.queryAlignment, result.referenceAlignment,
                             6, 8, 2, 30, 1));
}

TEST(BiWfaAlignment, ConfigurableExactLeafThresholdPreservesOptimalScore) {
    std::string reference(2000, 'A');
    std::string query = reference;
    for (std::size_t i = 3; i < query.size(); i += 5) {
        query[i] = 'C';
    }

    anchorwave::WfaExecutionOptions legacyOptions;
    legacyOptions.biWfaLeafScore = 250;
    anchorwave::WfaAlignmentResult legacy;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithBiWfa(
                      query, reference, 256ULL * 1024ULL * 1024ULL,
                      -6, -8, -2, -75, -1, legacy, legacyOptions));
    EXPECT_EQ(250, legacy.configuredBiWfaLeafScore);
    EXPECT_GT(legacy.biWfaLeafAlignments, 0ULL);

    anchorwave::WfaExecutionOptions largerLeafOptions;
    largerLeafOptions.biWfaLeafScore = 4000;
    anchorwave::WfaAlignmentResult largerLeaf;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithBiWfa(
                      query, reference, 256ULL * 1024ULL * 1024ULL,
                      -6, -8, -2, -75, -1,
                      largerLeaf, largerLeafOptions));
    EXPECT_EQ(4000, largerLeaf.configuredBiWfaLeafScore);
    EXPECT_EQ(legacy.score, largerLeaf.score);
    EXPECT_EQ(query, withoutGaps(largerLeaf.queryAlignment));
    EXPECT_EQ(reference, withoutGaps(largerLeaf.referenceAlignment));
    EXPECT_LE(largerLeaf.biWfaLeafAlignments,
              legacy.biWfaLeafAlignments);
    EXPECT_GE(largerLeaf.biWfaMaxLeafScore,
              legacy.biWfaMaxLeafScore);
    EXPECT_EQ(largerLeaf.score,
              scoreAlignment(largerLeaf.queryAlignment,
                             largerLeaf.referenceAlignment,
                             6, 8, 2, 75, 1));
}

TEST(WfaExecution, ConfiguresTailThreadsWithoutChangingTheDefault) {
    anchorwave::WfaAlignmentResult defaultResult;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      "ACGTACGT", "ACGTACGT", 64ULL * 1024ULL * 1024ULL,
                      -6, -8, -2, -75, -1, defaultResult));
    EXPECT_EQ(1, defaultResult.configuredMaxThreads);

    anchorwave::WfaExecutionOptions tailOptions;
    tailOptions.maxNumThreads = 8;
    tailOptions.minOffsetsPerThread = 256;
    anchorwave::WfaAlignmentResult tailResult;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      "ACGTACGT", "ACGTACGT", 64ULL * 1024ULL * 1024ULL,
                      -6, -8, -2, -75, -1, tailResult, tailOptions));
    EXPECT_EQ(8, tailResult.configuredMaxThreads);
    EXPECT_EQ(defaultResult.score, tailResult.score);

    tailOptions.maxNumThreads = 0;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      "ACGT", "ACGT", 64ULL * 1024ULL * 1024ULL,
                      -6, -8, -2, -75, -1, tailResult, tailOptions));
    EXPECT_EQ(1, tailResult.configuredMaxThreads);
}

TEST(WfaExecution, OneTwoFourAndEightThreadsKeepTheExactScore) {
    std::mt19937 generator(20260810);
    std::uniform_int_distribution<int> base(0, 3);
    constexpr char bases[] = "ACGT";
    std::string reference;
    reference.reserve(1500);
    for (int i = 0; i < 1500; ++i) {
        reference.push_back(bases[base(generator)]);
    }
    std::string query = reference;
    for (std::size_t i = 1; i < query.size(); i += 4) {
        char replacement = query[i];
        while (replacement == query[i]) {
            replacement = bases[base(generator)];
        }
        query[i] = replacement;
    }

    int64_t exactScore = std::numeric_limits<int64_t>::min();
    const int threadCounts[] = {1, 2, 4, 8};
    for (const int threads : threadCounts) {
        anchorwave::WfaExecutionOptions options;
        options.maxNumThreads = threads;
        options.minOffsetsPerThread = 1;
        anchorwave::WfaAlignmentResult result;
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
                  anchorwave::alignWithStandardWfa(
                          query, reference, 512ULL * 1024ULL * 1024ULL,
                          -6, -8, -2, -75, -1, result, options));
        if (exactScore == std::numeric_limits<int64_t>::min()) {
            exactScore = result.score;
        }
        EXPECT_EQ(exactScore, result.score) << "threads=" << threads;
        EXPECT_EQ(result.score,
                  scoreAlignment(result.queryAlignment,
                                 result.referenceAlignment,
                                 6, 8, 2, 75, 1));
        EXPECT_EQ(threads, result.configuredMaxThreads);
        EXPECT_EQ(anchorwave::wfaParallelSupportEnabled() ? threads : 1,
                  result.maxThreadsUsed)
                << "threads=" << threads;
    }
}

TEST(WfaExecution, CooperativeMemoryProbeReportsAndCanVetoGrowth) {
    const uint64_t budget = 64ULL * 1024ULL * 1024ULL;
    std::string reference(1200, 'A');
    std::string query = reference;
    for (std::size_t i = 0; i < query.size(); i += 12) {
        query[i] = 'C';
    }

    MemoryProbeState observed;
    anchorwave::WfaExecutionOptions observeOptions;
    observeOptions.memoryProbeScoreInterval = 1;
    observeOptions.memoryProbe = &recordMemoryProbe;
    observeOptions.memoryProbeArguments = &observed;
    anchorwave::WfaAlignmentResult completed;
    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      query, reference, budget,
                      -6, -8, -2, -75, -1,
                      completed, observeOptions));
    EXPECT_GT(observed.calls, 0ULL);
    EXPECT_EQ(budget, observed.observedLimitBytes);
    EXPECT_EQ(observed.calls, completed.memoryProbeCount);
    EXPECT_GE(completed.memoryPeakBytes, observed.peakBytes);

    MemoryProbeState vetoed;
    vetoed.rejectAboveBytes = 0;
    anchorwave::WfaExecutionOptions vetoOptions;
    vetoOptions.biWfaLeafScore = 250;
    vetoOptions.memoryProbeScoreInterval = 1;
    vetoOptions.memoryProbe = &recordMemoryProbe;
    vetoOptions.memoryProbeArguments = &vetoed;
    anchorwave::WfaAlignmentResult stopped;
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::MemoryLimit,
              anchorwave::alignWithBiWfa(
                      query, reference, budget,
                      -6, -8, -2, -75, -1,
                      stopped, vetoOptions));
    EXPECT_GT(vetoed.calls, 0ULL);
    EXPECT_EQ(budget, vetoed.observedLimitBytes);
    EXPECT_TRUE(stopped.queryAlignment.empty());
    EXPECT_TRUE(stopped.referenceAlignment.empty());
}

TEST(WfaExecution, CooperativeRuntimeLimitReturnsNoPartialAlignment) {
    const std::string query(12000, 'A');
    const std::string reference(12000, 'T');
    anchorwave::WfaExecutionOptions options;
    options.memoryProbeScoreInterval = 1;
    options.maximumRuntimeMilliseconds = 1;
    anchorwave::WfaAlignmentResult stopped;
    const auto status = anchorwave::alignWithStandardWfa(
            query, reference, kDefaultWorkerMemoryBudget,
            -6, -8, -2, -75, -1, stopped, options);
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::TimeLimit, status);
    EXPECT_TRUE(stopped.timeLimitExceeded);
    EXPECT_TRUE(stopped.queryAlignment.empty());
    EXPECT_TRUE(stopped.referenceAlignment.empty());

    anchorwave::WfaAlignmentResult recovery;
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed,
              anchorwave::alignWithStandardWfa(
                      "ACGT", "ACGT", kDefaultWorkerMemoryBudget,
                      -6, -8, -2, -75, -1, recovery));
}

TEST(BiWfaAlignment, MatchesStandardWfaOptimalScoreOnSeededVariants) {
    std::mt19937 generator(20260807);
    std::uniform_int_distribution<int> base(0, 3);
    std::uniform_int_distribution<int> action(0, 99);
    constexpr char bases[] = "ACGT";
    for (int example = 0; example < 100; ++example) {
        std::string reference;
        const int referenceLength = 20 + example % 81;
        reference.reserve(referenceLength);
        for (int i = 0; i < referenceLength; ++i) {
            reference.push_back(bases[base(generator)]);
        }

        std::string query;
        query.reserve(reference.size() + 16);
        for (char nucleotide : reference) {
            const int mutation = action(generator);
            if (mutation < 8) {
                continue;  // deletion from the query
            }
            if (mutation < 18) {
                query.push_back(bases[base(generator)]);  // insertion
            }
            query.push_back(mutation >= 18 && mutation < 33
                            ? bases[base(generator)] : nucleotide);
        }

        anchorwave::WfaAlignmentResult result;
        const auto status = anchorwave::alignWithBiWfa(
                query, reference, kDefaultWorkerMemoryBudget,
                -6, -8, -2, -75, -1, result);
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed, status)
                << "example " << example;
        EXPECT_EQ(standardWfaScore(query, reference), result.score)
                << "example " << example;
        EXPECT_EQ(query, withoutGaps(result.queryAlignment));
        EXPECT_EQ(reference, withoutGaps(result.referenceAlignment));
    }
}

TEST(BiWfaAlignment, HandlesStrongLengthImbalanceWithoutAdmissionRules) {
    const std::string query = "ACGTACGTACGTACGTACGTACGT";
    const std::string reference = "ACGT";
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithBiWfa(
            query, reference, kDefaultWorkerMemoryBudget,
            -6, -8, -2, -75, -1, result);

    ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed, status);
    EXPECT_EQ(query, withoutGaps(result.queryAlignment));
    EXPECT_EQ(reference, withoutGaps(result.referenceAlignment));
}

TEST(BiWfaAlignment, ReusesTheThreadLocalAlignerAcrossCalls) {
    const std::string reference = "ACGTTGCACCGGTTAACCGGTTACGT";
    const std::string query = "ACGTTGCAACCGGTTAACCGGTTACGT";
    for (int iteration = 0; iteration < 25; ++iteration) {
        anchorwave::WfaAlignmentResult result;
        const auto status = anchorwave::alignWithBiWfa(
                query, reference, kDefaultWorkerMemoryBudget,
                -6, -8, -2, -75, -1, result);
        ASSERT_EQ(anchorwave::WfaAlignmentStatus::Completed, status);
        EXPECT_EQ(query, withoutGaps(result.queryAlignment));
        EXPECT_EQ(reference, withoutGaps(result.referenceAlignment));
    }
}

TEST(BiWfaAlignment, IsolatedThreadLocalAlignersWorkWithSixWorkers) {
    constexpr int workerCount = 6;
    std::vector<int> workerSucceeded(workerCount, 0);
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (int worker = 0; worker < workerCount; ++worker) {
        workers.emplace_back([worker, &workerSucceeded]() {
            const std::string reference =
                    "ACGTTGCACCGGTTAACCGGTTACGT";
            const std::string query =
                    "ACGTTGCA" + std::string(worker + 1, 'A') +
                    "CCGGTTAACCGGTTACGT";
            for (int iteration = 0; iteration < 10; ++iteration) {
                anchorwave::WfaAlignmentResult result;
                const auto status = anchorwave::alignWithBiWfa(
                        query, reference, kDefaultWorkerMemoryBudget,
                        -6, -8, -2, -75, -1, result);
                if (status != anchorwave::WfaAlignmentStatus::Completed ||
                    withoutGaps(result.queryAlignment) != query ||
                    withoutGaps(result.referenceAlignment) != reference) {
                    return;
                }
            }
            workerSucceeded[worker] = 1;
        });
    }
    for (std::thread &worker : workers) {
        worker.join();
    }
    EXPECT_EQ(std::vector<int>(workerCount, 1), workerSucceeded);
}

TEST(BiWfaAlignment, RejectsAnInvalidMemoryBudget) {
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithBiWfa(
            "ACGT", "ACGT", 0, -6, -8, -2, -75, -1, result);

    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Failed, status);
    EXPECT_EQ(0ULL, result.memoryBudgetBytes);
    EXPECT_TRUE(result.queryAlignment.empty());
    EXPECT_TRUE(result.referenceAlignment.empty());
}

TEST(BiWfaAlignment, ReportsMemoryLimitWithoutPartialAlignment) {
    // A one-byte aggregate budget is intentionally impossible for this
    // high-score alignment and exercises the WF_STATUS_OOM path.
    const std::string query(4096, 'A');
    const std::string reference(4096, 'T');
    anchorwave::WfaAlignmentResult result;
    const auto status = anchorwave::alignWithBiWfa(
            query, reference, 1, -6, -8, -2, -75, -1, result);

    EXPECT_EQ(anchorwave::WfaAlignmentStatus::MemoryLimit, status);
    EXPECT_EQ(1ULL, result.memoryBudgetBytes);
    EXPECT_TRUE(result.queryAlignment.empty());
    EXPECT_TRUE(result.referenceAlignment.empty());

    // The thread-local object is reaped after OOM and must remain reusable.
    const auto recoveryStatus = anchorwave::alignWithBiWfa(
            "ACGT", "ACGT", 1, -6, -8, -2, -75, -1, result);
    EXPECT_EQ(anchorwave::WfaAlignmentStatus::Completed, recoveryStatus);
    EXPECT_EQ("ACGT", result.queryAlignment);
    EXPECT_EQ("ACGT", result.referenceAlignment);
}

}  // namespace
