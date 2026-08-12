#include "src/service/AnchorTaskExecutor.h"
#include "src/service/ParallelGapBatchScheduler.h"

#include "gtest/gtest.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

TEST(AnchorTaskExecutorTest, ExecutesEveryTaskWithinGlobalThreadLimit) {
    anchorwave::AnchorTaskExecutor executor(4);
    std::mutex mutex;
    std::condition_variable condition;
    bool release = false;
    std::atomic<int> entered(0);
    std::atomic<int> completed(0);

    for (int task = 0; task < 12; ++task) {
        executor.submit(static_cast<uint64_t>(task), [&]() {
            ++entered;
            std::unique_lock<std::mutex> lock(mutex);
            condition.wait(lock, [&]() { return release; });
            ++completed;
        });
    }

    while (entered.load() < 4) {
        std::this_thread::yield();
    }
    {
        std::lock_guard<std::mutex> lock(mutex);
        release = true;
    }
    condition.notify_all();
    executor.waitForIdle();

    EXPECT_EQ(12, completed.load());
    EXPECT_EQ(4, executor.workerCount());
    EXPECT_EQ(4, executor.peakActiveTasks());
}

TEST(AnchorTaskExecutorTest, PullsLargestQueuedCostFirst) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::mutex mutex;
    std::condition_variable condition;
    bool blockerEntered = false;
    bool releaseBlocker = false;
    std::vector<std::string> order;

    executor.submit(1000, [&]() {
        std::unique_lock<std::mutex> lock(mutex);
        blockerEntered = true;
        order.push_back("blocker");
        condition.notify_all();
        condition.wait(lock, [&]() { return releaseBlocker; });
    });
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return blockerEntered; });
    }

    executor.submit(1, [&]() {
        std::lock_guard<std::mutex> lock(mutex);
        order.push_back("low");
    });
    executor.submit(100, [&]() {
        std::lock_guard<std::mutex> lock(mutex);
        order.push_back("high");
    });
    executor.submit(10, [&]() {
        std::lock_guard<std::mutex> lock(mutex);
        order.push_back("medium");
    });
    {
        std::lock_guard<std::mutex> lock(mutex);
        releaseBlocker = true;
    }
    condition.notify_all();
    executor.waitForIdle();

    ASSERT_EQ(4U, order.size());
    EXPECT_EQ("blocker", order[0]);
    EXPECT_EQ("high", order[1]);
    EXPECT_EQ("medium", order[2]);
    EXPECT_EQ("low", order[3]);
}

TEST(AnchorTaskExecutorTest, PropagatesTaskExceptionAndRemainsUsable) {
    anchorwave::AnchorTaskExecutor executor(2);
    executor.submit(1, []() {
        throw std::runtime_error("anchor task failed");
    });
    EXPECT_THROW(executor.waitForIdle(), std::runtime_error);

    std::atomic<int> completed(0);
    executor.submit(1, [&]() { ++completed; });
    EXPECT_NO_THROW(executor.waitForIdle());
    EXPECT_EQ(1, completed.load());
}

TEST(AnchorTaskExecutorTest, UsesQuadraticSaturatingCost) {
    EXPECT_EQ(0ULL, anchorwave::anchorTaskEstimatedCost(0));
    EXPECT_EQ(100ULL, anchorwave::anchorTaskEstimatedCost(10));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::anchorTaskEstimatedCost(
                      std::numeric_limits<std::size_t>::max()));
}

TEST(AnchorTaskExecutorTest, WorkerCooperativelyCompletesNestedTaskGroup) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::atomic<int> childrenCompleted(0);

    for (int parent = 0; parent < 2; ++parent) {
        executor.submit(100, [&]() {
            anchorwave::AnchorTaskGroup group;
            for (int child = 0; child < 4; ++child) {
                executor.submit(group, 10, [&]() {
                    ++childrenCompleted;
                });
            }
            executor.waitForGroup(group);
        });
    }

    executor.waitForIdle();
    EXPECT_EQ(8, childrenCompleted.load());
    EXPECT_LE(executor.peakActiveTasks(), 2);
}

TEST(AnchorTaskExecutorTest, UsesSaturatingGapCost) {
    EXPECT_EQ(35u, anchorwave::anchorGapTaskEstimatedCost(5, 7));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(),
              anchorwave::anchorGapTaskEstimatedCost(
                      std::numeric_limits<std::size_t>::max(), 2));
}

TEST(AnchorTaskExecutorTest, CanBoundPendingNestedTasks) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::atomic<int> completed(0);
    executor.submit(100, [&]() {
        anchorwave::AnchorTaskGroup group;
        for (int child = 0; child < 5; ++child) {
            executor.submit(group, static_cast<uint64_t>(child), [&]() {
                ++completed;
            });
            executor.waitUntilGroupSizeAtMost(group, 1);
            EXPECT_GE(completed.load(), child);
        }
        executor.waitForGroup(group);
    });
    executor.waitForIdle();
    EXPECT_EQ(5, completed.load());
}

TEST(AnchorTaskExecutorTest, DeferrableTaskParksAndReleasesItsWorker) {
    anchorwave::AnchorTaskExecutor executor(1);
    anchorwave::AnchorTaskGroup group;
    std::mutex mutex;
    std::condition_variable condition;
    bool blockerEntered = false;
    bool releaseBlocker = false;
    std::vector<std::string> order;
    std::atomic<int> attempts(0);

    executor.submit(1000, [&]() {
        std::unique_lock<std::mutex> lock(mutex);
        blockerEntered = true;
        condition.notify_all();
        condition.wait(lock, [&]() { return releaseBlocker; });
    });
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return blockerEntered; });
    }

    executor.submitDeferrable(group, 100, [&]() {
        const int attempt = ++attempts;
        EXPECT_TRUE(anchorwave::currentAnchorTaskCanDefer());
        EXPECT_EQ(static_cast<uint64_t>(attempt - 1),
                  anchorwave::currentAnchorTaskDeferralCount());
        {
            std::lock_guard<std::mutex> lock(mutex);
            order.push_back(attempt == 1 ? "defer-1" : "defer-2");
        }
        if (attempt == 1) {
            throw anchorwave::AlignmentTaskDeferred(
                    std::chrono::milliseconds(40));
        }
    });
    executor.submit(group, 1, [&]() {
        EXPECT_FALSE(anchorwave::currentAnchorTaskCanDefer());
        std::lock_guard<std::mutex> lock(mutex);
        order.push_back("ordinary");
    });
    {
        std::lock_guard<std::mutex> lock(mutex);
        releaseBlocker = true;
    }
    condition.notify_all();

    executor.waitForGroup(group);
    executor.waitForIdle();
    ASSERT_EQ(3U, order.size());
    EXPECT_EQ("defer-1", order[0]);
    EXPECT_EQ("ordinary", order[1]);
    EXPECT_EQ("defer-2", order[2]);
    EXPECT_EQ(2, attempts.load());
}

TEST(AnchorTaskExecutorTest, RetainsPreparedStateAcrossDeferral) {
    anchorwave::AnchorTaskExecutor executor(1);
    anchorwave::AnchorTaskGroup group;
    std::atomic<int> attempts(0);
    executor.submitDeferrable(group, 10, [&]() {
        ++attempts;
        auto retained = anchorwave::currentAnchorTaskRetainedState();
        if (!retained) {
            retained = std::make_shared<int>(73);
            anchorwave::setCurrentAnchorTaskRetainedState(retained);
            throw anchorwave::AlignmentTaskDeferred(
                    std::chrono::milliseconds(1));
        }
        EXPECT_EQ(73, *std::static_pointer_cast<int>(retained));
    });
    executor.waitForGroup(group);
    EXPECT_EQ(2, attempts.load());
}

TEST(AnchorTaskExecutorTest, FullWorkerBatchIsNotPrematurelyTail) {
    anchorwave::AnchorTaskExecutor executor(16);
    std::mutex mutex;
    std::condition_variable condition;
    bool release = false;
    std::atomic<int> entered(0);
    for (int task = 0; task < 16; ++task) {
        executor.submit(100, [&]() {
            ++entered;
            std::unique_lock<std::mutex> lock(mutex);
            condition.wait(lock, [&]() { return release; });
        });
    }
    while (entered.load() < 16) {
        std::this_thread::yield();
    }
    const auto load = executor.loadSnapshot();
    EXPECT_EQ(16ULL, load.outstandingTasks);
    EXPECT_FALSE(load.tailPhase);
    {
        std::lock_guard<std::mutex> lock(mutex);
        release = true;
    }
    condition.notify_all();
    executor.waitForIdle();
}

TEST(AnchorTaskExecutorTest, PlannedRollingBacklogParticipatesInTailState) {
    anchorwave::AnchorTaskExecutor executor(16);
    const std::vector<uint64_t> costs{10, 20, 30, 40, 50};
    executor.registerPlannedTaskCosts(costs);

    const auto load = executor.loadSnapshot();
    EXPECT_EQ(5ULL, load.plannedTasks);
    EXPECT_EQ(5ULL, load.outstandingTasks);
    EXPECT_EQ(150ULL, load.outstandingEstimatedCost);
    EXPECT_EQ(50ULL, load.criticalEstimatedCost);
    EXPECT_FALSE(load.tailPhase);

    executor.cancelPlannedTaskCosts(costs);
    EXPECT_NO_THROW(executor.waitForIdle());
    EXPECT_EQ(0ULL, executor.loadSnapshot().plannedTasks);
}

TEST(AnchorTaskExecutorTest,
     OrderedHeadSeparatesAdmissionTailFromGlobalFutureBacklog) {
    anchorwave::AnchorTaskExecutor executor(16);
    const std::vector<uint64_t> costs(100, 10);
    executor.registerPlannedTaskCosts(costs);

    auto load = executor.loadSnapshot();
    EXPECT_EQ(100ULL, load.globalFutureTasks);
    EXPECT_EQ(0ULL, load.schedulableTasks);
    EXPECT_FALSE(load.globalTailPhase);
    EXPECT_FALSE(load.admissionTailPhase);

    executor.orderedHeadWaitStarted();
    load = executor.loadSnapshot();
    EXPECT_EQ(1ULL, load.blockedOrderedHeads);
    EXPECT_FALSE(load.globalTailPhase);
    EXPECT_TRUE(load.admissionTailPhase);

    executor.orderedHeadWaitFinished();
    EXPECT_FALSE(executor.loadSnapshot().admissionTailPhase);
    executor.cancelPlannedTaskCosts(costs);
}

TEST(AnchorTaskExecutorTest,
     OtherRunnableWorkPreventsOrderedAdmissionTail) {
    anchorwave::AnchorTaskExecutor executor(4);
    std::mutex mutex;
    std::condition_variable condition;
    bool release = false;
    std::atomic<int> entered(0);
    executor.orderedHeadWaitStarted();
    for (int task = 0; task < 3; ++task) {
        executor.submit(10, [&]() {
            ++entered;
            std::unique_lock<std::mutex> lock(mutex);
            condition.wait(lock, [&]() { return release; });
        });
    }
    while (entered.load() < 3) {
        std::this_thread::yield();
    }
    const auto load = executor.loadSnapshot();
    EXPECT_EQ(3ULL, load.schedulableTasks);
    EXPECT_FALSE(load.admissionTailPhase);

    executor.orderedHeadWaitFinished();
    {
        std::lock_guard<std::mutex> lock(mutex);
        release = true;
    }
    condition.notify_all();
    executor.waitForIdle();
}

namespace {

anchorwave::AlignmentGapDescriptor gapDescriptor(
        std::size_t anchorIndex, uint64_t estimatedCost = 1) {
    anchorwave::AlignmentGapDescriptor descriptor;
    descriptor.anchorIndex = anchorIndex;
    descriptor.estimatedCost = estimatedCost;
    descriptor.schedulingPriorityCost = estimatedCost;
    return descriptor;
}

}  // namespace

TEST(ParallelGapBatchSchedulerTest, RequiredResultDoesNotWaitForWholeWindow) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::mutex mutex;
    std::condition_variable condition;
    bool secondEntered = false;
    bool releaseSecond = false;

    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    descriptors.push_back(gapDescriptor(1, 100));
    descriptors.push_back(gapDescriptor(2, 100));
    anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
            descriptors, 2, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                if (descriptor.anchorIndex == 2) {
                    std::unique_lock<std::mutex> lock(mutex);
                    secondEntered = true;
                    condition.notify_all();
                    condition.wait(lock, [&]() { return releaseSecond; });
                }
                return static_cast<int>(descriptor.anchorIndex * 10);
            });
    scheduler.start();

    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return secondEntered; });
    }
    const std::shared_ptr<const int> first = scheduler.resultBeforeAnchor(1);
    ASSERT_TRUE(first);
    EXPECT_EQ(10, *first);
    {
        std::lock_guard<std::mutex> lock(mutex);
        EXPECT_FALSE(releaseSecond);
        releaseSecond = true;
    }
    condition.notify_all();
    const std::shared_ptr<const int> second = scheduler.resultBeforeAnchor(2);
    ASSERT_TRUE(second);
    EXPECT_EQ(20, *second);
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     PrefetchesLongestGapOutsideTheOrderedFrontier) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::mutex mutex;
    std::condition_variable condition;
    bool blockerEntered = false;
    bool releaseBlocker = false;
    std::vector<std::size_t> started;

    executor.submit(std::numeric_limits<uint64_t>::max(), [&]() {
        std::unique_lock<std::mutex> lock(mutex);
        blockerEntered = true;
        condition.notify_all();
        condition.wait(lock, [&]() { return releaseBlocker; });
    });
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return blockerEntered; });
    }

    std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
            gapDescriptor(1, 1),
            gapDescriptor(2, 2),
            gapDescriptor(3, 1000)};
    anchorwave::ParallelGapBatchSchedulerCore<std::size_t> scheduler(
            descriptors, 2, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                std::lock_guard<std::mutex> lock(mutex);
                started.push_back(descriptor.anchorIndex);
                condition.notify_all();
                return descriptor.anchorIndex;
            });
    scheduler.start();
    {
        std::lock_guard<std::mutex> lock(mutex);
        releaseBlocker = true;
    }
    condition.notify_all();
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return !started.empty(); });
        ASSERT_EQ(3U, started.front());
    }

    for (std::size_t anchorIndex = 1; anchorIndex <= 3; ++anchorIndex) {
        const auto result = scheduler.resultBeforeAnchor(anchorIndex);
        ASSERT_TRUE(result);
        EXPECT_EQ(anchorIndex, *result);
        EXPECT_LE(scheduler.pendingResultCount(),
                  scheduler.maximumOutstandingResultCount());
    }
    EXPECT_EQ(3U, scheduler.peakPendingResultCount());
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     SinglePendingSlotPreservesStrictAnchorOrder) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::mutex mutex;
    std::vector<std::size_t> started;
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
            gapDescriptor(1, 1),
            gapDescriptor(2, 1000),
            gapDescriptor(3, 1000000)};
    anchorwave::ParallelGapBatchSchedulerCore<std::size_t> scheduler(
            descriptors, 1, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                std::lock_guard<std::mutex> lock(mutex);
                started.push_back(descriptor.anchorIndex);
                return descriptor.anchorIndex;
            });
    scheduler.start();
    for (std::size_t anchorIndex = 1; anchorIndex <= 3; ++anchorIndex) {
        const auto result = scheduler.resultBeforeAnchor(anchorIndex);
        ASSERT_TRUE(result);
        EXPECT_EQ(anchorIndex, *result);
        EXPECT_LE(scheduler.pendingResultCount(), 1U);
    }
    executor.waitForIdle();
    EXPECT_EQ((std::vector<std::size_t>{1, 2, 3}), started);
    EXPECT_EQ(1U, scheduler.peakPendingResultCount());
}

TEST(ParallelGapBatchSchedulerTest,
     BlockedOrderedHeadMakesDormantFutureAnAdmissionTail) {
    anchorwave::AnchorTaskExecutor executor(4);
    std::mutex mutex;
    std::condition_variable condition;
    bool firstEntered = false;
    bool releaseFirst = false;
    std::atomic<int> laterCompleted(0);

    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    for (std::size_t anchorIndex = 1; anchorIndex <= 100; ++anchorIndex) {
        descriptors.push_back(gapDescriptor(anchorIndex, 100));
    }
    anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
            descriptors, 4, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                if (descriptor.anchorIndex == 1) {
                    std::unique_lock<std::mutex> lock(mutex);
                    firstEntered = true;
                    condition.notify_all();
                    condition.wait(lock, [&]() { return releaseFirst; });
                } else {
                    ++laterCompleted;
                }
                return static_cast<int>(descriptor.anchorIndex);
            },
            nullptr,
            [](const int &) { return 11ULL; });
    scheduler.start();

    std::future<std::shared_ptr<const int>> first = std::async(
            std::launch::async, [&]() { return scheduler.resultBeforeAnchor(1); });
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return firstEntered; });
    }
    anchorwave::AnchorTaskExecutor::LoadSnapshot blocked;
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::seconds(2);
    do {
        blocked = executor.loadSnapshot();
        if (laterCompleted.load() >= 7 &&
            scheduler.pendingResultCount() == 8 &&
            scheduler.completedResultCount() == 7 &&
            blocked.globalFutureTasks == 92 &&
            blocked.admissionTailPhase) {
            break;
        }
        std::this_thread::yield();
    } while (std::chrono::steady_clock::now() < deadline);
    ASSERT_EQ(7, laterCompleted.load());
    EXPECT_EQ(8U, scheduler.pendingResultCount());
    EXPECT_EQ(8U, scheduler.submittedDescriptorCount());
    EXPECT_EQ(1U, scheduler.inFlightResultCount());
    EXPECT_EQ(7U, scheduler.completedResultCount());
    EXPECT_EQ(77ULL, scheduler.completedResultBytes());
    EXPECT_EQ(0U, scheduler.deferredResultCount());
    EXPECT_EQ(8U, scheduler.peakPendingResultCount());
    EXPECT_EQ(7U, scheduler.peakCompletedResultCount());
    EXPECT_EQ(77ULL, scheduler.peakCompletedResultBytes());
    EXPECT_EQ(92ULL, blocked.globalFutureTasks);
    EXPECT_FALSE(blocked.globalTailPhase);
    EXPECT_TRUE(blocked.admissionTailPhase);
    EXPECT_EQ(1ULL, blocked.blockedOrderedHeads);
    EXPECT_LE(blocked.schedulableTasks, 2ULL);

    {
        std::lock_guard<std::mutex> lock(mutex);
        releaseFirst = true;
    }
    condition.notify_all();
    const std::shared_ptr<const int> result = first.get();
    ASSERT_TRUE(result);
    EXPECT_EQ(1, *result);
    EXPECT_EQ(0ULL, executor.loadSnapshot().blockedOrderedHeads);
}

TEST(ParallelGapBatchSchedulerTest,
     CompletionBackfillStaysWithinTheOutstandingResultCeiling) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::atomic<uint64_t> submitted(0);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    for (std::size_t anchorIndex = 1; anchorIndex <= 5; ++anchorIndex) {
        descriptors.push_back(gapDescriptor(anchorIndex, anchorIndex));
    }
    anchorwave::ParallelGapBatchSchedulerCore<std::size_t> scheduler(
            descriptors, 2, executor,
            [](const anchorwave::AlignmentGapDescriptor &descriptor) {
                return descriptor.anchorIndex;
            },
            &submitted);
    scheduler.start();

    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::seconds(2);
    while (scheduler.submittedDescriptorCount() < 4 &&
           std::chrono::steady_clock::now() < deadline) {
        std::this_thread::yield();
    }
    EXPECT_EQ(4U, scheduler.maximumOutstandingResultCount());
    EXPECT_EQ(4U, scheduler.pendingResultCount());
    EXPECT_EQ(4U, scheduler.submittedDescriptorCount());
    EXPECT_EQ(4ULL, submitted.load());
    {
        const auto load = executor.loadSnapshot();
        EXPECT_EQ(1ULL, load.plannedTasks);
        EXPECT_GE(load.outstandingTasks, 1ULL);
    }
    for (std::size_t anchorIndex = 1; anchorIndex <= 5; ++anchorIndex) {
        const std::shared_ptr<const std::size_t> result =
                scheduler.resultBeforeAnchor(anchorIndex);
        ASSERT_TRUE(result);
        EXPECT_EQ(anchorIndex, *result);
        EXPECT_LE(scheduler.pendingResultCount(), 4U);
    }
    EXPECT_EQ(5U, scheduler.submittedDescriptorCount());
    EXPECT_EQ(5ULL, submitted.load());
    EXPECT_EQ(4U, scheduler.peakPendingResultCount());
    EXPECT_EQ(0U, scheduler.pendingResultCount());
    EXPECT_EQ(0U, scheduler.inFlightResultCount());
    EXPECT_EQ(0U, scheduler.completedResultCount());
    EXPECT_EQ(0ULL, scheduler.completedResultBytes());
    EXPECT_EQ(0U, scheduler.deferredResultCount());
    executor.waitForIdle();
    EXPECT_LE(executor.peakActiveTasks(), 2);
}

TEST(ParallelGapBatchSchedulerTest,
     IdleWorkersPermitOneBoundedEmergencyBackfillBehindDeferredHead) {
    anchorwave::AnchorTaskExecutor executor(4);
    std::mutex mutex;
    std::condition_variable condition;
    bool firstEntered = false;
    bool releaseBlockedTasks = false;
    std::atomic<int> ninthCompleted(0);

    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    for (std::size_t anchorIndex = 1; anchorIndex <= 20; ++anchorIndex) {
        descriptors.push_back(gapDescriptor(anchorIndex, 100));
    }
    anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
            descriptors, 4, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                if (descriptor.anchorIndex == 1) {
                    std::unique_lock<std::mutex> lock(mutex);
                    firstEntered = true;
                    condition.notify_all();
                    condition.wait(lock, [&]() { return releaseBlockedTasks; });
                } else if (descriptor.anchorIndex <= 8) {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (!releaseBlockedTasks) {
                        throw anchorwave::AlignmentTaskDeferred(
                                std::chrono::milliseconds(20));
                    }
                } else if (descriptor.anchorIndex == 9) {
                    ++ninthCompleted;
                }
                return static_cast<int>(descriptor.anchorIndex);
            });
    scheduler.start();

    std::future<std::shared_ptr<const int>> first = std::async(
            std::launch::async,
            [&]() { return scheduler.resultBeforeAnchor(1); });
    {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [&]() { return firstEntered; });
    }
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::seconds(3);
    while ((scheduler.emergencyBackfillCount() == 0 ||
            ninthCompleted.load() == 0) &&
           std::chrono::steady_clock::now() < deadline) {
        std::this_thread::yield();
    }
    EXPECT_EQ(9U, scheduler.maximumEmergencyOutstandingResultCount());
    EXPECT_EQ(1U, scheduler.emergencyBackfillCount());
    EXPECT_EQ(9U, scheduler.peakPendingResultCount());
    EXPECT_EQ(1, ninthCompleted.load());

    {
        std::lock_guard<std::mutex> lock(mutex);
        releaseBlockedTasks = true;
    }
    condition.notify_all();
    ASSERT_TRUE(first.get());
    for (std::size_t anchorIndex = 2; anchorIndex <= 20; ++anchorIndex) {
        const auto result = scheduler.resultBeforeAnchor(anchorIndex);
        ASSERT_TRUE(result);
        EXPECT_EQ(static_cast<int>(anchorIndex), *result);
        EXPECT_LE(scheduler.pendingResultCount(), 9U);
    }
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     DestructionCancelsDescriptorsOutsideTheRollingWindow) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    for (std::size_t anchorIndex = 1; anchorIndex <= 6; ++anchorIndex) {
        descriptors.push_back(gapDescriptor(anchorIndex, 10));
    }
    {
        anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
                descriptors, 1, executor,
                [](const anchorwave::AlignmentGapDescriptor &descriptor) {
                    return static_cast<int>(descriptor.anchorIndex);
                });
        scheduler.start();
        EXPECT_EQ(5ULL, executor.loadSnapshot().plannedTasks);
    }
    EXPECT_NO_THROW(executor.waitForIdle());
    EXPECT_EQ(0ULL, executor.loadSnapshot().plannedTasks);
}

TEST(ParallelGapBatchSchedulerTest, PreservesDeferralUntilResultIsReady) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::atomic<int> attempts(0);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
            gapDescriptor(7, 10)};
    anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
            descriptors, 1, executor,
            [&](const anchorwave::AlignmentGapDescriptor &) {
                if (++attempts == 1) {
                    throw anchorwave::AlignmentTaskDeferred(
                            std::chrono::milliseconds(1));
                }
                return 73;
            });
    scheduler.start();
    const std::shared_ptr<const int> result = scheduler.resultBeforeAnchor(7);
    ASSERT_TRUE(result);
    EXPECT_EQ(73, *result);
    EXPECT_EQ(2, attempts.load());
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     MemoryDeferralBackfillsWithoutDuplicatingTheDeferredDescriptor) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::atomic<uint64_t> submitted(0);
    std::mutex attemptsMutex;
    std::vector<int> attempts(7, 0);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    for (std::size_t anchorIndex = 1; anchorIndex <= 6; ++anchorIndex) {
        descriptors.push_back(gapDescriptor(anchorIndex, 10));
    }

    anchorwave::ParallelGapBatchSchedulerCore<std::size_t> scheduler(
            descriptors, 2, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                int attempt = 0;
                {
                    std::lock_guard<std::mutex> lock(attemptsMutex);
                    attempt = ++attempts[descriptor.anchorIndex];
                }
                if (descriptor.anchorIndex == 1 && attempt == 1) {
                    throw anchorwave::AlignmentTaskDeferred(
                            std::chrono::milliseconds(400));
                }
                return descriptor.anchorIndex;
            },
            &submitted,
            [](const std::size_t &) { return 13ULL; });
    scheduler.start();

    const auto backfillDeadline = std::chrono::steady_clock::now() +
                                  std::chrono::milliseconds(250);
    while ((scheduler.peakDeferredResultCount() == 0 ||
            scheduler.submittedDescriptorCount() < 4 ||
            scheduler.completedResultCount() < 3) &&
           std::chrono::steady_clock::now() < backfillDeadline) {
        std::this_thread::yield();
    }
    EXPECT_EQ(1U, scheduler.peakDeferredResultCount());
    EXPECT_EQ(1U, scheduler.deferredResultCount());
    const std::size_t pendingBeforeConsumption =
            scheduler.pendingResultCount();
    EXPECT_GE(pendingBeforeConsumption, 4U);
    EXPECT_LE(pendingBeforeConsumption,
              scheduler.maximumEmergencyOutstandingResultCount());
    EXPECT_EQ(pendingBeforeConsumption,
              scheduler.submittedDescriptorCount());
    EXPECT_EQ(pendingBeforeConsumption, submitted.load());
    EXPECT_EQ(pendingBeforeConsumption - 1,
              scheduler.completedResultCount());
    EXPECT_EQ((pendingBeforeConsumption - 1) * 13ULL,
              scheduler.completedResultBytes());
    EXPECT_EQ(pendingBeforeConsumption - 4,
              scheduler.emergencyBackfillCount());

    std::vector<std::size_t> emitted;
    for (std::size_t anchorIndex = 1; anchorIndex <= 6; ++anchorIndex) {
        const auto result = scheduler.resultBeforeAnchor(anchorIndex);
        ASSERT_TRUE(result);
        emitted.push_back(*result);
        EXPECT_LE(scheduler.pendingResultCount(),
                  scheduler.maximumEmergencyOutstandingResultCount());
    }
    EXPECT_EQ((std::vector<std::size_t>{1, 2, 3, 4, 5, 6}), emitted);
    {
        std::lock_guard<std::mutex> lock(attemptsMutex);
        EXPECT_EQ(2, attempts[1]);
        for (std::size_t anchorIndex = 2; anchorIndex <= 6; ++anchorIndex) {
            EXPECT_EQ(1, attempts[anchorIndex]);
        }
    }
    EXPECT_EQ(6U, scheduler.submittedDescriptorCount());
    EXPECT_EQ(6ULL, submitted.load());
    EXPECT_EQ(0U, scheduler.pendingResultCount());
    EXPECT_EQ(0U, scheduler.inFlightResultCount());
    EXPECT_EQ(0U, scheduler.completedResultCount());
    EXPECT_EQ(0ULL, scheduler.completedResultBytes());
    EXPECT_EQ(0U, scheduler.deferredResultCount());
    EXPECT_LE(scheduler.peakPendingResultCount(),
              scheduler.maximumEmergencyOutstandingResultCount());
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     StressPreservesOrderedExactlyOnceResultsAcrossMixedDeferrals) {
    constexpr std::size_t descriptorCount = 200;
    constexpr uint64_t resultBytes = 17;
    anchorwave::AnchorTaskExecutor executor(4);
    std::mutex attemptsMutex;
    std::vector<int> attempts(descriptorCount + 1, 0);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors;
    descriptors.reserve(descriptorCount);
    for (std::size_t anchorIndex = 1;
         anchorIndex <= descriptorCount; ++anchorIndex) {
        descriptors.push_back(gapDescriptor(
                anchorIndex, (anchorIndex * 37) % 101 + 1));
    }

    anchorwave::ParallelGapBatchSchedulerCore<std::size_t> scheduler(
            descriptors, 4, executor,
            [&](const anchorwave::AlignmentGapDescriptor &descriptor) {
                int attempt = 0;
                {
                    std::lock_guard<std::mutex> lock(attemptsMutex);
                    attempt = ++attempts[descriptor.anchorIndex];
                }
                if (descriptor.anchorIndex % 9 == 0 && attempt == 1) {
                    throw anchorwave::AlignmentTaskDeferred(
                            std::chrono::milliseconds(1));
                }
                return descriptor.anchorIndex;
            },
            nullptr,
            [](const std::size_t &) { return resultBytes; });
    scheduler.start();

    for (std::size_t anchorIndex = 1;
         anchorIndex <= descriptorCount; ++anchorIndex) {
        const auto result = scheduler.resultBeforeAnchor(anchorIndex);
        ASSERT_TRUE(result);
        EXPECT_EQ(anchorIndex, *result);
        EXPECT_LE(scheduler.pendingResultCount(),
                  scheduler.maximumEmergencyOutstandingResultCount());
    }
    executor.waitForIdle();

    {
        std::lock_guard<std::mutex> lock(attemptsMutex);
        for (std::size_t anchorIndex = 1;
             anchorIndex <= descriptorCount; ++anchorIndex) {
            EXPECT_EQ(anchorIndex % 9 == 0 ? 2 : 1,
                      attempts[anchorIndex]);
        }
    }
    EXPECT_EQ(descriptorCount, scheduler.submittedDescriptorCount());
    EXPECT_EQ(0U, scheduler.pendingResultCount());
    EXPECT_EQ(0U, scheduler.inFlightResultCount());
    EXPECT_EQ(0U, scheduler.completedResultCount());
    EXPECT_EQ(0ULL, scheduler.completedResultBytes());
    EXPECT_EQ(0U, scheduler.deferredResultCount());
    EXPECT_LE(scheduler.peakPendingResultCount(),
              scheduler.maximumEmergencyOutstandingResultCount());
    EXPECT_LE(scheduler.peakInFlightResultCount(),
              scheduler.maximumEmergencyOutstandingResultCount());
    EXPECT_LE(scheduler.peakCompletedResultCount(),
              scheduler.maximumEmergencyOutstandingResultCount());
    EXPECT_LE(scheduler.peakCompletedResultBytes(),
              scheduler.maximumEmergencyOutstandingResultCount() *
                      resultBytes);
    EXPECT_LE(scheduler.peakDeferredResultCount(),
              scheduler.maximumEmergencyOutstandingResultCount());
}

TEST(ParallelGapBatchSchedulerTest,
     RejectsDescriptorsOutsideStrictAnchorOrder) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
            gapDescriptor(2, 10), gapDescriptor(1, 20)};
    anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
            descriptors, 1, executor,
            [](const anchorwave::AlignmentGapDescriptor &descriptor) {
                return static_cast<int>(descriptor.anchorIndex);
            });
    EXPECT_THROW(scheduler.start(), std::invalid_argument);
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     DeferredTaskRetainsItsPrivatePreparedSequenceCache) {
    anchorwave::AnchorTaskExecutor executor(1);
    std::atomic<int> attempts(0);
    std::atomic<bool> cacheWasRetained(false);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
            gapDescriptor(1, 100)};
    anchorwave::ParallelGapBatchSchedulerCore<std::size_t> scheduler(
            descriptors, 1, executor,
            [&](anchorwave::AlignmentGapDescriptor &descriptor) {
                const int attempt = ++attempts;
                if (!descriptor.preparedReferenceSequence) {
                    descriptor.preparedReferenceSequence =
                            std::make_shared<std::string>("ACGT");
                } else {
                    cacheWasRetained.store(
                            *descriptor.preparedReferenceSequence == "ACGT");
                }
                if (attempt == 1) {
                    throw anchorwave::AlignmentTaskDeferred(
                            std::chrono::milliseconds(1));
                }
                return descriptor.preparedReferenceSequence->size();
            });
    scheduler.start();
    const auto result = scheduler.resultBeforeAnchor(1);
    ASSERT_TRUE(result);
    EXPECT_EQ(4U, *result);
    EXPECT_EQ(2, attempts.load());
    EXPECT_TRUE(cacheWasRetained.load());
    executor.waitForIdle();
}

TEST(ParallelGapBatchSchedulerTest,
     CooperativelyCompletesWhenEveryWorkerOwnsAWindow) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::atomic<int> completedResults(0);

    for (int parent = 0; parent < 2; ++parent) {
        executor.submit(100, [&]() {
            std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
                    gapDescriptor(1, 10), gapDescriptor(2, 10)};
            anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
                    descriptors, 2, executor,
                    [](const anchorwave::AlignmentGapDescriptor &descriptor) {
                        return static_cast<int>(descriptor.anchorIndex);
                    });
            scheduler.start();
            for (std::size_t anchorIndex = 1; anchorIndex <= 2;
                 ++anchorIndex) {
                const std::shared_ptr<const int> result =
                        scheduler.resultBeforeAnchor(anchorIndex);
                ASSERT_TRUE(result);
                EXPECT_EQ(static_cast<int>(anchorIndex), *result);
                ++completedResults;
            }
        });
    }

    executor.waitForIdle();
    EXPECT_EQ(4, completedResults.load());
    EXPECT_LE(executor.peakActiveTasks(), 2);
}

TEST(ParallelGapBatchSchedulerTest, PropagatesIndividualTaskException) {
    anchorwave::AnchorTaskExecutor executor(2);
    std::vector<anchorwave::AlignmentGapDescriptor> descriptors{
            gapDescriptor(1, 100), gapDescriptor(2, 1)};
    {
        anchorwave::ParallelGapBatchSchedulerCore<int> scheduler(
                descriptors, 2, executor,
                [](const anchorwave::AlignmentGapDescriptor &descriptor) {
                    if (descriptor.anchorIndex == 1) {
                        throw std::runtime_error("gap failed");
                    }
                    return 2;
                });
        scheduler.start();
        EXPECT_THROW(scheduler.resultBeforeAnchor(1), std::runtime_error);
    }
    EXPECT_THROW(executor.waitForIdle(), std::runtime_error);
}
