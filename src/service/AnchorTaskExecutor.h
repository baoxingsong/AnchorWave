#pragma once

#include <condition_variable>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <memory>
#include <queue>
#include <set>
#include <thread>
#include <vector>

namespace anchorwave {

// A deferrable alignment task throws this when its preferred high-WFA memory
// is temporarily unavailable. The executor parks the task until retryAfter,
// releasing the worker to run other queued work in the meantime.
class AlignmentTaskDeferred : public std::exception {
public:
    explicit AlignmentTaskDeferred(
            std::chrono::milliseconds retryDelay =
                    std::chrono::milliseconds(250));
    const char *what() const noexcept override;
    std::chrono::milliseconds retryDelay() const;

private:
    std::chrono::milliseconds retryDelay_;
};

class AnchorTaskGroup {
public:
    AnchorTaskGroup() = default;

    AnchorTaskGroup(const AnchorTaskGroup &) = delete;
    AnchorTaskGroup &operator=(const AnchorTaskGroup &) = delete;

private:
    friend class AnchorTaskExecutor;

    std::atomic<std::size_t> remaining_{0};
    std::mutex mutex_;
    std::condition_variable completed_;
    std::exception_ptr firstException_;
};

// Fixed-size executor used by the anchor-organization stage. Workers pull the
// highest-cost available chromosome task from one shared queue, which avoids
// binding a slow chromosome to a particular worker. Callers keep task output
// private and merge it in deterministic chromosome order after waitForIdle().
class AnchorTaskExecutor {
public:
    explicit AnchorTaskExecutor(int maximumThreads);
    ~AnchorTaskExecutor();

    AnchorTaskExecutor(const AnchorTaskExecutor &) = delete;
    AnchorTaskExecutor &operator=(const AnchorTaskExecutor &) = delete;

    void submit(uint64_t estimatedCost, std::function<void()> function);
    void submit(AnchorTaskGroup &group, uint64_t estimatedCost,
                std::function<void()> function);
    void submitDeferrable(AnchorTaskGroup &group, uint64_t estimatedCost,
                          std::function<void()> function);
    // A rolling result window cannot submit its entire chromosome backlog at
    // once without retaining many completed alignment strings. Register the
    // not-yet-submitted costs so tail/makespan decisions still see the global
    // remaining work. submitPlannedDeferrable atomically transfers one cost
    // from that future backlog into the runnable queue.
    void registerPlannedTaskCosts(const std::vector<uint64_t> &estimatedCosts);
    void cancelPlannedTaskCosts(
            const std::vector<uint64_t> &estimatedCosts) noexcept;
    void submitPlannedDeferrable(
            AnchorTaskGroup &group, uint64_t estimatedCost,
            std::function<void()> function);
    void waitForIdle();
    void waitForGroup(AnchorTaskGroup &group);
    void waitUntilGroupSizeAtMost(AnchorTaskGroup &group,
                                  std::size_t maximumRemaining);

    // Ordered consumers (for example a chromosome MAF writer) can wait on
    // one gap while later results already occupy the rolling window.  Keep
    // this frontier state separate from the global future backlog so the
    // alignment selector can recognize an effective admission tail.
    void orderedHeadWaitStarted();
    void orderedHeadWaitFinished() noexcept;

    int workerCount() const;
    int peakActiveTasks() const;
    bool isTailPhase() const;

    struct LoadSnapshot {
        int workerCount = 1;
        int activeTasks = 0;
        uint64_t readyTasks = 0;
        uint64_t deferredTasks = 0;
        uint64_t plannedTasks = 0;
        uint64_t schedulableTasks = 0;
        uint64_t globalFutureTasks = 0;
        uint64_t blockedOrderedHeads = 0;
        uint64_t outstandingTasks = 0;
        uint64_t outstandingEstimatedCost = 0;
        uint64_t criticalEstimatedCost = 0;
        bool globalTailPhase = false;
        bool admissionTailPhase = false;
        bool tailPhase = false;
    };
    LoadSnapshot loadSnapshot() const;

private:
    struct Task {
        uint64_t estimatedCost = 0;
        uint64_t accountedCost = 0;
        uint64_t sequence = 0;
        AnchorTaskGroup *group = nullptr;
        std::function<void()> function;
        bool deferrable = false;
        uint64_t deferralCount = 0;
        std::chrono::steady_clock::time_point firstDeferred;
        std::chrono::steady_clock::time_point retryAfter;
        // Survives a deferral so expensive sequence profiling does not have to
        // be repeated at every retry.
        std::shared_ptr<void> retainedState;
    };

    struct TaskLess {
        bool operator()(const Task &left, const Task &right) const;
    };

    void workerLoop();
    void executeTask(Task task, bool countAsActiveTask);
    void recordException(AnchorTaskGroup *group,
                         std::exception_ptr exception);
    void submitToGroup(AnchorTaskGroup &group, uint64_t estimatedCost,
                       std::function<void()> function, bool deferrable,
                       bool consumePlannedCost = false);
    void promoteReadyDeferredLocked(
            std::chrono::steady_clock::time_point now);
    std::chrono::steady_clock::time_point nextDeferredRetryLocked() const;

    const int workerCount_;
    mutable std::mutex mutex_;
    std::condition_variable taskAvailable_;
    std::condition_variable idle_;
    std::priority_queue<Task, std::vector<Task>, TaskLess> tasks_;
    std::vector<Task> deferredTasks_;
    std::vector<std::thread> workers_;
    uint64_t nextSequence_ = 0;
    uint64_t outstandingEstimatedCost_ = 0;
    std::multiset<uint64_t> outstandingCosts_;
    uint64_t plannedEstimatedCost_ = 0;
    uint64_t blockedOrderedHeads_ = 0;
    std::multiset<uint64_t> plannedCosts_;
    int activeTasks_ = 0;
    int peakActiveTasks_ = 0;
    bool stopping_ = false;
    std::exception_ptr firstException_;
};

// Context queries used by alignSlidingWindow(). Only explicitly deferrable
// fine-grained gap tasks return true; chromosome parent tasks retain the
// blocking compatibility path.
bool currentAnchorTaskCanDefer();
bool currentAnchorTaskShouldDrainMemory();
uint64_t currentAnchorTaskDeferralCount();
AnchorTaskExecutor::LoadSnapshot currentAnchorTaskLoadSnapshot();

// One opaque state slot associated with the current deferrable task.  It is
// retained while the task is parked and destroyed when the task completes.
std::shared_ptr<void> currentAnchorTaskRetainedState();
void setCurrentAnchorTaskRetainedState(std::shared_ptr<void> state);

// The current longest-path implementation is quadratic in anchor count, so
// n^2 is a useful first-order scheduling weight. Saturate on overflow.
uint64_t anchorTaskEstimatedCost(std::size_t anchorCount);
uint64_t anchorGapTaskEstimatedCost(std::size_t referenceLength,
                                    std::size_t queryLength);

}  // namespace anchorwave
