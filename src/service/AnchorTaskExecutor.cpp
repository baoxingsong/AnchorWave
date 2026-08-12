#include "AnchorTaskExecutor.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace anchorwave {

namespace {
thread_local AnchorTaskExecutor *currentAnchorExecutor = nullptr;
thread_local bool currentTaskDeferrable = false;
thread_local uint64_t currentTaskDeferralCount = 0;
thread_local std::chrono::steady_clock::time_point currentTaskFirstDeferred;
thread_local std::shared_ptr<void> *currentTaskRetainedState = nullptr;
constexpr auto kMemoryStarvationAge = std::chrono::seconds(30);

uint64_t saturatingAdd(uint64_t first, uint64_t second) {
    return second > std::numeric_limits<uint64_t>::max() - first
           ? std::numeric_limits<uint64_t>::max() : first + second;
}
}

AlignmentTaskDeferred::AlignmentTaskDeferred(
        std::chrono::milliseconds retryDelay)
        : retryDelay_(std::max(retryDelay, std::chrono::milliseconds(1))) {}

const char *AlignmentTaskDeferred::what() const noexcept {
    return "alignment task deferred for preferred memory";
}

std::chrono::milliseconds AlignmentTaskDeferred::retryDelay() const {
    return retryDelay_;
}

bool AnchorTaskExecutor::TaskLess::operator()(
        const Task &left, const Task &right) const {
    if (left.estimatedCost != right.estimatedCost) {
        return left.estimatedCost < right.estimatedCost;
    }
    // std::priority_queue returns the "largest" element. Prefer the task that
    // was submitted first when estimated costs are equal.
    return left.sequence > right.sequence;
}

AnchorTaskExecutor::AnchorTaskExecutor(int maximumThreads)
        : workerCount_(maximumThreads) {
    if (maximumThreads <= 0) {
        throw std::invalid_argument(
                "anchor task executor thread count must be positive");
    }
    workers_.reserve(static_cast<std::size_t>(workerCount_));
    for (int worker = 0; worker < workerCount_; ++worker) {
        workers_.emplace_back(&AnchorTaskExecutor::workerLoop, this);
    }
}

AnchorTaskExecutor::~AnchorTaskExecutor() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stopping_ = true;
    }
    taskAvailable_.notify_all();
    for (std::thread &worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

void AnchorTaskExecutor::submit(
        uint64_t estimatedCost, std::function<void()> function) {
    if (!function) {
        throw std::invalid_argument("anchor task must not be empty");
    }
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (stopping_) {
            throw std::runtime_error("cannot submit to a stopped anchor executor");
        }
        Task task;
        task.estimatedCost = estimatedCost;
        task.accountedCost = estimatedCost;
        task.sequence = nextSequence_++;
        task.function = std::move(function);
        tasks_.push(std::move(task));
        outstandingEstimatedCost_ = saturatingAdd(
                outstandingEstimatedCost_, estimatedCost);
        outstandingCosts_.insert(estimatedCost);
    }
    taskAvailable_.notify_one();
}

void AnchorTaskExecutor::submit(
        AnchorTaskGroup &group, uint64_t estimatedCost,
        std::function<void()> function) {
    submitToGroup(group, estimatedCost, std::move(function), false, false);
}

void AnchorTaskExecutor::submitDeferrable(
        AnchorTaskGroup &group, uint64_t estimatedCost,
        std::function<void()> function) {
    submitToGroup(group, estimatedCost, std::move(function), true, false);
}

void AnchorTaskExecutor::registerPlannedTaskCosts(
        const std::vector<uint64_t> &estimatedCosts) {
    if (estimatedCosts.empty()) {
        return;
    }
    std::lock_guard<std::mutex> lock(mutex_);
    if (stopping_) {
        throw std::runtime_error(
                "cannot register work with a stopped anchor executor");
    }
    for (const uint64_t cost : estimatedCosts) {
        plannedCosts_.insert(cost);
        plannedEstimatedCost_ = saturatingAdd(plannedEstimatedCost_, cost);
    }
}

void AnchorTaskExecutor::cancelPlannedTaskCosts(
        const std::vector<uint64_t> &estimatedCosts) noexcept {
    if (estimatedCosts.empty()) {
        return;
    }
    std::lock_guard<std::mutex> lock(mutex_);
    for (const uint64_t cost : estimatedCosts) {
        const auto found = plannedCosts_.find(cost);
        if (found == plannedCosts_.end()) {
            continue;
        }
        plannedCosts_.erase(found);
        plannedEstimatedCost_ = cost > plannedEstimatedCost_
                                ? 0 : plannedEstimatedCost_ - cost;
    }
    if (tasks_.empty() && deferredTasks_.empty() && activeTasks_ == 0 &&
        plannedCosts_.empty()) {
        idle_.notify_all();
    }
}

void AnchorTaskExecutor::submitPlannedDeferrable(
        AnchorTaskGroup &group, uint64_t estimatedCost,
        std::function<void()> function) {
    submitToGroup(group, estimatedCost, std::move(function), true, true);
}

void AnchorTaskExecutor::submitToGroup(
        AnchorTaskGroup &group, uint64_t estimatedCost,
        std::function<void()> function, bool deferrable,
        bool consumePlannedCost) {
    if (!function) {
        throw std::invalid_argument("anchor task must not be empty");
    }
    group.remaining_.fetch_add(1, std::memory_order_relaxed);
    try {
        std::lock_guard<std::mutex> lock(mutex_);
        if (stopping_) {
            throw std::runtime_error(
                    "cannot submit to a stopped anchor executor");
        }
        const auto planned = consumePlannedCost
                             ? plannedCosts_.find(estimatedCost)
                             : plannedCosts_.end();
        if (consumePlannedCost && planned == plannedCosts_.end()) {
            throw std::logic_error(
                    "submitted rolling task was not registered as planned");
        }
        Task task;
        task.estimatedCost = estimatedCost;
        task.accountedCost = estimatedCost;
        task.sequence = nextSequence_++;
        task.group = &group;
        task.function = std::move(function);
        task.deferrable = deferrable;
        tasks_.push(std::move(task));
        outstandingEstimatedCost_ = saturatingAdd(
                outstandingEstimatedCost_, estimatedCost);
        outstandingCosts_.insert(estimatedCost);
        if (consumePlannedCost) {
            plannedCosts_.erase(planned);
            plannedEstimatedCost_ = estimatedCost > plannedEstimatedCost_
                                    ? 0
                                    : plannedEstimatedCost_ - estimatedCost;
        }
    } catch (...) {
        group.remaining_.fetch_sub(1, std::memory_order_relaxed);
        throw;
    }
    taskAvailable_.notify_one();
}

void AnchorTaskExecutor::waitForIdle() {
    std::exception_ptr exception;
    {
        std::unique_lock<std::mutex> lock(mutex_);
        idle_.wait(lock, [this]() {
            return tasks_.empty() && deferredTasks_.empty() &&
                   activeTasks_ == 0 && plannedCosts_.empty();
        });
        exception = firstException_;
        firstException_ = nullptr;
    }
    if (exception) {
        std::rethrow_exception(exception);
    }
}

void AnchorTaskExecutor::waitForGroup(AnchorTaskGroup &group) {
    waitUntilGroupSizeAtMost(group, 0);

    std::exception_ptr exception;
    {
        std::lock_guard<std::mutex> lock(group.mutex_);
        exception = group.firstException_;
    }
    if (exception) {
        std::rethrow_exception(exception);
    }
}

void AnchorTaskExecutor::waitUntilGroupSizeAtMost(
        AnchorTaskGroup &group, std::size_t maximumRemaining) {
    if (currentAnchorExecutor != this) {
        std::unique_lock<std::mutex> groupLock(group.mutex_);
        group.completed_.wait(groupLock, [&group, maximumRemaining]() {
            return group.remaining_.load(std::memory_order_acquire) <=
                   maximumRemaining;
        });
    } else {
        // A chromosome task may create independent gap tasks. Waiting without
        // helping would deadlock when every worker is occupied by a chromosome
        // parent. The waiting worker therefore executes queued work itself;
        // this adds no operating-system threads and preserves the global cap.
        while (group.remaining_.load(std::memory_order_acquire) >
               maximumRemaining) {
            Task task;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                while (true) {
                    promoteReadyDeferredLocked(
                            std::chrono::steady_clock::now());
                    if (stopping_ || !tasks_.empty() ||
                        group.remaining_.load(std::memory_order_acquire) <=
                                maximumRemaining) {
                        break;
                    }
                    if (deferredTasks_.empty()) {
                        taskAvailable_.wait(lock);
                    } else {
                        taskAvailable_.wait_until(
                                lock, nextDeferredRetryLocked());
                    }
                }
                if (group.remaining_.load(std::memory_order_acquire) <=
                    maximumRemaining) {
                    break;
                }
                if (tasks_.empty()) {
                    continue;
                }
                task = tasks_.top();
                tasks_.pop();
            }
            executeTask(std::move(task), false);
        }
    }
}

int AnchorTaskExecutor::workerCount() const {
    return workerCount_;
}

int AnchorTaskExecutor::peakActiveTasks() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return peakActiveTasks_;
}

bool AnchorTaskExecutor::isTailPhase() const {
    return loadSnapshot().tailPhase;
}

AnchorTaskExecutor::LoadSnapshot AnchorTaskExecutor::loadSnapshot() const {
    std::lock_guard<std::mutex> lock(mutex_);
    LoadSnapshot snapshot;
    snapshot.workerCount = workerCount_;
    snapshot.activeTasks = activeTasks_;
    snapshot.readyTasks = tasks_.size();
    snapshot.deferredTasks = deferredTasks_.size();
    snapshot.plannedTasks = plannedCosts_.size();
    snapshot.schedulableTasks = tasks_.size() + deferredTasks_.size() +
                                static_cast<std::size_t>(activeTasks_);
    snapshot.globalFutureTasks = plannedCosts_.size();
    snapshot.blockedOrderedHeads = blockedOrderedHeads_;
    const std::size_t outstanding = tasks_.size() + deferredTasks_.size() +
                                    static_cast<std::size_t>(activeTasks_) +
                                    plannedCosts_.size();
    snapshot.outstandingTasks = outstanding;
    snapshot.outstandingEstimatedCost = saturatingAdd(
            outstandingEstimatedCost_, plannedEstimatedCost_);
    const uint64_t activeCritical = outstandingCosts_.empty()
                                    ? 0 : *outstandingCosts_.rbegin();
    const uint64_t plannedCritical = plannedCosts_.empty()
                                     ? 0 : *plannedCosts_.rbegin();
    snapshot.criticalEstimatedCost = std::max(activeCritical,
                                              plannedCritical);
    // A worker-sized initial batch must not be called a tail on a 96/128-core
    // node. Enter draining only when few jobs remain, or when a dominant
    // critical-path job accounts for most remaining work and the ready set is
    // already smaller than half the machine.
    const std::size_t countTailThreshold = std::max<std::size_t>(
            2, static_cast<std::size_t>(workerCount_) / 8);
    const bool countTail = outstanding <= countTailThreshold;
    const uint64_t critical = snapshot.criticalEstimatedCost;
    const bool criticalTail = critical > 0 &&
            outstanding <= std::max<std::size_t>(
                    2, static_cast<std::size_t>(workerCount_) / 2) &&
            snapshot.outstandingEstimatedCost <= saturatingAdd(
                    critical, critical / 2);
    snapshot.globalTailPhase = countTail || criticalTail;
    // A large global future list is not runnable while completed results fill
    // an ordered rolling window behind its first missing result.  If every
    // currently schedulable task can be explained by an ordered head plus its
    // waiting parent, the machine is in an admission tail even though the
    // chromosome still contains many future descriptors.  This is the point
    // where a fast, memory-hungry exact engine should not be penalized by the
    // dormant future backlog.
    const uint64_t orderedFrontierAllowance = saturatingAdd(
            snapshot.blockedOrderedHeads, 1);
    const bool orderedAdmissionTail = snapshot.blockedOrderedHeads > 0 &&
            snapshot.schedulableTasks <= orderedFrontierAllowance;
    snapshot.admissionTailPhase = snapshot.globalTailPhase ||
                                  orderedAdmissionTail;
    // Preserve the historical field as the global completion-tail signal.
    snapshot.tailPhase = snapshot.globalTailPhase;
    return snapshot;
}

void AnchorTaskExecutor::orderedHeadWaitStarted() {
    std::lock_guard<std::mutex> lock(mutex_);
    ++blockedOrderedHeads_;
}

void AnchorTaskExecutor::orderedHeadWaitFinished() noexcept {
    std::lock_guard<std::mutex> lock(mutex_);
    if (blockedOrderedHeads_ > 0) {
        --blockedOrderedHeads_;
    }
}

void AnchorTaskExecutor::promoteReadyDeferredLocked(
        std::chrono::steady_clock::time_point now) {
    auto item = deferredTasks_.begin();
    while (item != deferredTasks_.end()) {
        if (item->retryAfter <= now || stopping_) {
            tasks_.push(std::move(*item));
            item = deferredTasks_.erase(item);
        } else {
            ++item;
        }
    }
}

std::chrono::steady_clock::time_point
AnchorTaskExecutor::nextDeferredRetryLocked() const {
    auto next = std::chrono::steady_clock::time_point::max();
    for (const Task &task : deferredTasks_) {
        next = std::min(next, task.retryAfter);
    }
    return next;
}

void AnchorTaskExecutor::workerLoop() {
    currentAnchorExecutor = this;
    while (true) {
        Task task;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            while (true) {
                promoteReadyDeferredLocked(std::chrono::steady_clock::now());
                if (stopping_ || !tasks_.empty()) {
                    break;
                }
                if (deferredTasks_.empty()) {
                    taskAvailable_.wait(lock);
                } else {
                    taskAvailable_.wait_until(
                            lock, nextDeferredRetryLocked());
                }
            }
            if (stopping_ && tasks_.empty() && deferredTasks_.empty()) {
                currentAnchorExecutor = nullptr;
                return;
            }
            task = tasks_.top();
            tasks_.pop();
            ++activeTasks_;
            peakActiveTasks_ = std::max(peakActiveTasks_, activeTasks_);
        }

        executeTask(std::move(task), true);
    }
}

void AnchorTaskExecutor::executeTask(Task task, bool countAsActiveTask) {
    const bool previousDeferrable = currentTaskDeferrable;
    const uint64_t previousDeferralCount = currentTaskDeferralCount;
    const auto previousFirstDeferred = currentTaskFirstDeferred;
    std::shared_ptr<void> *const previousRetainedState =
            currentTaskRetainedState;
    currentTaskDeferrable = task.deferrable;
    currentTaskDeferralCount = task.deferralCount;
    currentTaskFirstDeferred = task.firstDeferred;
    currentTaskRetainedState = &task.retainedState;
    bool deferred = false;
    try {
        task.function();
    } catch (const AlignmentTaskDeferred &request) {
        bool canDefer = false;
        if (task.deferrable) {
            const auto now = std::chrono::steady_clock::now();
            ++task.deferralCount;
            if (task.firstDeferred ==
                std::chrono::steady_clock::time_point()) {
                task.firstDeferred = now;
            }
            task.retryAfter = now + request.retryDelay();
            if (now - task.firstDeferred >= kMemoryStarvationAge) {
                task.estimatedCost = std::numeric_limits<uint64_t>::max();
            }
            {
                std::lock_guard<std::mutex> lock(mutex_);
                canDefer = !stopping_;
                if (canDefer) {
                    deferredTasks_.push_back(std::move(task));
                }
            }
            if (canDefer) {
                taskAvailable_.notify_all();
                deferred = true;
            }
        }
        if (!canDefer) {
            recordException(task.group, std::current_exception());
        }
    } catch (...) {
        recordException(task.group, std::current_exception());
    }
    currentTaskDeferrable = previousDeferrable;
    currentTaskDeferralCount = previousDeferralCount;
    currentTaskFirstDeferred = previousFirstDeferred;
    currentTaskRetainedState = previousRetainedState;

    if (!deferred && task.group != nullptr) {
        // Publish completion while holding the group mutex. waitForGroup()
        // takes the same mutex before returning, so the group cannot be
        // destroyed while this worker still accesses its condition variable.
        {
            std::lock_guard<std::mutex> groupLock(task.group->mutex_);
            task.group->remaining_.fetch_sub(1, std::memory_order_acq_rel);
            task.group->completed_.notify_all();
        }
        // Wake cooperative parent workers whose predicate also observes group
        // completion. This closes the notification race between nested groups.
        taskAvailable_.notify_all();
    }

    if (!deferred) {
        std::lock_guard<std::mutex> lock(mutex_);
        outstandingEstimatedCost_ = task.accountedCost >
                                    outstandingEstimatedCost_
                                    ? 0
                                    : outstandingEstimatedCost_ -
                                      task.accountedCost;
        const auto cost = outstandingCosts_.find(task.accountedCost);
        if (cost != outstandingCosts_.end()) {
            outstandingCosts_.erase(cost);
        }
    }

    if (countAsActiveTask) {
        std::lock_guard<std::mutex> lock(mutex_);
        --activeTasks_;
        if (tasks_.empty() && deferredTasks_.empty() && activeTasks_ == 0 &&
            plannedCosts_.empty()) {
            idle_.notify_all();
        }
    }
}

bool currentAnchorTaskCanDefer() {
    return currentAnchorExecutor != nullptr && currentTaskDeferrable;
}

bool currentAnchorTaskShouldDrainMemory() {
    if (currentAnchorExecutor == nullptr || !currentTaskDeferrable) {
        return false;
    }
    if (currentAnchorExecutor->loadSnapshot().admissionTailPhase) {
        return true;
    }
    return currentTaskFirstDeferred !=
                   std::chrono::steady_clock::time_point() &&
           std::chrono::steady_clock::now() - currentTaskFirstDeferred >=
                   kMemoryStarvationAge;
}

uint64_t currentAnchorTaskDeferralCount() {
    return currentTaskDeferralCount;
}

AnchorTaskExecutor::LoadSnapshot currentAnchorTaskLoadSnapshot() {
    return currentAnchorExecutor == nullptr
           ? AnchorTaskExecutor::LoadSnapshot{}
           : currentAnchorExecutor->loadSnapshot();
}

std::shared_ptr<void> currentAnchorTaskRetainedState() {
    return currentTaskRetainedState == nullptr
           ? std::shared_ptr<void>() : *currentTaskRetainedState;
}

void setCurrentAnchorTaskRetainedState(std::shared_ptr<void> state) {
    if (currentTaskRetainedState != nullptr) {
        *currentTaskRetainedState = std::move(state);
    }
}

void AnchorTaskExecutor::recordException(
        AnchorTaskGroup *group, std::exception_ptr exception) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!firstException_) {
            firstException_ = exception;
        }
    }
    if (group != nullptr) {
        std::lock_guard<std::mutex> lock(group->mutex_);
        if (!group->firstException_) {
            group->firstException_ = exception;
        }
    }
}

uint64_t anchorTaskEstimatedCost(std::size_t anchorCount) {
    const uint64_t count = static_cast<uint64_t>(anchorCount);
    if (count != 0 &&
        count > std::numeric_limits<uint64_t>::max() / count) {
        return std::numeric_limits<uint64_t>::max();
    }
    return count * count;
}

uint64_t anchorGapTaskEstimatedCost(
        std::size_t referenceLength, std::size_t queryLength) {
    const uint64_t reference = static_cast<uint64_t>(referenceLength);
    const uint64_t query = static_cast<uint64_t>(queryLength);
    if (query != 0 &&
        reference > std::numeric_limits<uint64_t>::max() / query) {
        return std::numeric_limits<uint64_t>::max();
    }
    return reference * query;
}

}  // namespace anchorwave
