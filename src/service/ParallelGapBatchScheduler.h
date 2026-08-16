#pragma once

#include "AlignmentGapPlanner.h"
#include "AnchorTaskExecutor.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace anchorwave {

// Testable rolling-window core used by the genome-alignment gap scheduler.
// Each descriptor owns an independent task group and shared future. Waiting
// for one required result therefore retains AnchorTaskExecutor's cooperative
// nested-task execution without imposing a barrier on unrelated gaps.
//
// The initial runnable window contains at most maximumPendingResults tasks.
// Half protects the earliest anchor-order frontier and half prefetches the
// largest predicted tasks. Completion or memory deferral immediately frees a
// runnable slot, so a replacement may be submitted before ordered output
// consumes the old slot. Normal operation retains at most two windows. When
// ordered output hides a dormant future backlog and the process memory
// scheduler approves, a wider look-ahead window may be exposed. This
// changes queue visibility only: no more than one worker-width is runnable and
// every alignment attempt still passes the process-wide memory admission gate.
// The owning pipeline additionally limits retained completed-result bytes.
template <typename Result>
class ParallelGapBatchSchedulerCore {
public:
    using ResultPointer = std::shared_ptr<const Result>;
    using WorkFunction =
            std::function<Result(AlignmentGapDescriptor &)>;
    using ResultBytesFunction = std::function<uint64_t(const Result &)>;
    using EmergencyBackfillFunction =
            std::function<bool(uint64_t, std::size_t)>;

    ParallelGapBatchSchedulerCore(
            std::vector<AlignmentGapDescriptor> descriptors,
            std::size_t maximumPendingResults,
            AnchorTaskExecutor &executor,
            WorkFunction work,
            std::atomic<uint64_t> *submittedTaskCount = nullptr,
            ResultBytesFunction resultBytes = ResultBytesFunction(),
            EmergencyBackfillFunction emergencyBackfill =
                    EmergencyBackfillFunction(),
            uint64_t maximumCompletedResultBytes =
                    std::numeric_limits<uint64_t>::max())
            : descriptors_(std::move(descriptors)),
              maximumPendingResults_(maximumPendingResults),
              maximumOutstandingResults_(
                      maximumPendingResults == 1
                      ? 1 : saturatingDouble(maximumPendingResults)),
              maximumEmergencyOutstandingResults_(
                      maximumPendingResults == 1
                      ? 1 : saturatingDouble(saturatingDouble(
                                saturatingDouble(maximumPendingResults)))),
              executor_(executor),
              work_(std::move(work)),
              resultBytes_(std::move(resultBytes)),
              emergencyBackfill_(std::move(emergencyBackfill)),
              maximumCompletedResultBytes_(maximumCompletedResultBytes),
              submittedTaskCount_(submittedTaskCount),
              longestPlanned_(DescriptorPriorityLess{&descriptors_}) {
        if (maximumPendingResults_ == 0) {
            throw std::invalid_argument(
                    "parallel gap rolling window must be positive");
        }
        if (!work_) {
            throw std::invalid_argument(
                    "parallel gap work function must not be empty");
        }
    }

    ParallelGapBatchSchedulerCore(
            const ParallelGapBatchSchedulerCore &) = delete;
    ParallelGapBatchSchedulerCore &operator=(
            const ParallelGapBatchSchedulerCore &) = delete;

    ~ParallelGapBatchSchedulerCore() {
        {
            std::lock_guard<std::mutex> lock(stateMutex_);
            shuttingDown_ = true;
            backfillEnabled_ = false;
        }
        waitForPendingNoThrow();
        cancelUnsubmittedNoThrow();
    }

    // Kept separate from construction so a submission failure unwinds through
    // a fully constructed object whose destructor can safely drain tasks that
    // were submitted before the failure.
    void start() {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (started_) {
            throw std::logic_error(
                    "parallel gap rolling window already started");
        }
        started_ = true;
        states_.assign(descriptors_.size(), DescriptorState::Planned);
        std::vector<uint64_t> costs;
        costs.reserve(descriptors_.size());
        descriptorByAnchor_.reserve(descriptors_.size());
        for (std::size_t index = 0; index < descriptors_.size(); ++index) {
            const AlignmentGapDescriptor &descriptor = descriptors_[index];
            if (index > 0 &&
                descriptor.anchorIndex <= descriptors_[index - 1].anchorIndex) {
                throw std::invalid_argument(
                        "parallel gap descriptors must have strictly increasing anchor indices");
            }
            if (!descriptorByAnchor_.emplace(
                        descriptor.anchorIndex, index).second) {
                throw std::invalid_argument(
                        "parallel gap descriptors contain a duplicate anchor index");
            }
            costs.push_back(descriptor.estimatedCost);
            longestPlanned_.push(index);
        }
        executor_.registerPlannedTaskCosts(costs);
        plannedCostsRegistered_ = true;
        try {
            fillWindowLocked();
        } catch (...) {
            backfillEnabled_ = false;
            cancelUnsubmittedLockedNoThrow();
            throw;
        }
    }

    ResultPointer resultBeforeAnchor(std::size_t anchorIndex) {
        std::shared_ptr<ResultSlot> slot;
        std::size_t descriptorIndex = 0;
        {
            std::lock_guard<std::mutex> lock(stateMutex_);
            rethrowAsyncFailureLocked();
            auto found = pending_.find(anchorIndex);
            if (found == pending_.end()) {
                const auto descriptor = descriptorByAnchor_.find(anchorIndex);
                if (descriptor == descriptorByAnchor_.end()) {
                    return ResultPointer();
                }
                descriptorIndex = descriptor->second;
                if (states_[descriptorIndex] == DescriptorState::Planned) {
                    if (pending_.size() >=
                            maximumEmergencyOutstandingResults_) {
                        throw std::logic_error(
                                "ordered gap frontier exhausted its reserved scheduling slots");
                    }
                    submitOneLocked(descriptorIndex);
                    found = pending_.find(anchorIndex);
                }
                if (found == pending_.end()) {
                    throw std::logic_error(
                            "parallel gap descriptor lost before ordered consumption");
                }
            } else {
                descriptorIndex = descriptorByAnchor_.at(anchorIndex);
            }
            slot = found->second;
        }

        ResultPointer result;
        std::exception_ptr failure;
        {
            // This group contains exactly one logical result. If this call
            // occurs inside a chromosome worker, waitForGroup cooperatively
            // executes queued work and therefore cannot deadlock the fixed-size
            // executor.
            OrderedHeadWaitGuard orderedHead(executor_);
            try {
                executor_.waitForGroup(slot->group);
                result = slot->ready.get();
            } catch (...) {
                failure = std::current_exception();
            }
        }

        {
            std::lock_guard<std::mutex> lock(stateMutex_);
            const auto found = pending_.find(anchorIndex);
            if (found == pending_.end() || found->second != slot) {
                throw std::logic_error(
                        "parallel gap result was consumed more than once");
            }
            pending_.erase(found);
            consumeCompletedLocked(descriptorIndex, *slot);
            if (failure) {
                backfillEnabled_ = false;
            } else {
                fillWindowLocked();
                rethrowAsyncFailureLocked();
            }
        }
        if (failure) {
            std::rethrow_exception(failure);
        }
        return result;
    }

    std::size_t pendingResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return pending_.size();
    }

    std::size_t peakPendingResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return peakPendingResultCount_;
    }

    std::size_t submittedDescriptorCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return submittedDescriptorCount_;
    }

    std::size_t inFlightResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return inFlightResults_;
    }

    std::size_t completedResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return completedResults_;
    }

    uint64_t completedResultBytes() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return completedResultBytes_;
    }

    std::size_t deferredResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return deferredResults_;
    }

    std::size_t maximumOutstandingResultCount() const {
        return maximumOutstandingResults_;
    }

    std::size_t maximumEmergencyOutstandingResultCount() const {
        return maximumEmergencyOutstandingResults_;
    }

    std::size_t emergencyBackfillCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return emergencyBackfillCount_;
    }

    std::size_t peakInFlightResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return peakInFlightResults_;
    }

    std::size_t peakCompletedResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return peakCompletedResults_;
    }

    uint64_t peakCompletedResultBytes() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return peakCompletedResultBytes_;
    }

    std::size_t peakDeferredResultCount() const {
        std::lock_guard<std::mutex> lock(stateMutex_);
        return peakDeferredResults_;
    }

private:
    enum class DescriptorState : uint8_t {
        Planned,
        InFlight,
        Deferred,
        Completed,
        Consumed
    };

    struct DescriptorPriorityLess {
        const std::vector<AlignmentGapDescriptor> *descriptors = nullptr;

        bool operator()(std::size_t left, std::size_t right) const {
            const AlignmentGapDescriptor &leftDescriptor =
                    (*descriptors)[left];
            const AlignmentGapDescriptor &rightDescriptor =
                    (*descriptors)[right];
            const uint64_t leftPriority =
                    leftDescriptor.schedulingPriorityCost != 0
                    ? leftDescriptor.schedulingPriorityCost
                    : leftDescriptor.estimatedCost;
            const uint64_t rightPriority =
                    rightDescriptor.schedulingPriorityCost != 0
                    ? rightDescriptor.schedulingPriorityCost
                    : rightDescriptor.estimatedCost;
            if (leftPriority != rightPriority) {
                return leftPriority < rightPriority;
            }
            return leftDescriptor.anchorIndex > rightDescriptor.anchorIndex;
        }
    };

    class OrderedHeadWaitGuard {
    public:
        explicit OrderedHeadWaitGuard(AnchorTaskExecutor &executor)
                : executor_(executor) {
            executor_.orderedHeadWaitStarted();
        }

        ~OrderedHeadWaitGuard() {
            executor_.orderedHeadWaitFinished();
        }

        OrderedHeadWaitGuard(const OrderedHeadWaitGuard &) = delete;
        OrderedHeadWaitGuard &operator=(const OrderedHeadWaitGuard &) = delete;

    private:
        AnchorTaskExecutor &executor_;
    };

    struct ResultSlot {
        ResultSlot()
                : ready(promise.get_future().share()) {}

        AnchorTaskGroup group;
        std::promise<ResultPointer> promise;
        std::shared_future<ResultPointer> ready;
        uint64_t resultBytes = 0;
    };

    enum class TaskTerminalState : uint8_t {
        Completed,
        Failed
    };

    static std::size_t saturatingDouble(std::size_t value) {
        return value > std::numeric_limits<std::size_t>::max() / 2
               ? std::numeric_limits<std::size_t>::max() : value * 2;
    }

    static std::size_t saturatingAddSize(std::size_t first,
                                         std::size_t second) {
        return second > std::numeric_limits<std::size_t>::max() - first
               ? std::numeric_limits<std::size_t>::max() : first + second;
    }

    static uint64_t saturatingAdd(uint64_t first, uint64_t second) {
        return second > std::numeric_limits<uint64_t>::max() - first
               ? std::numeric_limits<uint64_t>::max() : first + second;
    }

    static uint64_t saturatingMultiply(uint64_t first, uint64_t second) {
        if (first == 0 || second == 0) {
            return 0;
        }
        return first > std::numeric_limits<uint64_t>::max() / second
               ? std::numeric_limits<uint64_t>::max() : first * second;
    }

    static uint64_t intervalLength(uint64_t start, uint64_t end) {
        return end < start ? 0 : saturatingAdd(end - start, 1);
    }

    static uint64_t conservativeResultBytes(
            const AlignmentGapDescriptor &descriptor) {
        // A global pairwise alignment has at most q+r columns and stores two
        // gapped strings. Include the Result object itself; callers may provide
        // a more precise ResultBytesFunction for non-alignment result types.
        const uint64_t inputBases = saturatingAdd(
                intervalLength(descriptor.referenceStart,
                               descriptor.referenceEnd),
                intervalLength(descriptor.queryStart,
                               descriptor.queryEnd));
        return saturatingAdd(
                static_cast<uint64_t>(sizeof(Result)),
                saturatingMultiply(2, inputBases));
    }

    void taskAttemptStarted(std::size_t descriptorIndex) noexcept {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (descriptorIndex >= states_.size()) {
            recordAsyncLogicFailureLocked(
                    "parallel gap attempt used an invalid descriptor index");
            return;
        }
        DescriptorState &state = states_[descriptorIndex];
        if (state == DescriptorState::Deferred) {
            if (deferredResults_ == 0) {
                recordAsyncLogicFailureLocked(
                        "parallel gap deferred-result count underflow");
                return;
            }
            --deferredResults_;
            ++inFlightResults_;
            peakInFlightResults_ = std::max(
                    peakInFlightResults_, inFlightResults_);
            state = DescriptorState::InFlight;
        } else if (state != DescriptorState::InFlight) {
            recordAsyncLogicFailureLocked(
                    "parallel gap descriptor started in an invalid state");
        }
    }

    void taskDeferred(std::size_t descriptorIndex) noexcept {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (descriptorIndex >= states_.size() ||
            states_[descriptorIndex] != DescriptorState::InFlight ||
            inFlightResults_ == 0) {
            recordAsyncLogicFailureLocked(
                    "parallel gap descriptor deferred in an invalid state");
            return;
        }
        states_[descriptorIndex] = DescriptorState::Deferred;
        --inFlightResults_;
        ++deferredResults_;
        peakDeferredResults_ = std::max(
                peakDeferredResults_, deferredResults_);
        tryBackfillLockedNoThrow();
    }

    void taskTerminated(std::size_t descriptorIndex,
                        const std::shared_ptr<ResultSlot> &slot,
                        uint64_t resultBytes,
                        TaskTerminalState terminalState) noexcept {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (descriptorIndex >= states_.size() ||
            states_[descriptorIndex] != DescriptorState::InFlight ||
            inFlightResults_ == 0 || !slot) {
            recordAsyncLogicFailureLocked(
                    "parallel gap descriptor completed in an invalid state");
            return;
        }
        states_[descriptorIndex] = DescriptorState::Completed;
        --inFlightResults_;
        ++completedResults_;
        slot->resultBytes = resultBytes;
        completedResultBytes_ = saturatingAdd(completedResultBytes_,
                                               resultBytes);
        peakCompletedResults_ = std::max(
                peakCompletedResults_, completedResults_);
        peakCompletedResultBytes_ = std::max(
                peakCompletedResultBytes_, completedResultBytes_);
        if (terminalState == TaskTerminalState::Failed) {
            backfillEnabled_ = false;
            return;
        }
        tryBackfillLockedNoThrow();
    }

    void consumeCompletedLocked(std::size_t descriptorIndex,
                                const ResultSlot &slot) {
        if (descriptorIndex >= states_.size() ||
            states_[descriptorIndex] != DescriptorState::Completed ||
            completedResults_ == 0) {
            throw std::logic_error(
                    "parallel gap result consumed before task completion");
        }
        states_[descriptorIndex] = DescriptorState::Consumed;
        --completedResults_;
        completedResultBytes_ = slot.resultBytes > completedResultBytes_
                                ? 0
                                : completedResultBytes_ - slot.resultBytes;
    }

    void tryBackfillLockedNoThrow() noexcept {
        if (!backfillEnabled_ || shuttingDown_ || asyncFailure_) {
            return;
        }
        try {
            fillWindowLocked();
        } catch (...) {
            if (!asyncFailure_) {
                asyncFailure_ = std::current_exception();
            }
            backfillEnabled_ = false;
        }
    }

    void recordAsyncLogicFailureLocked(const char *message) noexcept {
        if (!asyncFailure_) {
            try {
                asyncFailure_ = std::make_exception_ptr(
                        std::logic_error(message));
            } catch (...) {
                asyncFailure_ = std::current_exception();
            }
        }
        backfillEnabled_ = false;
    }

    void rethrowAsyncFailureLocked() const {
        if (asyncFailure_) {
            std::rethrow_exception(asyncFailure_);
        }
    }

    void fillWindowLocked() {
        if (!started_ || shuttingDown_ || !backfillEnabled_) {
            return;
        }
        // Completed alignments remain resident after their per-attempt memory
        // token is released.  Stop speculative refill once that retained
        // output reaches the owning pipeline's budget; consuming an ordered
        // result calls fillWindowLocked() again and resumes progress.
        if (completedResults_ > 0 &&
            completedResultBytes_ >= maximumCompletedResultBytes_) {
            return;
        }
        const std::size_t availableRunnableSlots =
                inFlightResults_ < maximumPendingResults_
                ? maximumPendingResults_ - inFlightResults_ : 0;
        std::size_t outstandingCeiling = maximumOutstandingResults_;
        if (maximumPendingResults_ > 1 &&
            pending_.size() >= maximumOutstandingResults_ &&
            pending_.size() < maximumEmergencyOutstandingResults_) {
            const AnchorTaskExecutor::LoadSnapshot load =
                    executor_.loadSnapshot();
            // Ordered output can fill every normal 2T result slot behind one
            // slow head while dormant future work could use idle CPUs. Expose
            // up to 8T descriptors only when the caller accepts the currently
            // retained result bytes and the process-wide memory state. The
            // runnable count remains T, so this is look-ahead rather than
            // thread oversubscription.
            if (emergencyBackfill_ &&
                emergencyBackfill_(completedResultBytes_, pending_.size()) &&
                load.deferredTasks == 0 &&
                load.globalFutureTasks > 0) {
                outstandingCeiling = maximumEmergencyOutstandingResults_;
            }
        }
        const std::size_t availableOutstandingSlots =
                pending_.size() < outstandingCeiling
                ? outstandingCeiling - pending_.size() : 0;
        std::size_t submissionsRemaining = std::min(
                availableRunnableSlots, availableOutstandingSlots);
        if (submissionsRemaining == 0) {
            return;
        }

        const std::size_t frontierSlots = maximumPendingResults_ == 1
                                          ? 1
                                          : (maximumPendingResults_ + 1) / 2;
        std::size_t frontierDescriptors = 0;
        for (std::size_t index = 0;
             index < descriptors_.size() &&
             frontierDescriptors < frontierSlots &&
             submissionsRemaining > 0;
             ++index) {
            if (states_[index] == DescriptorState::Consumed) {
                continue;
            }
            ++frontierDescriptors;
            if (states_[index] == DescriptorState::Planned) {
                submitOneLocked(index);
                --submissionsRemaining;
            }
        }
        while (submissionsRemaining > 0) {
            while (!longestPlanned_.empty() &&
                   states_[longestPlanned_.top()] !=
                           DescriptorState::Planned) {
                longestPlanned_.pop();
            }
            if (longestPlanned_.empty()) {
                break;
            }
            const std::size_t index = longestPlanned_.top();
            longestPlanned_.pop();
            submitOneLocked(index);
            --submissionsRemaining;
        }
    }

    void submitOneLocked(std::size_t descriptorIndex) {
        if (descriptorIndex >= descriptors_.size() ||
            states_[descriptorIndex] != DescriptorState::Planned) {
            throw std::logic_error(
                    "parallel gap descriptor submitted more than once");
        }
        if (pending_.size() >= maximumEmergencyOutstandingResults_) {
            throw std::logic_error(
                    "parallel gap outstanding-result ceiling exceeded");
        }
        const bool emergencySubmission =
                pending_.size() >= maximumOutstandingResults_;
        AlignmentGapDescriptor descriptor =
                descriptors_[descriptorIndex];
        const std::shared_ptr<ResultSlot> slot =
                std::make_shared<ResultSlot>();
        const auto inserted = pending_.emplace(descriptor.anchorIndex, slot);
        if (!inserted.second) {
            throw std::logic_error(
                    "parallel gap descriptor already has a pending result");
        }

        const WorkFunction work = work_;
        const ResultBytesFunction resultBytes = resultBytes_;
        states_[descriptorIndex] = DescriptorState::InFlight;
        ++inFlightResults_;
        peakInFlightResults_ = std::max(
                peakInFlightResults_, inFlightResults_);
        try {
            executor_.submitPlannedDeferrable(
                    slot->group, descriptor.estimatedCost,
                    [this, slot, descriptor, descriptorIndex, work,
                     resultBytes]() mutable {
                        taskAttemptStarted(descriptorIndex);
                        try {
                            Result value = work(descriptor);
                            const uint64_t bytes = resultBytes
                                    ? resultBytes(value)
                                    : conservativeResultBytes(descriptor);
                            slot->promise.set_value(
                                    std::make_shared<const Result>(
                                            std::move(value)));
                            taskTerminated(
                                    descriptorIndex, slot, bytes,
                                    TaskTerminalState::Completed);
                        } catch (const AlignmentTaskDeferred &) {
                            taskDeferred(descriptorIndex);
                            // A deferred task will execute this same closure
                            // again. Its promise must remain unsatisfied until
                            // an actual result or terminal failure occurs.
                            throw;
                        } catch (...) {
                            const std::exception_ptr failure =
                                    std::current_exception();
                            try {
                                slot->promise.set_exception(failure);
                            } catch (...) {
                                // Preserve the original alignment exception;
                                // the executor records it for this group and
                                // the outer waitForIdle().
                            }
                            taskTerminated(
                                    descriptorIndex, slot, 0,
                                    TaskTerminalState::Failed);
                            std::rethrow_exception(failure);
                        }
                    });
        } catch (...) {
            states_[descriptorIndex] = DescriptorState::Planned;
            --inFlightResults_;
            pending_.erase(inserted.first);
            throw;
        }
        ++submittedDescriptorCount_;
        if (emergencySubmission) {
            ++emergencyBackfillCount_;
        }
        peakPendingResultCount_ = std::max(
                peakPendingResultCount_, pending_.size());
        if (submittedTaskCount_ != nullptr) {
            submittedTaskCount_->fetch_add(1, std::memory_order_relaxed);
        }
    }

    void waitForPendingNoThrow() noexcept {
        std::vector<std::shared_ptr<ResultSlot>> pending;
        {
            std::lock_guard<std::mutex> lock(stateMutex_);
            pending.reserve(pending_.size());
            for (const auto &entry : pending_) {
                pending.push_back(entry.second);
            }
        }
        for (const std::shared_ptr<ResultSlot> &slot : pending) {
            try {
                executor_.waitForGroup(slot->group);
            } catch (...) {
                // AnchorTaskExecutor also retains the first terminal failure;
                // the owning stage's waitForIdle() propagates it after this
                // scheduler has safely released every task group.
            }
        }
    }

    void cancelUnsubmittedLockedNoThrow() noexcept {
        if (!plannedCostsRegistered_) {
            return;
        }
        std::vector<uint64_t> remaining;
        remaining.reserve(descriptors_.size());
        for (std::size_t index = 0; index < descriptors_.size(); ++index) {
            if (states_[index] == DescriptorState::Planned) {
                remaining.push_back(descriptors_[index].estimatedCost);
            }
        }
        executor_.cancelPlannedTaskCosts(remaining);
        plannedCostsRegistered_ = false;
    }

    void cancelUnsubmittedNoThrow() noexcept {
        std::lock_guard<std::mutex> lock(stateMutex_);
        cancelUnsubmittedLockedNoThrow();
    }

    std::vector<AlignmentGapDescriptor> descriptors_;
    const std::size_t maximumPendingResults_;
    const std::size_t maximumOutstandingResults_;
    const std::size_t maximumEmergencyOutstandingResults_;
    AnchorTaskExecutor &executor_;
    WorkFunction work_;
    ResultBytesFunction resultBytes_;
    EmergencyBackfillFunction emergencyBackfill_;
    const uint64_t maximumCompletedResultBytes_;
    std::atomic<uint64_t> *submittedTaskCount_;
    mutable std::mutex stateMutex_;
    std::size_t peakPendingResultCount_ = 0;
    std::size_t inFlightResults_ = 0;
    std::size_t completedResults_ = 0;
    uint64_t completedResultBytes_ = 0;
    std::size_t deferredResults_ = 0;
    std::size_t peakInFlightResults_ = 0;
    std::size_t peakCompletedResults_ = 0;
    uint64_t peakCompletedResultBytes_ = 0;
    std::size_t peakDeferredResults_ = 0;
    std::unordered_map<std::size_t, std::shared_ptr<ResultSlot>> pending_;
    std::unordered_map<std::size_t, std::size_t> descriptorByAnchor_;
    std::vector<DescriptorState> states_;
    std::priority_queue<std::size_t, std::vector<std::size_t>,
                        DescriptorPriorityLess> longestPlanned_;
    std::size_t submittedDescriptorCount_ = 0;
    std::size_t emergencyBackfillCount_ = 0;
    std::exception_ptr asyncFailure_;
    bool started_ = false;
    bool plannedCostsRegistered_ = false;
    bool backfillEnabled_ = true;
    bool shuttingDown_ = false;
};

}  // namespace anchorwave
