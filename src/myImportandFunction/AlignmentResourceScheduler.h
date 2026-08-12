#pragma once

#include <condition_variable>
#include <chrono>
#include <cstdint>
#include <deque>
#include <functional>
#include <mutex>
#include <unordered_map>

namespace anchorwave {

struct AlignmentResourcePlan {
    int requestedMaxThreads = 1;
    int effectiveMaxThreads = 1;
    uint64_t maxProcessMemoryBytes = 0;
    uint64_t baselineResidentBytes = 0;
    uint64_t perWorkerMemoryBytes = 0;
    uint64_t perWorkerReservationBytes = 0;
    uint64_t safetyReserveBytes = 0;
    uint64_t taskMemoryCapacityBytes = 0;
};

// Convert a command-line memory limit expressed in GiB to bytes. Throws
// std::invalid_argument/std::overflow_error for invalid input.
uint64_t memoryLimitBytesFromGiB(double gibibytes);

// Return the current process resident set when the host provides a supported
// interface (Mach on macOS or /proc on Linux), otherwise zero.
uint64_t currentProcessResidentBytes();

// Build a process-wide reservation plan. A zero maxProcessMemoryBytes retains
// legacy thread-only scheduling. For an explicit cap, baseline resident memory
// and a safety reserve are removed from -M; the rest is shared dynamically by
// task-specific predicted-memory reservations.
AlignmentResourcePlan makeAlignmentResourcePlan(
        int requestedMaxThreads,
        uint64_t maxProcessMemoryBytes,
        int64_t windowWidth,
        uint64_t baselineResidentBytes);

struct AlignmentMemorySchedulerStats {
    uint64_t peakReservedBytes = 0;
    uint64_t peakProjectedProcessBytes = 0;
    uint64_t peakObservedResidentBytes = 0;
    uint64_t reservationCount = 0;
    uint64_t waitedReservationCount = 0;
    uint64_t temporaryReservationDeferrals = 0;
    uint64_t preferredReservationCount = 0;
    uint64_t preferredWaitCount = 0;
    uint64_t impossibleReservationCount = 0;
    int peakConcurrentReservations = 0;
};

struct AlignmentResourceSnapshot {
    bool enabled = false;
    uint64_t maxProcessMemoryBytes = 0;
    uint64_t safetyReserveBytes = 0;
    uint64_t taskMemoryCapacityBytes = 0;
    uint64_t activeReservedBytes = 0;
    uint64_t observedResidentBytes = 0;
    // Resident bytes observed above the baseline and all currently active
    // reservations. This covers completed alignment strings, WFA caches and
    // other post-baseline state whose owning task has already released its
    // reservation.
    uint64_t untrackedResidentBytes = 0;
    uint64_t immediatelyAvailableBytes = 0;
    uint64_t preferredQueuedBytes = 0;
    int activeReservations = 0;
    uint64_t preferredRequests = 0;
};

// Runtime memory feasibility is more precise than a binary fits/does-not-fit
// result.  In particular, a request that fitted the baseline captured at
// startup can become impossible after retained result/output RSS raises the
// process floor.  CoolingRss remains temporary; StableRuntimeFloor is emitted
// only after a bounded observation interval with no tracked releaser and no
// meaningful RSS decline.
enum class MemoryAvailability {
    Ready,
    TrackedWait,
    CoolingRss,
    StableRuntimeFloor,
    StaticInfeasible
};

struct MemoryAvailabilityEstimate {
    MemoryAvailability availability = MemoryAvailability::StaticInfeasible;
    // Ready=0; TrackedWait=predicted release time; CoolingRss=remaining
    // bounded observation grace; infeasible states=infinity.
    double waitMinutes = 0.0;
    uint64_t observedResidentBytes = 0;
    uint64_t projectedProcessBytes = 0;
    uint64_t shortfallBytes = 0;
};

// One state object belongs to one logical task/candidate reservation and is
// retained across executor deferrals.  Keeping it outside the global
// scheduler prevents unrelated requests of the same size from sharing an RSS
// aging decision and makes the classifier deterministic under a fake clock.
struct AlignmentMemoryPressureState {
    uint64_t requestedBytes = 0;
    uint64_t lowestObservedResidentBytes = 0;
    bool initialized = false;
    std::chrono::steady_clock::time_point firstBlocked;
    std::chrono::steady_clock::time_point lastMeaningfulRssDrop;

    void reset() {
        requestedBytes = 0;
        lowestObservedResidentBytes = 0;
        initialized = false;
        firstBlocked = std::chrono::steady_clock::time_point();
        lastMeaningfulRssDrop = std::chrono::steady_clock::time_point();
    }
};

class AlignmentMemoryScheduler;
struct AlignmentMemoryTryResult;

enum class AlignmentMemoryAdmission {
    Acquired,
    TemporarilyUnavailable,
    PermanentlyInfeasible
};

// Move-only RAII token returned by AlignmentMemoryScheduler. For the
// nonblocking/preferred APIs, inspect AlignmentMemoryTryResult::admission to
// distinguish structural infeasibility from current RSS pressure.
class AlignmentMemoryReservation {
public:
    AlignmentMemoryReservation() = default;
    ~AlignmentMemoryReservation();
    AlignmentMemoryReservation(AlignmentMemoryReservation &&other) noexcept;
    AlignmentMemoryReservation &operator=(
            AlignmentMemoryReservation &&other) noexcept;

    AlignmentMemoryReservation(const AlignmentMemoryReservation &) = delete;
    AlignmentMemoryReservation &operator=(
            const AlignmentMemoryReservation &) = delete;

    explicit operator bool() const;
    uint64_t reservedBytes() const;

private:
    friend class AlignmentMemoryScheduler;
    friend AlignmentMemoryReservation acquireAlignmentMemory(
            uint64_t predictedPeakBytes, double predictedRuntimeMinutes);
    friend AlignmentMemoryTryResult tryAcquireAlignmentMemory(
            uint64_t predictedPeakBytes, double predictedRuntimeMinutes);
    friend AlignmentMemoryTryResult acquirePreferredAlignmentMemory(
            uint64_t predictedPeakBytes, double predictedRuntimeMinutes);
    AlignmentMemoryReservation(AlignmentMemoryScheduler *scheduler,
                               uint64_t reservationId,
                               uint64_t reservedBytes, bool valid);
    void release();

    AlignmentMemoryScheduler *scheduler_ = nullptr;
    uint64_t reservationId_ = 0;
    uint64_t reservedBytes_ = 0;
    bool valid_ = false;
};

struct AlignmentMemoryTryResult {
    AlignmentMemoryAdmission admission =
            AlignmentMemoryAdmission::PermanentlyInfeasible;
    AlignmentMemoryReservation reservation;

    explicit operator bool() const {
        return admission == AlignmentMemoryAdmission::Acquired &&
               static_cast<bool>(reservation);
    }
};

// Process-wide memory admission controller for individual sequence-alignment
// attempts. Each engine's working-memory budget is independently capped at
// w^2; this scheduler admits its guarded prediction plus transient task bytes
// under -M. The safety reserve is subtracted before any task receives a token.
class AlignmentMemoryScheduler {
public:
    using ResidentMemoryReader = std::function<uint64_t()>;
    using Clock = std::chrono::steady_clock;
    using NowReader = std::function<Clock::time_point()>;

    explicit AlignmentMemoryScheduler(
            const AlignmentResourcePlan &plan,
            ResidentMemoryReader residentMemoryReader =
                    currentProcessResidentBytes,
            NowReader nowReader = NowReader());
    AlignmentMemoryScheduler(const AlignmentMemoryScheduler &) = delete;
    AlignmentMemoryScheduler &operator=(
            const AlignmentMemoryScheduler &) = delete;

    AlignmentMemoryReservation acquire(
            uint64_t predictedPeakBytes,
            double predictedRuntimeMinutes = 0.0);
    // Nonblocking admission used by deferrable alignment tasks. Temporary
    // unavailability is distinct from a request that can never fit under -M.
    AlignmentMemoryTryResult tryAcquire(
            uint64_t predictedPeakBytes,
            double predictedRuntimeMinutes = 0.0);
    // FIFO draining admission for an aged or tail-phase high-WFA request.
    // While the oldest preferred request waits, younger requests cannot
    // consume memory required to satisfy it. Once tracked reservations drain,
    // persistent RSS pressure is returned as TemporarilyUnavailable rather
    // than being mislabeled as a structurally impossible task.
    AlignmentMemoryTryResult acquirePreferred(
            uint64_t predictedPeakBytes,
            double predictedRuntimeMinutes = 0.0);
    // Predict the earliest tracked-reservation completion that creates enough
    // room.  Residual/untracked RSS remains charged throughout the simulation.
    double estimatedWaitMinutes(uint64_t predictedPeakBytes);
    // Structured availability used by the exact-tier dispatcher. Passing a
    // pressure state enables the 30-second stable-RSS-floor classifier;
    // nullptr retains the legacy non-aging query used by telemetry.
    MemoryAvailabilityEstimate memoryAvailability(
            uint64_t predictedPeakBytes,
            AlignmentMemoryPressureState *pressureState = nullptr);
    uint64_t maximumSingleTaskReservationBytes() const;
    AlignmentResourceSnapshot snapshot() const;
    AlignmentMemorySchedulerStats stats() const;
    bool enabled() const;

private:
    friend class AlignmentMemoryReservation;
    struct ActiveReservation {
        uint64_t bytes = 0;
        std::chrono::steady_clock::time_point predictedFinish;
    };
    struct PreferredRequest {
        uint64_t ticket = 0;
        uint64_t bytes = 0;
    };

    void release(uint64_t reservationId, uint64_t reservedBytes);
    uint64_t observedResidentBytesLocked();
    uint64_t projectedProcessBytesLocked() const;
    bool fitsLocked(uint64_t predictedPeakBytes,
                    bool respectPreferredQueue) const;
    AlignmentMemoryReservation admitLocked(
            uint64_t predictedPeakBytes, double predictedRuntimeMinutes);

    const AlignmentResourcePlan plan_;
    ResidentMemoryReader residentMemoryReader_;
    NowReader nowReader_;
    mutable std::mutex mutex_;
    std::condition_variable memoryAvailable_;
    uint64_t activeReservedBytes_ = 0;
    uint64_t untrackedResidentBytes_ = 0;
    uint64_t nextReservationId_ = 1;
    uint64_t nextPreferredTicket_ = 1;
    int activeReservations_ = 0;
    std::unordered_map<uint64_t, ActiveReservation> activeReservationInfo_;
    std::deque<PreferredRequest> preferredRequests_;
    AlignmentMemorySchedulerStats stats_;
};

// Install one scheduler for the duration of a genoAli/proali sequence stage.
// alignSlidingWindow() uses it from all worker threads. Calls outside such a
// scope retain the legacy unconstrained behavior.
class ScopedAlignmentMemoryScheduler {
public:
    explicit ScopedAlignmentMemoryScheduler(
            AlignmentMemoryScheduler &scheduler);
    ~ScopedAlignmentMemoryScheduler();
    ScopedAlignmentMemoryScheduler(
            const ScopedAlignmentMemoryScheduler &) = delete;
    ScopedAlignmentMemoryScheduler &operator=(
            const ScopedAlignmentMemoryScheduler &) = delete;

private:
    AlignmentMemoryScheduler *scheduler_ = nullptr;
};

AlignmentMemoryReservation acquireAlignmentMemory(
        uint64_t predictedPeakBytes,
        double predictedRuntimeMinutes = 0.0);
AlignmentMemoryTryResult tryAcquireAlignmentMemory(
        uint64_t predictedPeakBytes,
        double predictedRuntimeMinutes = 0.0);
AlignmentMemoryTryResult acquirePreferredAlignmentMemory(
        uint64_t predictedPeakBytes,
        double predictedRuntimeMinutes = 0.0);
double estimatedAlignmentMemoryWaitMinutes(uint64_t predictedPeakBytes);
MemoryAvailabilityEstimate estimateAlignmentMemoryAvailability(
        uint64_t predictedPeakBytes,
        AlignmentMemoryPressureState *pressureState = nullptr);
uint64_t maximumAlignmentTaskReservationBytes();
AlignmentResourceSnapshot alignmentResourceSnapshot();
bool alignmentMemorySchedulingEnabled();

// Conservative transient storage for encoded inputs, aligned output strings,
// CIGAR conversion, and allocator rounding. This amount is part of each task's
// reservation, rather than an allowance outside -M.
uint64_t alignmentTaskTransientMemoryBytes(uint64_t queryLength,
                                           uint64_t referenceLength);

// A blocking, dynamically replenished worker-slot limiter. Completed workers
// immediately return their slot to queued chromosomes/collinear blocks.
class DynamicThreadLimiter {
public:
    explicit DynamicThreadLimiter(int maximumActiveThreads);
    DynamicThreadLimiter(const DynamicThreadLimiter &) = delete;
    DynamicThreadLimiter &operator=(const DynamicThreadLimiter &) = delete;

    void acquire();
    void release();
    void waitUntilIdle();

    int maximumActiveThreads() const;
    int peakActiveThreads() const;

private:
    const int maximumActiveThreads_;
    mutable std::mutex mutex_;
    std::condition_variable condition_;
    int activeThreads_ = 0;
    int peakActiveThreads_ = 0;
};

}  // namespace anchorwave
