#include "AlignmentResourceScheduler.h"

#include "WfaAlignment.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

namespace anchorwave {

namespace {
constexpr uint64_t kPerWorkerNonWfaHeadroom = 512ULL * 1024ULL * 1024ULL;
constexpr uint64_t kMinimumSafetyReserve = 1ULL * 1024ULL * 1024ULL * 1024ULL;
constexpr uint64_t kMaximumSafetyReserve = 32ULL * 1024ULL * 1024ULL * 1024ULL;
constexpr uint64_t kSafetyBytesPerWorker = 64ULL * 1024ULL * 1024ULL;
constexpr uint64_t kMinimumTaskCapacity = 128ULL * 1024ULL * 1024ULL;
constexpr uint64_t kTaskTransientBase = 16ULL * 1024ULL * 1024ULL;
constexpr uint64_t kTaskBytesPerInputBase = 6;
constexpr uint64_t kMinimumMeaningfulRssDrop =
        256ULL * 1024ULL * 1024ULL;
constexpr auto kStableRssFloorAge = std::chrono::seconds(30);
constexpr auto kStableRssNoProgressAge = std::chrono::seconds(5);
constexpr double kCoolingRssPollMinutes = 1.0 / 60.0;
// A fast exact candidate must get one bounded cooling opportunity before a
// memory-light, minutes-long exact fallback starts.  Round1 interval 1951
// showed that infinity here immediately selected a 927-second BiWFA although
// the 18.8-second KSW2-full became admissible as other tasks released pages.
constexpr double kActiveCoolingDominanceGraceMinutes = 0.5;
std::atomic<AlignmentMemoryScheduler *> installedMemoryScheduler(nullptr);

uint64_t saturatingAdd(uint64_t first, uint64_t second) {
    if (second > std::numeric_limits<uint64_t>::max() - first) {
        return std::numeric_limits<uint64_t>::max();
    }
    return first + second;
}
}

uint64_t memoryLimitBytesFromGiB(double gibibytes) {
    if (!std::isfinite(gibibytes) || gibibytes <= 0.0) {
        throw std::invalid_argument("memory limit must be a positive number");
    }
    constexpr long double bytesPerGiB = 1024.0L * 1024.0L * 1024.0L;
    const long double bytes =
            static_cast<long double>(gibibytes) * bytesPerGiB;
    if (bytes > static_cast<long double>(
                    std::numeric_limits<uint64_t>::max())) {
        throw std::overflow_error("memory limit is too large");
    }
    return static_cast<uint64_t>(bytes);
}

uint64_t currentProcessResidentBytes() {
#if defined(__APPLE__)
    mach_task_basic_info_data_t information{};
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    const kern_return_t status = task_info(
            mach_task_self(), MACH_TASK_BASIC_INFO,
            reinterpret_cast<task_info_t>(&information), &count);
    return status == KERN_SUCCESS
           ? static_cast<uint64_t>(information.resident_size)
           : 0;
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm");
    uint64_t totalPages = 0;
    uint64_t residentPages = 0;
    if (!(statm >> totalPages >> residentPages)) {
        return 0;
    }
    const long pageSize = sysconf(_SC_PAGESIZE);
    if (pageSize <= 0 ||
        residentPages > std::numeric_limits<uint64_t>::max() /
                        static_cast<uint64_t>(pageSize)) {
        return 0;
    }
    return residentPages * static_cast<uint64_t>(pageSize);
#else
    return 0;
#endif
}

AlignmentResourcePlan makeAlignmentResourcePlan(
        int requestedMaxThreads,
        uint64_t maxProcessMemoryBytes,
        int64_t windowWidth,
        uint64_t baselineResidentBytes) {
    if (requestedMaxThreads <= 0) {
        throw std::invalid_argument("maximum thread count must be positive");
    }
    const uint64_t perWorkerMemoryBytes = wfaMemoryBudgetBytes(windowWidth);
    if (perWorkerMemoryBytes == 0) {
        throw std::invalid_argument("alignment window must be positive");
    }

    AlignmentResourcePlan plan;
    plan.requestedMaxThreads = requestedMaxThreads;
    plan.effectiveMaxThreads = requestedMaxThreads;
    plan.maxProcessMemoryBytes = maxProcessMemoryBytes;
    plan.baselineResidentBytes = baselineResidentBytes;
    plan.perWorkerMemoryBytes = perWorkerMemoryBytes;
    plan.perWorkerReservationBytes =
            perWorkerMemoryBytes > std::numeric_limits<uint64_t>::max() -
                                   kPerWorkerNonWfaHeadroom
            ? std::numeric_limits<uint64_t>::max()
            : perWorkerMemoryBytes + kPerWorkerNonWfaHeadroom;

    if (maxProcessMemoryBytes == 0) {
        return plan;
    }
    if (maxProcessMemoryBytes < perWorkerMemoryBytes) {
        throw std::invalid_argument(
                "maximum process memory must be at least w^2 bytes");
    }
    if (baselineResidentBytes >= maxProcessMemoryBytes) {
        throw std::runtime_error(
                "baseline resident memory already reaches the process limit");
    }
    const uint64_t remaining = maxProcessMemoryBytes - baselineResidentBytes;
    // Keep the reserve inside -M, but do not strand 50 GiB merely because a
    // node exposes 1 TiB.  Per-worker headroom covers concurrent allocator and
    // output transients; the percentage term covers model residuals.
    const uint64_t workerReserve = requestedMaxThreads > 0 &&
            static_cast<uint64_t>(requestedMaxThreads) >
                    std::numeric_limits<uint64_t>::max() /
                            kSafetyBytesPerWorker
            ? std::numeric_limits<uint64_t>::max()
            : static_cast<uint64_t>(requestedMaxThreads) *
                    kSafetyBytesPerWorker;
    plan.safetyReserveBytes = std::min<uint64_t>(
            kMaximumSafetyReserve,
            std::max<uint64_t>(
                    kMinimumSafetyReserve,
                    saturatingAdd(maxProcessMemoryBytes / 40,
                                  workerReserve)));
    if (remaining <= plan.safetyReserveBytes ||
        remaining - plan.safetyReserveBytes < kMinimumTaskCapacity) {
        throw std::runtime_error(
                "not enough memory remains after the in-limit safety reserve");
    }
    plan.taskMemoryCapacityBytes = remaining - plan.safetyReserveBytes;
    // Memory now controls individual algorithm attempts, not an up-front
    // number of w^2 workers. Keep all requested workers available so many
    // small predicted tasks can run concurrently.
    plan.effectiveMaxThreads = requestedMaxThreads;
    return plan;
}

AlignmentMemoryReservation::AlignmentMemoryReservation(
        AlignmentMemoryScheduler *scheduler, uint64_t reservationId,
        uint64_t reservedBytes, bool valid)
        : scheduler_(scheduler), reservationId_(reservationId),
          reservedBytes_(reservedBytes), valid_(valid) {}

AlignmentMemoryReservation::~AlignmentMemoryReservation() {
    release();
}

AlignmentMemoryReservation::AlignmentMemoryReservation(
        AlignmentMemoryReservation &&other) noexcept
        : scheduler_(other.scheduler_), reservedBytes_(other.reservedBytes_),
          valid_(other.valid_) {
    reservationId_ = other.reservationId_;
    other.scheduler_ = nullptr;
    other.reservationId_ = 0;
    other.reservedBytes_ = 0;
    other.valid_ = false;
}

AlignmentMemoryReservation &AlignmentMemoryReservation::operator=(
        AlignmentMemoryReservation &&other) noexcept {
    if (this != &other) {
        release();
        scheduler_ = other.scheduler_;
        reservationId_ = other.reservationId_;
        reservedBytes_ = other.reservedBytes_;
        valid_ = other.valid_;
        other.scheduler_ = nullptr;
        other.reservationId_ = 0;
        other.reservedBytes_ = 0;
        other.valid_ = false;
    }
    return *this;
}

AlignmentMemoryReservation::operator bool() const {
    return valid_;
}

uint64_t AlignmentMemoryReservation::reservedBytes() const {
    return reservedBytes_;
}

void AlignmentMemoryReservation::release() {
    if (scheduler_ != nullptr && valid_) {
        scheduler_->release(reservationId_, reservedBytes_);
    }
    scheduler_ = nullptr;
    reservationId_ = 0;
    reservedBytes_ = 0;
    valid_ = false;
}

AlignmentMemoryScheduler::AlignmentMemoryScheduler(
        const AlignmentResourcePlan &plan,
        ResidentMemoryReader residentMemoryReader,
        NowReader nowReader)
        : plan_(plan), residentMemoryReader_(std::move(residentMemoryReader)),
          nowReader_(std::move(nowReader)) {
    if (!residentMemoryReader_) {
        residentMemoryReader_ = currentProcessResidentBytes;
    }
    if (!nowReader_) {
        nowReader_ = []() { return Clock::now(); };
    }
}

bool AlignmentMemoryScheduler::enabled() const {
    return plan_.maxProcessMemoryBytes != 0;
}

uint64_t AlignmentMemoryScheduler::observedResidentBytesLocked() {
    const uint64_t observed = residentMemoryReader_();
    stats_.peakObservedResidentBytes = std::max(
            stats_.peakObservedResidentBytes, observed);
    const uint64_t accounted = saturatingAdd(
            plan_.baselineResidentBytes, activeReservedBytes_);
    const uint64_t currentUntracked = observed > accounted
                                      ? observed - accounted : 0;
    // With no active algorithm, RSS above the baseline is directly
    // attributable to retained output/caches and may also shrink after the
    // consumer advances. While algorithms are active, keep the largest
    // observed residual so a later reservation cannot overlap it.
    if (activeReservations_ == 0) {
        untrackedResidentBytes_ = currentUntracked;
    } else {
        untrackedResidentBytes_ = std::max(
                untrackedResidentBytes_, currentUntracked);
    }
    return observed;
}

uint64_t AlignmentMemoryScheduler::projectedProcessBytesLocked() const {
    return saturatingAdd(
            saturatingAdd(plan_.baselineResidentBytes,
                          activeReservedBytes_),
            untrackedResidentBytes_);
}

bool AlignmentMemoryScheduler::fitsLocked(
        uint64_t predictedPeakBytes,
        bool respectPreferredQueue) const {
    const uint64_t observed = residentMemoryReader_
                              ? residentMemoryReader_()
                              : 0;
    const uint64_t projected = std::max(
            observed, projectedProcessBytesLocked());
    const uint64_t usableLimit = plan_.maxProcessMemoryBytes -
                                 plan_.safetyReserveBytes;
    if (projected >= usableLimit ||
        predictedPeakBytes > usableLimit - projected) {
        return false;
    }
    if (respectPreferredQueue && !preferredRequests_.empty()) {
        const uint64_t remainingAfterAdmission =
                usableLimit - projected - predictedPeakBytes;
        if (preferredRequests_.front().bytes > remainingAfterAdmission) {
            return false;
        }
    }
    return true;
}

AlignmentMemoryReservation AlignmentMemoryScheduler::admitLocked(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    const uint64_t reservationId = nextReservationId_++;
    const double safeMinutes = std::isfinite(predictedRuntimeMinutes) &&
                               predictedRuntimeMinutes > 0.0
                               ? predictedRuntimeMinutes : 1.0;
    const auto duration = std::chrono::duration_cast<
            std::chrono::steady_clock::duration>(
                    std::chrono::duration<double, std::ratio<60>>(
                            safeMinutes));
    activeReservedBytes_ = saturatingAdd(
            activeReservedBytes_, predictedPeakBytes);
    ++activeReservations_;
    ++stats_.reservationCount;
    activeReservationInfo_.emplace(
            reservationId,
            ActiveReservation{predictedPeakBytes,
                              nowReader_() + duration});
    stats_.peakReservedBytes = std::max(
            stats_.peakReservedBytes, activeReservedBytes_);
    stats_.peakConcurrentReservations = std::max(
            stats_.peakConcurrentReservations, activeReservations_);
    stats_.peakProjectedProcessBytes = std::max(
            stats_.peakProjectedProcessBytes,
            projectedProcessBytesLocked());
    return AlignmentMemoryReservation(
            this, reservationId, predictedPeakBytes, true);
}

AlignmentMemoryTryResult AlignmentMemoryScheduler::tryAcquire(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    if (!enabled()) {
        AlignmentMemoryTryResult result;
        result.admission = AlignmentMemoryAdmission::Acquired;
        result.reservation = AlignmentMemoryReservation(
                nullptr, 0, predictedPeakBytes, true);
        return result;
    }
    predictedPeakBytes = std::max<uint64_t>(predictedPeakBytes, 1);
    if (predictedPeakBytes > plan_.taskMemoryCapacityBytes) {
        std::lock_guard<std::mutex> lock(mutex_);
        ++stats_.impossibleReservationCount;
        return AlignmentMemoryTryResult{
                AlignmentMemoryAdmission::PermanentlyInfeasible,
                AlignmentMemoryReservation()};
    }

    std::lock_guard<std::mutex> lock(mutex_);
    observedResidentBytesLocked();
    if (fitsLocked(predictedPeakBytes, true)) {
        AlignmentMemoryTryResult result;
        result.admission = AlignmentMemoryAdmission::Acquired;
        result.reservation = admitLocked(
                predictedPeakBytes, predictedRuntimeMinutes);
        return result;
    }
    // A request smaller than taskMemoryCapacityBytes is structurally viable.
    // Elevated RSS with no active reservation may be retained WFA cache or
    // output backlog, so it is temporary pressure rather than proof that this
    // candidate can never run.
    ++stats_.temporaryReservationDeferrals;
    return AlignmentMemoryTryResult{
            AlignmentMemoryAdmission::TemporarilyUnavailable,
            AlignmentMemoryReservation()};
}

AlignmentMemoryReservation AlignmentMemoryScheduler::acquire(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    if (!enabled()) {
        return AlignmentMemoryReservation(
                nullptr, 0, predictedPeakBytes, true);
    }
    predictedPeakBytes = std::max<uint64_t>(predictedPeakBytes, 1);
    if (predictedPeakBytes > plan_.taskMemoryCapacityBytes) {
        std::lock_guard<std::mutex> lock(mutex_);
        ++stats_.impossibleReservationCount;
        return AlignmentMemoryReservation();
    }

    std::unique_lock<std::mutex> lock(mutex_);
    bool waited = false;
    while (true) {
        observedResidentBytesLocked();
        if (fitsLocked(predictedPeakBytes, true)) {
            if (waited) {
                ++stats_.waitedReservationCount;
            }
            return admitLocked(predictedPeakBytes, predictedRuntimeMinutes);
        }
        if (activeReservations_ == 0 && preferredRequests_.empty()) {
            ++stats_.impossibleReservationCount;
            return AlignmentMemoryReservation();
        }
        waited = true;
        memoryAvailable_.wait(lock);
    }
}

AlignmentMemoryTryResult AlignmentMemoryScheduler::acquirePreferred(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    if (!enabled()) {
        AlignmentMemoryTryResult result;
        result.admission = AlignmentMemoryAdmission::Acquired;
        result.reservation = AlignmentMemoryReservation(
                nullptr, 0, predictedPeakBytes, true);
        return result;
    }
    predictedPeakBytes = std::max<uint64_t>(predictedPeakBytes, 1);
    if (predictedPeakBytes > plan_.taskMemoryCapacityBytes) {
        std::lock_guard<std::mutex> lock(mutex_);
        ++stats_.impossibleReservationCount;
        return AlignmentMemoryTryResult{
                AlignmentMemoryAdmission::PermanentlyInfeasible,
                AlignmentMemoryReservation()};
    }

    std::unique_lock<std::mutex> lock(mutex_);
    const uint64_t ticket = nextPreferredTicket_++;
    preferredRequests_.push_back(PreferredRequest{ticket,
                                                  predictedPeakBytes});
    ++stats_.preferredReservationCount;
    bool waited = false;
    while (true) {
        observedResidentBytesLocked();
        if (!preferredRequests_.empty() &&
            preferredRequests_.front().ticket == ticket &&
            fitsLocked(predictedPeakBytes, false)) {
            preferredRequests_.pop_front();
            if (waited) {
                ++stats_.preferredWaitCount;
                ++stats_.waitedReservationCount;
            }
            AlignmentMemoryTryResult result;
            result.admission = AlignmentMemoryAdmission::Acquired;
            result.reservation = admitLocked(
                    predictedPeakBytes, predictedRuntimeMinutes);
            memoryAvailable_.notify_all();
            return result;
        }
        if (activeReservations_ == 0 &&
            !fitsLocked(predictedPeakBytes, false)) {
            const auto found = std::find_if(
                    preferredRequests_.begin(), preferredRequests_.end(),
                    [ticket](const PreferredRequest &request) {
                        return request.ticket == ticket;
                    });
            if (found != preferredRequests_.end()) {
                preferredRequests_.erase(found);
            }
            ++stats_.temporaryReservationDeferrals;
            memoryAvailable_.notify_all();
            return AlignmentMemoryTryResult{
                    AlignmentMemoryAdmission::TemporarilyUnavailable,
                    AlignmentMemoryReservation()};
        }
        waited = true;
        memoryAvailable_.wait(lock);
    }
}

MemoryAvailabilityEstimate AlignmentMemoryScheduler::memoryAvailability(
        uint64_t predictedPeakBytes,
        AlignmentMemoryPressureState *pressureState) {
    MemoryAvailabilityEstimate result;
    predictedPeakBytes = std::max<uint64_t>(predictedPeakBytes, 1);
    if (!enabled()) {
        result.availability = MemoryAvailability::Ready;
        result.waitMinutes = 0.0;
        if (pressureState != nullptr) {
            pressureState->reset();
        }
        return result;
    }

    std::lock_guard<std::mutex> lock(mutex_);
    const auto now = nowReader_();
    result.observedResidentBytes = observedResidentBytesLocked();
    result.projectedProcessBytes = std::max(
            result.observedResidentBytes, projectedProcessBytesLocked());
    const uint64_t usableLimit = plan_.maxProcessMemoryBytes -
                                 plan_.safetyReserveBytes;
    if (predictedPeakBytes > plan_.taskMemoryCapacityBytes) {
        result.availability = MemoryAvailability::StaticInfeasible;
        result.waitMinutes = std::numeric_limits<double>::infinity();
        result.shortfallBytes = predictedPeakBytes -
                                plan_.taskMemoryCapacityBytes;
        if (pressureState != nullptr) {
            pressureState->reset();
        }
        return result;
    }

    const uint64_t queuedPreferredBytes = preferredRequests_.empty()
                                          ? 0
                                          : preferredRequests_.front().bytes;
    const uint64_t required = saturatingAdd(predictedPeakBytes,
                                             queuedPreferredBytes);
    const uint64_t immediatelyAvailable =
            result.projectedProcessBytes < usableLimit
            ? usableLimit - result.projectedProcessBytes : 0;
    result.shortfallBytes = required > immediatelyAvailable
                            ? required - immediatelyAvailable : 0;
    if (result.shortfallBytes == 0) {
        result.availability = MemoryAvailability::Ready;
        result.waitMinutes = 0.0;
        if (pressureState != nullptr) {
            pressureState->reset();
        }
        return result;
    }

    // Simulate exactly the tracked-reservation releases used by admission.
    // Untracked RSS and the oldest preferred request remain charged throughout
    // the simulation.  A finite result is therefore a real tracked wait, not
    // an assumption that allocator caches will disappear.
    const uint64_t simulatedUntracked = untrackedResidentBytes_;
    std::vector<ActiveReservation> active;
    active.reserve(activeReservationInfo_.size());
    for (const auto &item : activeReservationInfo_) {
        active.push_back(item.second);
    }
    std::sort(active.begin(), active.end(),
              [](const ActiveReservation &left,
                 const ActiveReservation &right) {
                  return left.predictedFinish < right.predictedFinish;
              });
    uint64_t simulatedReserved = activeReservedBytes_;
    for (const ActiveReservation &reservation : active) {
        simulatedReserved = reservation.bytes > simulatedReserved
                            ? 0 : simulatedReserved - reservation.bytes;
        const uint64_t projected = saturatingAdd(
                saturatingAdd(plan_.baselineResidentBytes,
                              simulatedReserved),
                simulatedUntracked);
        if (projected < usableLimit && required <= usableLimit - projected) {
            const auto remaining = reservation.predictedFinish > now
                                   ? reservation.predictedFinish - now
                                   : Clock::duration::zero();
            result.availability = MemoryAvailability::TrackedWait;
            result.waitMinutes = std::chrono::duration<
                    double, std::ratio<60>>(remaining).count();
            if (pressureState != nullptr) {
                pressureState->reset();
            }
            return result;
        }
    }

    // Active algorithms or an older preferred ticket are possible releasers.
    // Do not age the request towards a stable floor while such an event still
    // exists, even if the current conservative residual-RSS projection says
    // that the tracked release alone is insufficient.
    if (activeReservations_ > 0 || !preferredRequests_.empty()) {
        result.availability = MemoryAvailability::CoolingRss;
        // The tracked releases alone are insufficient under the conservative
        // residual-RSS projection, but completing aligners can also destroy
        // large private WFA allocators and return their resident pages. Give a
        // fast exact engine one bounded grace window so the dominance policy
        // can park instead of immediately starting a minutes-long succinct
        // mode. A genuinely long tracked wait was already returned above and
        // is still compared normally (for example interval 1973 at 19 min).
        result.waitMinutes = kActiveCoolingDominanceGraceMinutes;
        if (pressureState != nullptr) {
            pressureState->reset();
        }
        return result;
    }

    if (pressureState == nullptr) {
        result.availability = MemoryAvailability::CoolingRss;
        result.waitMinutes = kCoolingRssPollMinutes;
        return result;
    }

    if (!pressureState->initialized ||
        pressureState->requestedBytes != predictedPeakBytes) {
        pressureState->reset();
        pressureState->initialized = true;
        pressureState->requestedBytes = predictedPeakBytes;
        pressureState->lowestObservedResidentBytes =
                result.observedResidentBytes;
        pressureState->firstBlocked = now;
        pressureState->lastMeaningfulRssDrop = now;
    } else {
        const uint64_t meaningfulDrop = std::max<uint64_t>(
                kMinimumMeaningfulRssDrop,
                plan_.maxProcessMemoryBytes / 100);
        if (result.observedResidentBytes <
                    pressureState->lowestObservedResidentBytes &&
            pressureState->lowestObservedResidentBytes -
                    result.observedResidentBytes >= meaningfulDrop) {
            pressureState->lowestObservedResidentBytes =
                    result.observedResidentBytes;
            pressureState->lastMeaningfulRssDrop = now;
        }
    }

    const auto blockedAge = now - pressureState->firstBlocked;
    const auto noProgressAge = now -
                               pressureState->lastMeaningfulRssDrop;
    if (blockedAge >= kStableRssFloorAge &&
        noProgressAge >= kStableRssNoProgressAge) {
        result.availability = MemoryAvailability::StableRuntimeFloor;
        result.waitMinutes = std::numeric_limits<double>::infinity();
        return result;
    }

    result.availability = MemoryAvailability::CoolingRss;
    const auto remaining = blockedAge < kStableRssFloorAge
                           ? kStableRssFloorAge - blockedAge
                           : Clock::duration::zero();
    result.waitMinutes = std::max(
            kCoolingRssPollMinutes,
            std::chrono::duration<double, std::ratio<60>>(
                    remaining).count());
    return result;
}

double AlignmentMemoryScheduler::estimatedWaitMinutes(
        uint64_t predictedPeakBytes) {
    return memoryAvailability(predictedPeakBytes, nullptr).waitMinutes;
}

uint64_t AlignmentMemoryScheduler::maximumSingleTaskReservationBytes() const {
    return enabled() ? plan_.taskMemoryCapacityBytes
                     : plan_.perWorkerMemoryBytes;
}

AlignmentResourceSnapshot AlignmentMemoryScheduler::snapshot() const {
    AlignmentResourceSnapshot result;
    std::lock_guard<std::mutex> lock(mutex_);
    result.enabled = enabled();
    result.maxProcessMemoryBytes = plan_.maxProcessMemoryBytes;
    result.safetyReserveBytes = plan_.safetyReserveBytes;
    result.taskMemoryCapacityBytes = plan_.taskMemoryCapacityBytes;
    result.activeReservedBytes = activeReservedBytes_;
    result.untrackedResidentBytes = untrackedResidentBytes_;
    result.activeReservations = activeReservations_;
    result.preferredRequests = preferredRequests_.size();
    for (const PreferredRequest &request : preferredRequests_) {
        result.preferredQueuedBytes = saturatingAdd(
                result.preferredQueuedBytes, request.bytes);
    }
    result.observedResidentBytes = residentMemoryReader_
                                   ? residentMemoryReader_() : 0;
    if (!enabled()) {
        result.immediatelyAvailableBytes =
                std::numeric_limits<uint64_t>::max();
        return result;
    }
    const uint64_t usableLimit = plan_.maxProcessMemoryBytes -
                                 plan_.safetyReserveBytes;
    const uint64_t projected = std::max(
            result.observedResidentBytes,
            projectedProcessBytesLocked());
    result.immediatelyAvailableBytes = projected < usableLimit
                                       ? usableLimit - projected : 0;
    return result;
}

void AlignmentMemoryScheduler::release(
        uint64_t reservationId, uint64_t reservedBytes) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (activeReservations_ <= 0 ||
            reservedBytes > activeReservedBytes_) {
            std::terminate();
        }
        const auto found = activeReservationInfo_.find(reservationId);
        if (found == activeReservationInfo_.end() ||
            found->second.bytes != reservedBytes) {
            std::terminate();
        }
        activeReservationInfo_.erase(found);
        activeReservedBytes_ -= reservedBytes;
        --activeReservations_;
        observedResidentBytesLocked();
    }
    memoryAvailable_.notify_all();
}

AlignmentMemorySchedulerStats AlignmentMemoryScheduler::stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stats_;
}

ScopedAlignmentMemoryScheduler::ScopedAlignmentMemoryScheduler(
        AlignmentMemoryScheduler &scheduler)
        : scheduler_(&scheduler) {
    AlignmentMemoryScheduler *expected = nullptr;
    if (!installedMemoryScheduler.compare_exchange_strong(
                expected, scheduler_, std::memory_order_release,
                std::memory_order_relaxed)) {
        scheduler_ = nullptr;
        throw std::logic_error(
                "an alignment memory scheduler is already installed");
    }
}

ScopedAlignmentMemoryScheduler::~ScopedAlignmentMemoryScheduler() {
    if (scheduler_ == nullptr) {
        return;
    }
    AlignmentMemoryScheduler *expected = scheduler_;
    if (!installedMemoryScheduler.compare_exchange_strong(
                expected, nullptr, std::memory_order_acq_rel,
                std::memory_order_relaxed)) {
        std::terminate();
    }
}

AlignmentMemoryReservation acquireAlignmentMemory(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    if (scheduler == nullptr) {
        return AlignmentMemoryReservation(
                nullptr, 0, predictedPeakBytes, true);
    }
    return scheduler->acquire(predictedPeakBytes, predictedRuntimeMinutes);
}

AlignmentMemoryTryResult tryAcquireAlignmentMemory(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    if (scheduler == nullptr) {
        AlignmentMemoryTryResult result;
        result.admission = AlignmentMemoryAdmission::Acquired;
        result.reservation = AlignmentMemoryReservation(
                nullptr, 0, predictedPeakBytes, true);
        return result;
    }
    return scheduler->tryAcquire(predictedPeakBytes,
                                 predictedRuntimeMinutes);
}

AlignmentMemoryTryResult acquirePreferredAlignmentMemory(
        uint64_t predictedPeakBytes, double predictedRuntimeMinutes) {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    if (scheduler == nullptr) {
        AlignmentMemoryTryResult result;
        result.admission = AlignmentMemoryAdmission::Acquired;
        result.reservation = AlignmentMemoryReservation(
                nullptr, 0, predictedPeakBytes, true);
        return result;
    }
    return scheduler->acquirePreferred(predictedPeakBytes,
                                       predictedRuntimeMinutes);
}

double estimatedAlignmentMemoryWaitMinutes(uint64_t predictedPeakBytes) {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    return scheduler == nullptr ? 0.0
                                : scheduler->estimatedWaitMinutes(
                                          predictedPeakBytes);
}

MemoryAvailabilityEstimate estimateAlignmentMemoryAvailability(
        uint64_t predictedPeakBytes,
        AlignmentMemoryPressureState *pressureState) {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    if (scheduler == nullptr) {
        MemoryAvailabilityEstimate result;
        result.availability = MemoryAvailability::Ready;
        result.waitMinutes = 0.0;
        if (pressureState != nullptr) {
            pressureState->reset();
        }
        return result;
    }
    return scheduler->memoryAvailability(predictedPeakBytes, pressureState);
}

uint64_t maximumAlignmentTaskReservationBytes() {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    return scheduler == nullptr ? 0
                                : scheduler->maximumSingleTaskReservationBytes();
}

AlignmentResourceSnapshot alignmentResourceSnapshot() {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    return scheduler == nullptr ? AlignmentResourceSnapshot{}
                                : scheduler->snapshot();
}

bool alignmentMemorySchedulingEnabled() {
    AlignmentMemoryScheduler *scheduler = installedMemoryScheduler.load(
            std::memory_order_acquire);
    return scheduler != nullptr && scheduler->enabled();
}

uint64_t alignmentTaskTransientMemoryBytes(
        uint64_t queryLength, uint64_t referenceLength) {
    const uint64_t sequenceBytes = saturatingAdd(queryLength,
                                                 referenceLength);
    const uint64_t scaled = sequenceBytes >
                                    std::numeric_limits<uint64_t>::max() /
                                            kTaskBytesPerInputBase
                            ? std::numeric_limits<uint64_t>::max()
                            : sequenceBytes * kTaskBytesPerInputBase;
    return saturatingAdd(kTaskTransientBase, scaled);
}

DynamicThreadLimiter::DynamicThreadLimiter(int maximumActiveThreads)
        : maximumActiveThreads_(maximumActiveThreads) {
    if (maximumActiveThreads <= 0) {
        throw std::invalid_argument(
                "maximum active thread count must be positive");
    }
}

void DynamicThreadLimiter::acquire() {
    std::unique_lock<std::mutex> lock(mutex_);
    condition_.wait(lock, [this]() {
        return activeThreads_ < maximumActiveThreads_;
    });
    ++activeThreads_;
    peakActiveThreads_ = std::max(peakActiveThreads_, activeThreads_);
}

void DynamicThreadLimiter::release() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (activeThreads_ <= 0) {
            throw std::logic_error("alignment worker slot released twice");
        }
        --activeThreads_;
    }
    condition_.notify_all();
}

void DynamicThreadLimiter::waitUntilIdle() {
    std::unique_lock<std::mutex> lock(mutex_);
    condition_.wait(lock, [this]() { return activeThreads_ == 0; });
}

int DynamicThreadLimiter::maximumActiveThreads() const {
    return maximumActiveThreads_;
}

int DynamicThreadLimiter::peakActiveThreads() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return peakActiveThreads_;
}

}  // namespace anchorwave
