#include "WfaAlignment.h"

extern "C" {
#include "../../WFA2-lib/wavefront/wfa.h"
}

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#if defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(__GLIBC__)
#include <malloc.h>
#endif

namespace anchorwave {
namespace {

constexpr uint64_t kBiWfaComponents = 3;
// A cached aligner is private to one worker and one WFA memory mode. With
// dozens of workers, the old 512-MiB ceiling allowed completed attempts to
// retain tens of GiB after their global reservation was released. Most tasks
// also use a different memory cap, which invalidates that cache on the next
// call. Keep only genuinely small reusable working sets.
constexpr uint64_t kRetainedMemoryCeiling = 64ULL * 1024ULL * 1024ULL;
constexpr int kHistoricalBiWfaLeafScore = 250;
constexpr int kMaximumBiWfaLeafScore = 16384;
// Reserve half of one BiWFA component and price a two-piece affine high-WFA
// leaf at 64 bytes/score^2. Combined, this is sqrt(componentBudget/128).
constexpr uint64_t kBiWfaLeafBytesPerScoreSquared = 128;
constexpr int kMaximumInnerWfaThreads = 128;

struct WfaConfig {
    int32_t mismatch = 0;
    int32_t gapOpen1 = 0;
    int32_t gapExtend1 = 0;
    int32_t gapOpen2 = 0;
    int32_t gapExtend2 = 0;
    uint64_t memoryLimitBytes = 0;
    uint64_t aggregateMemoryLimitBytes = 0;
    wavefront_memory_t memoryMode = wavefront_memory_high;
    bool singletrack = false;
    int maxNumThreads = 1;
    int minOffsetsPerThread = 500;
    int biWfaLeafScore = kHistoricalBiWfaLeafScore;
    int memoryProbeScoreInterval = 3000;
    WfaMemoryProbe memoryProbe = nullptr;
    void *memoryProbeArguments = nullptr;

    bool operator==(const WfaConfig &other) const {
        return mismatch == other.mismatch &&
               gapOpen1 == other.gapOpen1 &&
               gapExtend1 == other.gapExtend1 &&
               gapOpen2 == other.gapOpen2 &&
               gapExtend2 == other.gapExtend2 &&
               memoryLimitBytes == other.memoryLimitBytes &&
               aggregateMemoryLimitBytes == other.aggregateMemoryLimitBytes &&
               memoryMode == other.memoryMode &&
               singletrack == other.singletrack &&
               maxNumThreads == other.maxNumThreads &&
               minOffsetsPerThread == other.minOffsetsPerThread &&
               biWfaLeafScore == other.biWfaLeafScore &&
               memoryProbeScoreInterval == other.memoryProbeScoreInterval &&
               memoryProbe == other.memoryProbe;
    }
};

struct CompositeProbeContext {
    WfaMemoryProbe externalProbe = nullptr;
    void *externalArguments = nullptr;
    std::chrono::steady_clock::time_point deadline;
    bool deadlineEnabled = false;
    bool timeLimitExceeded = false;
};

bool compositeResourceProbe(void *arguments,
                            uint64_t memoryUsedBytes,
                            uint64_t memoryLimitBytes,
                            int score) {
    CompositeProbeContext *const context =
            static_cast<CompositeProbeContext *>(arguments);
    if (context == nullptr) {
        return true;
    }
    if (context->externalProbe != nullptr &&
        !context->externalProbe(context->externalArguments,
                                memoryUsedBytes, memoryLimitBytes, score)) {
        return false;
    }
    if (context->deadlineEnabled &&
        std::chrono::steady_clock::now() >= context->deadline) {
        context->timeLimitExceeded = true;
        return false;
    }
    return true;
}

int32_t positivePenalty(int32_t penalty) {
    if (penalty == std::numeric_limits<int32_t>::min()) {
        return std::numeric_limits<int32_t>::max();
    }
    return penalty < 0 ? -penalty : penalty;
}

int normalizedMaxThreads(int requested) {
    return std::max(1, std::min(requested, kMaximumInnerWfaThreads));
}

int normalizedMinOffsetsPerThread(int requested) {
    return std::max(1, requested);
}

int normalizedProbeInterval(int requested) {
    return std::max(1, requested);
}

class ThreadWfaAligner {
public:
    ThreadWfaAligner() = default;
    ThreadWfaAligner(const ThreadWfaAligner &) = delete;
    ThreadWfaAligner &operator=(const ThreadWfaAligner &) = delete;

    ~ThreadWfaAligner() {
        reset();
    }

    wavefront_aligner_t *get(const WfaConfig &config) {
        if (aligner_ != nullptr && config == config_) {
            aligner_->system.memory_probe = config.memoryProbe;
            aligner_->system.memory_probe_arguments =
                    config.memoryProbeArguments;
            aligner_->system.probe_interval_global =
                    config.memoryProbeScoreInterval;
            return aligner_;
        }
        reset();

        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine_2p;
        attributes.affine2p_penalties.match = 0;
        attributes.affine2p_penalties.mismatch = config.mismatch;
        attributes.affine2p_penalties.gap_opening1 = config.gapOpen1;
        attributes.affine2p_penalties.gap_extension1 = config.gapExtend1;
        attributes.affine2p_penalties.gap_opening2 = config.gapOpen2;
        attributes.affine2p_penalties.gap_extension2 = config.gapExtend2;
        attributes.alignment_form.span = alignment_end2end;
        attributes.alignment_scope = compute_alignment;
        attributes.memory_mode = config.memoryMode;
        attributes.singletrack = config.singletrack;
        attributes.heuristic.strategy = wf_heuristic_none;
        attributes.system.max_memory_abort = config.memoryLimitBytes;
        attributes.system.biwfa_max_memory_abort =
                config.aggregateMemoryLimitBytes;
        attributes.system.biwfa_leaf_score = config.biWfaLeafScore;
        attributes.system.probe_interval_global =
                config.memoryProbeScoreInterval;
        attributes.system.memory_probe = config.memoryProbe;
        attributes.system.memory_probe_arguments =
                config.memoryProbeArguments;
        attributes.system.verbose = 0;
        attributes.system.max_num_threads = config.maxNumThreads;
        attributes.system.min_offsets_per_thread =
                config.minOffsetsPerThread;

        aligner_ = wavefront_aligner_new(&attributes);
        config_ = config;
        return aligner_;
    }

    void reapIfLarge(uint64_t aggregateBudget, uint64_t memoryUsed) {
        if (aligner_ == nullptr) {
            return;
        }
        const uint64_t retainedTarget = std::min(aggregateBudget / 4,
                                                 kRetainedMemoryCeiling);
        if (memoryUsed > retainedTarget) {
            // WFA2's reap operation returns slabs to the aligner's private
            // allocator, but that allocator deliberately keeps the pages for
            // reuse.  Across several worker threads and memory modes those
            // inactive pages remain in process RSS after the scheduler token
            // is released, so a later large exact request can wait forever
            // even though no alignment is active.  A large working set is
            // unlikely to be reusable by the common short gaps; destroy the
            // complete aligner so its private allocator is returned to the
            // system and keep caching only genuinely small working sets.
            reset();
        }
    }

    void reap() {
        // Failure paths immediately try another exact engine.  Retaining the
        // failed engine's private allocator would overlap it with the
        // fallback reservation and invalidate the global -M accounting.
        reset();
    }

private:
    void reset() {
        if (aligner_ != nullptr) {
            wavefront_aligner_delete(aligner_);
            aligner_ = nullptr;
        }
    }

    wavefront_aligner_t *aligner_ = nullptr;
    WfaConfig config_;
};

// Keep independent per-thread caches. Switching one cache among memory modes
// would destroy and recreate aligners at every fallback, negating much of the
// benefit of the medium/low modes on workloads with many short intervals.
thread_local ThreadWfaAligner standardThreadAligner;
thread_local ThreadWfaAligner singletrackThreadAligner;
thread_local ThreadWfaAligner mediumThreadAligner;
thread_local ThreadWfaAligner lowThreadAligner;
thread_local ThreadWfaAligner bidirectionalThreadAligner;

bool decodeCigar(const std::string &query,
                 const std::string &reference,
                 const cigar_t &cigar,
                 std::string &queryAlignment,
                 std::string &referenceAlignment) {
    queryAlignment.clear();
    referenceAlignment.clear();
    queryAlignment.reserve(query.size() + reference.size());
    referenceAlignment.reserve(query.size() + reference.size());

    std::size_t referencePosition = 0;
    std::size_t queryPosition = 0;
    for (int i = cigar.begin_offset; i < cigar.end_offset; ++i) {
        switch (cigar.operations[i]) {
            case 'M':
            case 'X':
                if (queryPosition >= query.size() ||
                    referencePosition >= reference.size()) {
                    return false;
                }
                queryAlignment.push_back(query[queryPosition++]);
                referenceAlignment.push_back(reference[referencePosition++]);
                break;
            case 'I':
                if (queryPosition >= query.size()) {
                    return false;
                }
                queryAlignment.push_back(query[queryPosition++]);
                referenceAlignment.push_back('-');
                break;
            case 'D':
                if (referencePosition >= reference.size()) {
                    return false;
                }
                queryAlignment.push_back('-');
                referenceAlignment.push_back(reference[referencePosition++]);
                break;
            default:
                return false;
        }
    }

    return queryPosition == query.size() &&
           referencePosition == reference.size() &&
           queryAlignment.size() == referenceAlignment.size();
}

WfaAlignmentStatus alignWithWfaMode(
        const std::string &query,
        const std::string &reference,
        uint64_t aggregateMemoryBudgetBytes,
        uint64_t alignerMemoryLimitBytes,
        wavefront_memory_t memoryMode,
        WfaAlgorithm algorithm,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        ThreadWfaAligner &threadAligner,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    result = WfaAlignmentResult{};
    result.memoryBudgetBytes = aggregateMemoryBudgetBytes;
    result.algorithm = algorithm;
    result.configuredMaxThreads = normalizedMaxThreads(
            options.maxNumThreads);
    if (aggregateMemoryBudgetBytes == 0 || alignerMemoryLimitBytes == 0 ||
        query.size() > INT_MAX || reference.size() > INT_MAX) {
        result.wfaStatus = WF_STATUS_UNATTAINABLE;
        return WfaAlignmentStatus::Failed;
    }

    WfaConfig config;
    config.mismatch = positivePenalty(mismatchingPenalty);
    config.gapOpen1 = positivePenalty(openGapPenalty1);
    config.gapExtend1 = positivePenalty(extendGapPenalty1);
    config.gapOpen2 = positivePenalty(openGapPenalty2);
    config.gapExtend2 = positivePenalty(extendGapPenalty2);
    config.memoryLimitBytes = alignerMemoryLimitBytes;
    config.aggregateMemoryLimitBytes = aggregateMemoryBudgetBytes;
    config.memoryMode = memoryMode;
    config.singletrack = algorithm == WfaAlgorithm::Singletrack;
    config.maxNumThreads = result.configuredMaxThreads;
    config.minOffsetsPerThread = normalizedMinOffsetsPerThread(
            options.minOffsetsPerThread);
    config.biWfaLeafScore = options.biWfaLeafScore > 0
            ? std::min(options.biWfaLeafScore,
                       kMaximumBiWfaLeafScore)
            : biWfaLeafScoreFromMemoryBudgetBytes(
                    aggregateMemoryBudgetBytes);
    config.memoryProbeScoreInterval = normalizedProbeInterval(
            options.memoryProbeScoreInterval);
    CompositeProbeContext probeContext;
    probeContext.externalProbe = options.memoryProbe;
    probeContext.externalArguments = options.memoryProbeArguments;
    if (options.maximumRuntimeMilliseconds > 0) {
        probeContext.deadlineEnabled = true;
        probeContext.deadline = std::chrono::steady_clock::now() +
                std::chrono::milliseconds(
                        options.maximumRuntimeMilliseconds);
        // Limit wall-clock overshoot without materially changing the WFA
        // recurrence cost. BiWFA subsidiaries also share this aggregate
        // callback and deadline.
        config.memoryProbeScoreInterval = std::min(
                config.memoryProbeScoreInterval, 256);
        config.memoryProbe = &compositeResourceProbe;
        config.memoryProbeArguments = &probeContext;
    } else {
        config.memoryProbe = options.memoryProbe;
        config.memoryProbeArguments = options.memoryProbeArguments;
    }
    result.configuredBiWfaLeafScore =
            algorithm == WfaAlgorithm::Bidirectional
            ? config.biWfaLeafScore : 0;
    wavefront_aligner_t *const aligner = threadAligner.get(config);
    if (aligner == nullptr) {
        result.wfaStatus = WF_STATUS_OOM;
        return WfaAlignmentStatus::MemoryLimit;
    }

    result.wfaStatus = wavefront_align(
            aligner,
            reference.data(), static_cast<int>(reference.size()),
            query.data(), static_cast<int>(query.size()));
    // The callback argument can point to this stack frame. Cached aligners
    // must never retain it between attempts; get() installs the next call's
    // callback and argument again before wavefront_align().
    aligner->system.memory_probe = nullptr;
    aligner->system.memory_probe_arguments = nullptr;
    result.memoryUsedBytes = aligner->align_status.memory_used;
    result.memoryPeakBytes = aligner->align_status.memory_peak;
    result.memoryProbeCount = aligner->align_status.memory_probe_count;
    result.maxThreadsUsed = aligner->align_status.max_num_threads_used;
    result.timeLimitExceeded = probeContext.timeLimitExceeded;
    if (aligner->bialigner != nullptr) {
        result.biWfaLeafAlignments =
                aligner->bialigner->num_leaf_alignments;
        result.biWfaMaxLeafScore = aligner->bialigner->max_leaf_score;
        result.maxThreadsUsed =
                aligner->bialigner->max_num_threads_used;
    }

    if (result.wfaStatus != WF_STATUS_ALG_COMPLETED) {
        // A fallback can start immediately after this call. Do not allow the
        // failed standard-WFA working set to overlap the BiWFA allocation.
        threadAligner.reap();
        if (result.timeLimitExceeded) {
            return WfaAlignmentStatus::TimeLimit;
        }
        return result.wfaStatus == WF_STATUS_OOM
               ? WfaAlignmentStatus::MemoryLimit
               : WfaAlignmentStatus::Failed;
    }

    result.score = aligner->cigar->score;
    if (!decodeCigar(query, reference, *aligner->cigar,
                     result.queryAlignment, result.referenceAlignment)) {
        result.wfaStatus = WF_STATUS_UNATTAINABLE;
        result.queryAlignment.clear();
        result.referenceAlignment.clear();
        threadAligner.reap();
        return WfaAlignmentStatus::Failed;
    }

    threadAligner.reapIfLarge(aggregateMemoryBudgetBytes,
                              result.memoryUsedBytes);
    return WfaAlignmentStatus::Completed;
}

}  // namespace

uint64_t wfaMemoryBudgetBytes(int64_t windowWidth) {
    if (windowWidth <= 0) {
        return 0;
    }
    const uint64_t width = static_cast<uint64_t>(windowWidth);
    if (width > std::numeric_limits<uint64_t>::max() / width) {
        return std::numeric_limits<uint64_t>::max();
    }
    return width * width;
}

uint64_t standardWfaTrialMemoryBudgetBytes(
        uint64_t workerMemoryBudgetBytes) {
    return workerMemoryBudgetBytes;
}

void reapCurrentThreadWfaCaches() {
    standardThreadAligner.reap();
    singletrackThreadAligner.reap();
    mediumThreadAligner.reap();
    lowThreadAligner.reap();
    bidirectionalThreadAligner.reap();
}

bool releaseUnusedAlignmentMemoryToSystem() {
    reapCurrentThreadWfaCaches();
#if defined(__APPLE__)
    return malloc_zone_pressure_relief(nullptr, 0) != 0;
#elif defined(__GLIBC__)
    return malloc_trim(0) != 0;
#else
    return false;
#endif
}

uint64_t biWfaComponentMemoryLimitBytes(
        uint64_t workerMemoryBudgetBytes) {
    if (workerMemoryBudgetBytes == 0) {
        return 0;
    }
    return std::max<uint64_t>(
            1, workerMemoryBudgetBytes / kBiWfaComponents);
}

int biWfaLeafScoreFromMemoryBudgetBytes(
        uint64_t workerMemoryBudgetBytes) {
    const uint64_t componentBytes =
            biWfaComponentMemoryLimitBytes(workerMemoryBudgetBytes);
    if (componentBytes == 0) {
        return kHistoricalBiWfaLeafScore;
    }
    const long double scaled = static_cast<long double>(componentBytes) /
            static_cast<long double>(kBiWfaLeafBytesPerScoreSquared);
    int score = static_cast<int>(std::sqrt(scaled));
    score = std::max(score, kHistoricalBiWfaLeafScore);
    return std::min(score, kMaximumBiWfaLeafScore);
}

bool wfaParallelSupportEnabled() {
#ifdef WFA_PARALLEL
    return true;
#else
    return false;
#endif
}

WfaAlignmentStatus alignWithSingletrackWfa(
        const std::string &query,
        const std::string &reference,
        uint64_t memoryBudgetBytes,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    return alignWithWfaMode(
            query, reference,
            memoryBudgetBytes, memoryBudgetBytes,
            wavefront_memory_high, WfaAlgorithm::Singletrack,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            singletrackThreadAligner, result, options);
}

WfaAlignmentStatus alignWithStandardWfa(
        const std::string &query,
        const std::string &reference,
        uint64_t memoryBudgetBytes,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    return alignWithWfaMode(
            query, reference,
            memoryBudgetBytes, memoryBudgetBytes,
            wavefront_memory_high, WfaAlgorithm::Standard,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            standardThreadAligner, result, options);
}

WfaAlignmentStatus alignWithMediumWfa(
        const std::string &query,
        const std::string &reference,
        uint64_t memoryBudgetBytes,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    return alignWithWfaMode(
            query, reference,
            memoryBudgetBytes, memoryBudgetBytes,
            wavefront_memory_med, WfaAlgorithm::Medium,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            mediumThreadAligner, result, options);
}

WfaAlignmentStatus alignWithLowWfa(
        const std::string &query,
        const std::string &reference,
        uint64_t memoryBudgetBytes,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    return alignWithWfaMode(
            query, reference,
            memoryBudgetBytes, memoryBudgetBytes,
            wavefront_memory_low, WfaAlgorithm::Low,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            lowThreadAligner, result, options);
}

WfaAlignmentStatus alignWithBiWfa(
        const std::string &query,
        const std::string &reference,
        uint64_t memoryBudgetBytes,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    return alignWithWfaMode(
            query, reference,
            memoryBudgetBytes,
            biWfaComponentMemoryLimitBytes(memoryBudgetBytes),
            wavefront_memory_ultralow, WfaAlgorithm::Bidirectional,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            bidirectionalThreadAligner, result, options);
}

WfaAlignmentStatus alignWithAdaptiveWfa(
        const std::string &query,
        const std::string &reference,
        uint64_t workerMemoryBudgetBytes,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        WfaAlignmentResult &result,
        const WfaExecutionOptions &options) {
    const uint64_t standardBudget =
            standardWfaTrialMemoryBudgetBytes(workerMemoryBudgetBytes);
    WfaAlignmentStatus status = alignWithSingletrackWfa(
            query, reference, workerMemoryBudgetBytes,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            result, options);
    if (status == WfaAlignmentStatus::Completed) {
        return status;
    }
    if (status == WfaAlignmentStatus::TimeLimit) {
        return status;
    }

    status = alignWithStandardWfa(
            query, reference, standardBudget,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            result, options);
    if (status == WfaAlignmentStatus::Completed) {
        return status;
    }
    if (status == WfaAlignmentStatus::TimeLimit) {
        return status;
    }

    status = alignWithMediumWfa(
            query, reference, workerMemoryBudgetBytes,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            result, options);
    if (status == WfaAlignmentStatus::Completed) {
        return status;
    }
    if (status == WfaAlignmentStatus::TimeLimit) {
        return status;
    }

    status = alignWithLowWfa(
            query, reference, workerMemoryBudgetBytes,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            result, options);
    if (status == WfaAlignmentStatus::Completed) {
        return status;
    }
    if (status == WfaAlignmentStatus::TimeLimit) {
        return status;
    }

    return alignWithBiWfa(
            query, reference, workerMemoryBudgetBytes,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            result, options);
}

}  // namespace anchorwave
