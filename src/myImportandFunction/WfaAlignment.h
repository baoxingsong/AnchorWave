#pragma once

#include <cstdint>
#include <string>

namespace anchorwave {

enum class WfaAlignmentStatus {
    Completed,
    MemoryLimit,
    TimeLimit,
    Failed
};

enum class WfaAlgorithm {
    None,
    Singletrack,
    Standard,
    Medium,
    Low,
    Bidirectional
};

// Called at WFA score-probe boundaries. Returning false stops the current
// exact attempt with MemoryLimit. This is a cooperative growth hook, not an
// allocator-level hard cap: WFA can allocate between consecutive probes.
using WfaMemoryProbe = bool (*)(void *arguments,
                               uint64_t memoryUsedBytes,
                               uint64_t memoryLimitBytes,
                               int score);

struct WfaExecutionOptions {
    // Keep one inner WFA thread during ordinary genome-level parallelism.
    // A tail scheduler may explicitly assign 2/4/8 to a long final task.
    int maxNumThreads = 1;
    int minOffsetsPerThread = 500;

    // Zero derives a conservative exact high-WFA leaf threshold from this
    // attempt's aggregate memory budget. 250 reproduces the historical WFA2
    // setting exactly.
    int biWfaLeafScore = 0;

    int memoryProbeScoreInterval = 3000;
    WfaMemoryProbe memoryProbe = nullptr;
    void *memoryProbeArguments = nullptr;

    // Cooperative wall-clock deadline for this exact attempt. Zero disables
    // the runtime limit. The same WFA score probe used for memory admission
    // enforces this without adding a polling thread.
    uint64_t maximumRuntimeMilliseconds = 0;
};

struct WfaAlignmentResult {
    int64_t score = 0;
    std::string queryAlignment;
    std::string referenceAlignment;
    uint64_t memoryBudgetBytes = 0;
    uint64_t memoryUsedBytes = 0;
    uint64_t memoryPeakBytes = 0;
    uint64_t memoryProbeCount = 0;
    bool timeLimitExceeded = false;
    int wfaStatus = 0;
    int configuredMaxThreads = 1;
    int maxThreadsUsed = 1;
    int configuredBiWfaLeafScore = 0;
    uint64_t biWfaLeafAlignments = 0;
    int biWfaMaxLeafScore = 0;
    WfaAlgorithm algorithm = WfaAlgorithm::None;
};

// Preserve the historical -w resource model: a width of w represents an
// approximate per-AnchorWave-thread alignment budget of w^2 bytes.
uint64_t wfaMemoryBudgetBytes(int64_t windowWidth);

// Standard and Singletrack full-CIGAR WFA are high-memory fast paths.
// Singletrack retains only M wavefronts for traceback. Both modes use the
// caller's complete per-alignment -w^2 budget; process-wide -M scheduling
// controls how many such attempts may run concurrently.
uint64_t standardWfaTrialMemoryBudgetBytes(uint64_t workerMemoryBudgetBytes);

// Release every WFA cache owned by the calling worker. The resource scheduler
// uses this once before treating a stable RSS floor as unavailable memory.
void reapCurrentThreadWfaCaches();

// Return released alignment allocations to the operating system when the
// process is memory-idle but its allocator still reports a high RSS. This is
// a pressure-path operation, not part of the normal per-alignment fast path.
// Returns true when the platform allocator reports that it released pages.
bool releaseUnusedAlignmentMemoryToSystem();

// BiWFA owns three subsidiary WFA objects. Divide an already-computed byte
// budget among them; no sequence length, score, or -w admission rule is used.
uint64_t biWfaComponentMemoryLimitBytes(uint64_t workerMemoryBudgetBytes);

// Conservative execution rule for exact high-WFA leaves inside BiWFA. The
// result is monotonic in the aggregate attempt budget, never below the legacy
// value 250, and capped to avoid replacing BiWFA with one giant high WFA.
int biWfaLeafScoreFromMemoryBudgetBytes(uint64_t workerMemoryBudgetBytes);

// True only when this binary compiled WFA2 with its OpenMP kernels. A thread
// request above one remains advisory when this returns false.
bool wfaParallelSupportEnabled();

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
        const WfaExecutionOptions &options = WfaExecutionOptions());

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
        const WfaExecutionOptions &options = WfaExecutionOptions());

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
        const WfaExecutionOptions &options = WfaExecutionOptions());

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
        const WfaExecutionOptions &options = WfaExecutionOptions());

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
        const WfaExecutionOptions &options = WfaExecutionOptions());

// Try Singletrack and the exact WFA memory modes from fastest/highest-memory
// to progressively lower-memory, then use BiWFA. Every failed working set is
// released before the next mode starts and all modes remain exact.
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
        const WfaExecutionOptions &options = WfaExecutionOptions());

}  // namespace anchorwave
