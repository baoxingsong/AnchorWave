//
// Created by song on 8/5/18.
//

#ifndef HAVE_KALLOC
#define HAVE_KALLOC 1
#endif
#include "../../minimap2/kalloc.h"
#include "alignSlidingWindow.h"
#include "AlignmentAlgorithmSelector.h"
#include "AlignmentResourceScheduler.h"
#include "WfaAlignment.h"
#include "src/service/AnchorTaskExecutor.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

struct CandidateMemoryEnvelope {
    uint64_t algorithmBudgetBytes = 0;
    uint64_t reservationBytes = 0;
};

uint64_t saturatedAdd(uint64_t first, uint64_t second) {
    if (second > std::numeric_limits<uint64_t>::max() - first) {
        return std::numeric_limits<uint64_t>::max();
    }
    return first + second;
}

int64_t scoreGappedTwoPieceAffineAlignment(
        const std::string &queryAlignment,
        const std::string &referenceAlignment,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2);

uint64_t addFractionalMargin(uint64_t bytes, uint64_t numerator,
                             uint64_t denominator) {
    if (bytes == 0) {
        return 1;
    }
    if (bytes > std::numeric_limits<uint64_t>::max() / numerator) {
        return std::numeric_limits<uint64_t>::max();
    }
    return saturatedAdd(bytes * numerator / denominator,
                        bytes * numerator % denominator != 0 ? 1 : 0);
}

const anchorwave::AlgorithmCostEstimate *candidateEstimate(
        const anchorwave::AlignmentSelectionPlan &plan,
        anchorwave::AlignmentCandidate candidate) {
    for (const anchorwave::AlgorithmCostEstimate &estimate : plan.estimates) {
        if (estimate.candidate == candidate) {
            return &estimate;
        }
    }
    return nullptr;
}

bool isWfaCandidate(anchorwave::AlignmentCandidate candidate) {
    return candidate == anchorwave::AlignmentCandidate::SingletrackWfa ||
           candidate == anchorwave::AlignmentCandidate::StandardWfa ||
           candidate == anchorwave::AlignmentCandidate::MediumWfa ||
           candidate == anchorwave::AlignmentCandidate::LowWfa;
}

CandidateMemoryEnvelope candidateMemoryEnvelope(
        const anchorwave::AlignmentSelectionPlan &plan,
        anchorwave::AlignmentCandidate candidate,
        uint64_t workerBudget, uint64_t transientBytes,
        uint64_t predictedAlgorithmOverride = 0) {
    const anchorwave::AlgorithmCostEstimate *estimate = candidateEstimate(
            plan, candidate);
    uint64_t predictedAlgorithmBytes = estimate != nullptr
                                       ? estimate->memoryBytes
                                       : workerBudget;
    if (predictedAlgorithmOverride > 0) {
        predictedAlgorithmBytes = predictedAlgorithmOverride;
    }
    // B73/Mo17 telemetry found that WFA's allocator peak can exceed the
    // analytic working-set prediction by about 20%, and an abort followed by
    // medium WFA added ~20 seconds to the critical path. Give WFA 25%
    // algorithm headroom plus a 15% reservation guard for allocator/probe
    // overshoot. The AVX512 page-aware KSW2 model predicted isolated RSS
    // within 1%, so do not duplicate the process-wide safety reserve with a
    // per-attempt percentage there. Transient
    // input/output bytes are still charged explicitly, and the adaptive
    // process reserve remains inside -M for aggregate allocator residuals.
    const uint64_t predictionMarginPercent =
            isWfaCandidate(candidate) ? 125 : 100;
    uint64_t algorithmBudget = addFractionalMargin(
            predictedAlgorithmBytes, predictionMarginPercent, 100);
    // The engine working set is capped at w^2 except for the tier-wide
    // qlen*rlen <= w^2 quality contract. The reservation can be larger because
    // transient input/output and allocator guard bytes are separately charged
    // to process -M.
    uint64_t algorithmCeiling = workerBudget;
    const bool exactCandidate = std::find(
            plan.exactCandidates.begin(), plan.exactCandidates.end(),
            candidate) != plan.exactCandidates.end();
    if (plan.exactCellEnvelopeException && exactCandidate) {
        // The raw cell envelope is the quality exception, not a KSW2-full
        // exception. Let every Tier-1 implementation request its guarded
        // prediction; process-wide -M admission remains the hard bound.
        algorithmCeiling = std::max(workerBudget, algorithmBudget);
    }
    algorithmBudget = std::min(algorithmBudget, algorithmCeiling);
    const uint64_t reservationGuardPercent =
            isWfaCandidate(candidate) ? 115 : 100;
    uint64_t guardedAlgorithmBytes = addFractionalMargin(
            algorithmBudget, reservationGuardPercent, 100);
    return CandidateMemoryEnvelope{
            algorithmBudget,
            saturatedAdd(guardedAlgorithmBytes, transientBytes)};
}

std::size_t candidateIndex(anchorwave::AlignmentCandidate candidate) {
    return static_cast<std::size_t>(candidate);
}

struct PreparedAlignmentSelection {
    anchorwave::AlignmentSelectionPlan plan;
    const std::string *queryObject = nullptr;
    const std::string *referenceObject = nullptr;
    uint64_t intervalId = 0;
    uint64_t attemptCount = 0;
    uint64_t queryLength = 0;
    uint64_t referenceLength = 0;
    double exactRuntimeSpentSeconds = 0.0;
    double certifiedScoreRuntimeSeconds = 0.0;
    int64_t certifiedOptimalScore = 0;
    bool certifiedOptimalScoreReady = false;
    int64_t certifiedFinalBandWidth = -1;
    uint64_t certifiedTracebackAttempts = 0;
    bool exactTimeLimitReached = false;
    bool lastExactAttemptMemoryLimited = false;
    uint64_t lastExactAttemptMemoryBytes = 0;
    bool singletrackMemoryModelInvalidated = false;
    bool singletrackExpandedRetryPending = false;
    bool singletrackExpandedRetryUsed = false;
    bool selectionRecorded = false;
    bool planTraced = false;
    int64_t retainedBandedScore = std::numeric_limits<int64_t>::min();
    int64_t retainedBandedBandWidth = -1;
    std::string retainedBandedQueryAlignment;
    std::string retainedBandedReferenceAlignment;
    std::array<bool, anchorwave::kAlignmentCandidateCount> attempted{};
    std::array<bool, anchorwave::kAlignmentCandidateCount>
            runtimeFloorRecorded{};
    std::array<bool, anchorwave::kAlignmentCandidateCount>
            allocatorReclaimAttempted{};
    std::array<anchorwave::AlignmentMemoryPressureState,
               anchorwave::kAlignmentCandidateCount> memoryPressure{};
};

struct RuntimeCandidate {
    anchorwave::AlignmentCandidate candidate =
            anchorwave::AlignmentCandidate::SlidingWindow;
    const anchorwave::AlgorithmCostEstimate *estimate = nullptr;
    anchorwave::AlgorithmCostEstimate effectiveEstimate;
    CandidateMemoryEnvelope memoryEnvelope;
    double predictedWaitMinutes = 0.0;
    double expectedCompletionMinutes = 0.0;
    double memoryOpportunityMinutes = 0.0;
    double systemCost = 0.0;
    bool expandedSingletrackRetry = false;
    anchorwave::MemoryAvailabilityEstimate memoryAvailability;
};

double dynamicMemoryShadowMinutes(
        double runtimeMinutes,
        uint64_t reservationBytes,
        const anchorwave::AlignmentResourceSnapshot &resources,
        const anchorwave::AnchorTaskExecutor::LoadSnapshot &load) {
    if (load.admissionTailPhase || !(runtimeMinutes > 0.0)) {
        return 0.0;
    }
    const long double capacity = resources.enabled
            ? static_cast<long double>(std::max<uint64_t>(
                      1, resources.taskMemoryCapacityBytes))
            : static_cast<long double>(std::max<uint64_t>(
                      1, reservationBytes)) *
                      static_cast<long double>(std::max(1, load.workerCount));
    const int workerCount = std::max(1, load.workerCount);
    if (!resources.enabled || workerCount <= 1) {
        return 0.0;
    }
    // Convert memory into the number of average worker shares it occupies.
    // One share belongs to the task's own CPU slot; only the excess can displace
    // other runnable exact tasks. The previous formula used only the process
    // memory fraction, underpricing a 28-GiB task by roughly workerCount-fold on
    // a 34-worker/102-GiB run.
    const double memoryEquivalentWorkers = static_cast<double>(
            static_cast<long double>(reservationBytes) * workerCount /
            capacity);
    const double displacedWorkers = std::max(
            0.0, memoryEquivalentWorkers - 1.0);
    const double competingTasks = load.schedulableTasks > 0
            ? static_cast<double>(load.schedulableTasks - 1) : 0.0;
    const double backlogUtilization = std::min(
            1.0, competingTasks /
                         static_cast<double>(workerCount - 1));
    return runtimeMinutes * displacedWorkers * backlogUtilization;
}

}  // namespace

namespace {

class ScopedKswScratchArena {
public:
    ScopedKswScratchArena() : arena_(km_init()) {}
    ~ScopedKswScratchArena() { km_destroy(arena_); }

    ScopedKswScratchArena(const ScopedKswScratchArena &) = delete;
    ScopedKswScratchArena &operator=(const ScopedKswScratchArena &) = delete;

    void *get() const { return arena_; }
    void reset() {
        km_destroy(arena_);
        arena_ = km_init();
    }

private:
    void *arena_;
};

int32_t minimap2AlignmentWithArena(void *arena,
                           const std::string &_dna_q, const std::string &_dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int32_t &matchingScore, int32_t mismatchingPenalty,
                           int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2) {
    int8_t a = 0, b = mismatchingPenalty < 0 ? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    // Keep KSW2 on the same objective as WFA.  The fifth symbol represents
    // every non-ACGT base.  Treating it as a zero-cost wildcard (the old
    // matrix) made KSW2 and WFA return different "exact" scores on N-rich
    // genome intervals.  Equal ambiguous symbols match; every unequal pair
    // receives the configured mismatch penalty.
    int8_t mat[25] = {a, b, b, b, b,
                      b, a, b, b, b,
                      b, b, a, b, b,
                      b, b, b, a, b,
                      b, b, b, b, a};
    int tl = _dna_d.length(), ql = _dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = _dna_d.c_str();
    const char *qseq = _dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);

    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table

    ts = (uint8_t *) malloc(tl);
    qs = (uint8_t *) malloc(ql);

    int i;
    int kswFlags = 0;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
        if (ts[i] == 4) kswFlags |= KSW_EZ_GENERIC_SC;
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t) qseq[i]];
        if (qs[i] == 4) kswFlags |= KSW_EZ_GENERIC_SC;
    }

    //flag  0x01 score only

#ifdef __AVX512BW__
    //    std::cout << "using AVX512" << std::endl;
        ksw_extd2_avx512(arena, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, kswFlags, & ez);
#elif __AVX2__
    //    std::cout << "using AVX2" << std::endl;
    ksw_extd2_avx2(arena, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, kswFlags, & ez);
#else
    //    std::cout << "using SSE" << std::endl;
    ksw_extd2_sse(arena, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, kswFlags, &ez);
#endif
    _alignment_q.clear();
    _alignment_d.clear();
    _alignment_q.reserve(_dna_q.size() + _dna_d.size());
    _alignment_d.reserve(_alignment_q.capacity());
    int pattern_pos = 0, text_pos = 0;

    // Consume KSW2's packed CIGAR directly. The former text conversion,
    // reparsing and single-base append loop produced identical output with
    // substantially more allocation and string work.
    for (i = 0; i < ez.n_cigar; ++i) {
        const int cLen = static_cast<int>(ez.cigar[i] >> 4);
        const int operation = static_cast<int>(ez.cigar[i] & 0xf);
        if (operation == 0 || operation == 7 || operation == 8) {
            _alignment_q.append(_dna_q, text_pos, cLen);
            _alignment_d.append(_dna_d, pattern_pos, cLen);
            pattern_pos += cLen;
            text_pos += cLen;
        } else if (operation == 1) {
            _alignment_q.append(_dna_q, text_pos, cLen);
            _alignment_d.append(static_cast<std::size_t>(cLen), '-');
            text_pos += cLen;
        } else if (operation == 2) {
            _alignment_q.append(static_cast<std::size_t>(cLen), '-');
            _alignment_d.append(_dna_d, pattern_pos, cLen);
            pattern_pos += cLen;
        }
    }
    assert(text_pos == ql);
    assert(pattern_pos == tl);
    int32_t totalScore = ez.score;
    kfree(arena, ez.cigar);
    free(ts);
    free(qs);
    return totalScore;
}

int32_t alignmentMinimap2WithArena(void *arena,
                           const std::string &_dna_q, const std::string &_dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int32_t &matchingScore, int32_t mismatchingPenalty, int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2, int32_t &endPositionq, int32_t &endPositiont) {

    int8_t a = 0, b = mismatchingPenalty < 0 ? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, b,
                      b, a, b, b, b,
                      b, b, a, b, b,
                      b, b, b, a, b,
                      b, b, b, b, a};
    int tl = _dna_d.length(), ql = _dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = _dna_d.c_str();
    const char *qseq = _dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);

    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table

    ts = (uint8_t *) malloc(tl);
    qs = (uint8_t *) malloc(ql);

    int i;
    // The historical implementation first ran score-only KSW2 to decide
    // which sequence end was better, then recomputed the same DP matrix with
    // traceback on the selected prefixes. Ask KSW2 to retain traceback once
    // and backtrack directly from the better end instead.
    int kswFlags = KSW_EZ_SEMIGLOBAL_END;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
        if (ts[i] == 4) kswFlags |= KSW_EZ_GENERIC_SC;
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t) qseq[i]];
        if (qs[i] == 4) kswFlags |= KSW_EZ_GENERIC_SC;
    }

//flag  0x01 score only

#ifdef __AVX512BW__
    //    std::cout << "using AVX512" << std::endl;
    ksw_extd2_avx512(arena, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, kswFlags, & ez);
#elif __AVX2__
    //    std::cout << "using AVX2" << std::endl;
        ksw_extd2_avx2(arena, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, kswFlags, & ez);
#else
    //    std::cout << "using SSE" << std::endl;
    ksw_extd2_sse(arena, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, kswFlags, &ez);
#endif

    int32_t maxScore = ez.score;
    if (ez.mqe > ez.mte) {
        endPositiont = ez.mqe_t + 1;
        endPositionq = ql;
    } else {
        endPositiont = tl;
        endPositionq = ez.mte_q + 1;
    }

    _alignment_q.clear();
    _alignment_d.clear();
    _alignment_q.reserve(static_cast<std::size_t>(endPositionq) +
                         static_cast<std::size_t>(endPositiont));
    _alignment_d.reserve(_alignment_q.capacity());
    int queryPosition = 0;
    int targetPosition = 0;
    for (i = 0; i < ez.n_cigar; ++i) {
        const int length = static_cast<int>(ez.cigar[i] >> 4);
        const int operation = static_cast<int>(ez.cigar[i] & 0xf);
        if (operation == 0 || operation == 7 || operation == 8) {
            _alignment_q.append(_dna_q, queryPosition, length);
            _alignment_d.append(_dna_d, targetPosition, length);
            queryPosition += length;
            targetPosition += length;
        } else if (operation == 1) {
            _alignment_q.append(_dna_q, queryPosition, length);
            _alignment_d.append(static_cast<std::size_t>(length), '-');
            queryPosition += length;
        } else if (operation == 2) {
            _alignment_q.append(static_cast<std::size_t>(length), '-');
            _alignment_d.append(_dna_d, targetPosition, length);
            targetPosition += length;
        }
    }
    assert(queryPosition == endPositionq);
    assert(targetPosition == endPositiont);
    kfree(arena, ez.cigar);
    free(ts);
    free(qs);
    return maxScore;
}

}  // namespace

int32_t minimap2_alignment(const std::string &_dna_q, const std::string &_dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int32_t &matchingScore, int32_t mismatchingPenalty,
                           int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2) {
    return minimap2AlignmentWithArena(
            nullptr, _dna_q, _dna_d, _alignment_q, _alignment_d,
            matchingScore, mismatchingPenalty,
            _open_gap_penalty1, _extend_gap_penalty1,
            _open_gap_penalty2, _extend_gap_penalty2);
}

int32_t alignment_minimap2(const std::string &_dna_q, const std::string &_dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int32_t &matchingScore, int32_t mismatchingPenalty, int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2, int32_t &endPositionq, int32_t &endPositiont) {
    return alignmentMinimap2WithArena(
            nullptr, _dna_q, _dna_d, _alignment_q, _alignment_d,
            matchingScore, mismatchingPenalty,
            _open_gap_penalty1, _extend_gap_penalty1,
            _open_gap_penalty2, _extend_gap_penalty2,
            endPositionq, endPositiont);
}

int64_t alignSlidingWindow(const std::string &dna_q, const std::string &dna_d, int64_t _length_of_q, int64_t _length_of_d,
                           std::string &_alignment_q, std::string &_alignment_d, const int64_t &slidingWindowSize,
                           const int32_t &matchingScore, const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                           const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    //2^15 = 32768
    //of the maximum length of the windowSize of is about 32000/2 = 16000
    int32_t databaseStart = 1;
    int32_t databaseEnd = 0;
    int32_t queryStart = 1;
    int32_t queryEnd = 0;
    int64_t totalScore = 0;
    int32_t endPositionq;
    int32_t endPositiont;
    // The arena lives for the complete sliding task. KSW2 frees each window's
    // blocks back into it, so the next dependent window can reuse the same
    // already-faulted pages. Destroying the arena here returns all storage as
    // soon as the task result has been produced.
    ScopedKswScratchArena scratchArena;
    std::size_t arenaQueryLength = 0;
    std::size_t arenaReferenceLength = 0;

    // Each input base appears exactly once in the final pairwise alignment;
    // reserve the maximum possible gapped length up front so long fallback
    // intervals do not repeatedly copy an ever-growing output string.
    const std::size_t maximumAlignmentLength =
            static_cast<std::size_t>(_length_of_q) +
            static_cast<std::size_t>(_length_of_d);
    _alignment_q.clear();
    _alignment_d.clear();
    _alignment_q.reserve(maximumAlignmentLength);
    _alignment_d.reserve(maximumAlignmentLength);

    //if (sqrt(_length_of_d) * sqrt(_length_of_q) <= slidingWindowSize) {
    if ((_length_of_d * 1.0 / slidingWindowSize) * (_length_of_q * 1.0 / slidingWindowSize) <= 1) {
        totalScore += minimap2AlignmentWithArena(scratchArena.get(), dna_q, dna_d, _alignment_q, _alignment_d,
                                         matchingScore, mismatchingPenalty,
                                         openGapPenalty1, extendGapPenalty1,
                                         openGapPenalty2, extendGapPenalty2);
        databaseStart = _length_of_d + 1;
        queryStart = _length_of_q + 1;
    } else {
        while (databaseStart <= _length_of_d && queryStart <= _length_of_q) {
            databaseEnd = databaseStart + slidingWindowSize;
            queryEnd = queryStart + slidingWindowSize;
            if (databaseEnd > _length_of_d) {
                databaseEnd = _length_of_d;
            }

            if (queryEnd > _length_of_q) {
                queryEnd = _length_of_q;
            }

            std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd);
            std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd);
            std::string alignment_q = "";
            std::string alignment_d = "";

            // A different matrix shape touches a different sparse subset of
            // the rectangular traceback allocation. Keeping both subsets
            // resident would make the union exceed the w^2 task envelope.
            // Reuse pages only for consecutive equal-sized full windows and
            // return them before the shorter final window.
            if (arenaQueryLength != 0 &&
                (arenaQueryLength != qSeq.size() ||
                 arenaReferenceLength != dSeq.size())) {
                scratchArena.reset();
            }
            arenaQueryLength = qSeq.size();
            arenaReferenceLength = dSeq.size();

            if (slidingWindowSize > 1073741824) {
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            } else {
                if (extendGapPenalty2 + matchingScore < 0) {
                    totalScore += alignmentMinimap2WithArena(scratchArena.get(), qSeq, dSeq, alignment_q, alignment_d, matchingScore, mismatchingPenalty + matchingScore, openGapPenalty1 + matchingScore, extendGapPenalty1 + matchingScore, openGapPenalty2 + matchingScore, extendGapPenalty2 + matchingScore, endPositionq, endPositiont);
                } else {
                    totalScore += alignmentMinimap2WithArena(scratchArena.get(), qSeq, dSeq, alignment_q, alignment_d, matchingScore, mismatchingPenalty + matchingScore - 1, openGapPenalty1 + matchingScore - 1, extendGapPenalty1 + matchingScore - 1, openGapPenalty2 + matchingScore - 1, -1, endPositionq, endPositiont);
                }
            }

            const std::size_t consumedQuery = alignment_q.size() -
                    static_cast<std::size_t>(std::count(
                            alignment_q.begin(), alignment_q.end(), '-'));
            const std::size_t consumedReference = alignment_d.size() -
                    static_cast<std::size_t>(std::count(
                            alignment_d.begin(), alignment_d.end(), '-'));
            if (static_cast<std::size_t>(endPositionq) != consumedQuery ||
                static_cast<std::size_t>(endPositiont) !=
                        consumedReference) {
                throw std::runtime_error(
                        "KSW2 semiglobal traceback did not reconstruct its "
                        "reported endpoint");
            }

            _alignment_q += alignment_q;
            _alignment_d += alignment_d;
            queryStart += endPositionq;
            databaseStart += endPositiont;
        }
    }

    int32_t final_indel_length = 0;
    int32_t count_1 = _length_of_d - databaseStart;
    if (count_1 >= 0) {
        _alignment_q += std::string(count_1 + 1, '-');
        _alignment_d += dna_d.substr(databaseStart - 1, count_1 + 1);
        final_indel_length += count_1 + 1;
    }
    assert(_alignment_d.size() == _alignment_q.size());

    int32_t count_2 = _length_of_q - queryStart;
    if (count_2 >= 0) {
        _alignment_q += dna_q.substr(queryStart - 1, count_2 + 1);
        _alignment_d += std::string(count_2 + 1, '-');
        final_indel_length += count_2 + 1;
    }
    assert(_alignment_d.size() == _alignment_q.size());

    if (final_indel_length > 0) {
        totalScore += std::max(openGapPenalty1 + extendGapPenalty1 * final_indel_length, openGapPenalty2 + extendGapPenalty2 * final_indel_length);
    }

    // The semiglobal chunks use transformed penalties to choose their
    // endpoints, so summing the per-chunk KSW scores is not the score of the
    // emitted end-to-end alignment under the Tier-1 objective. Report the
    // latter so measured fallback quality is directly comparable with full
    // DP and with the independently rescored banded-KSW output.
    return scoreGappedTwoPieceAffineAlignment(
            _alignment_q, _alignment_d,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2);
}


namespace {

struct KswDeadlineContext {
    std::chrono::steady_clock::time_point deadline;
    bool stopRequested = false;
};

int kswDeadlineProbe(void *arguments, int, int) {
    KswDeadlineContext *const context =
            static_cast<KswDeadlineContext *>(arguments);
    if (context == nullptr) {
        return 1;
    }
    if (std::chrono::steady_clock::now() < context->deadline) {
        return 1;
    }
    context->stopRequested = true;
    return 0;
}

class ScopedKswDeadline {
public:
    explicit ScopedKswDeadline(double maximumRuntimeSeconds) {
        if (maximumRuntimeSeconds > 0.0 &&
            std::isfinite(maximumRuntimeSeconds)) {
            const auto duration = std::chrono::duration<double>(
                    maximumRuntimeSeconds);
            context_.deadline = std::chrono::steady_clock::now() +
                    std::chrono::duration_cast<
                            std::chrono::steady_clock::duration>(duration);
            ksw_set_progress_callback(&kswDeadlineProbe, &context_, 256);
            active_ = true;
        }
    }

    ScopedKswDeadline(const ScopedKswDeadline &) = delete;
    ScopedKswDeadline &operator=(const ScopedKswDeadline &) = delete;

    ~ScopedKswDeadline() {
        if (active_) {
            ksw_clear_progress_callback();
        }
    }

    bool stopped() const { return context_.stopRequested; }

private:
    KswDeadlineContext context_{};
    bool active_ = false;
};

void encodeKswSequence(const std::string &sequence,
                       std::vector<uint8_t> &encoded,
                       int &flags) {
    encoded.resize(sequence.size());
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence[i]) {
            case 'A': case 'a': encoded[i] = 0; break;
            case 'C': case 'c': encoded[i] = 1; break;
            case 'G': case 'g': encoded[i] = 2; break;
            case 'T': case 't': encoded[i] = 3; break;
            default:
                encoded[i] = 4;
                flags |= KSW_EZ_GENERIC_SC;
                break;
        }
    }
}

bool gappedAlignmentReconstructs(const std::string &alignment,
                                 const std::string &sequence) {
    std::size_t position = 0;
    for (const char base : alignment) {
        if (base == '-') {
            continue;
        }
        if (position >= sequence.size() || sequence[position] != base) {
            return false;
        }
        ++position;
    }
    return position == sequence.size();
}

int64_t scoreGappedTwoPieceAffineAlignment(
        const std::string &queryAlignment,
        const std::string &referenceAlignment,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2) {
    if (queryAlignment.size() != referenceAlignment.size()) {
        throw std::runtime_error(
                "KSW2 returned alignment rows with different lengths");
    }
    int64_t score = 0;
    std::size_t column = 0;
    while (column < queryAlignment.size()) {
        const bool queryGap = queryAlignment[column] == '-';
        const bool referenceGap = referenceAlignment[column] == '-';
        if (queryGap && referenceGap) {
            throw std::runtime_error(
                    "KSW2 returned an alignment column with two gaps");
        }
        if (queryGap || referenceGap) {
            const bool gapInQuery = queryGap;
            std::size_t end = column + 1;
            while (end < queryAlignment.size() &&
                   (queryAlignment[end] == '-') == gapInQuery &&
                   (referenceAlignment[end] == '-') != gapInQuery) {
                ++end;
            }
            const int64_t length = static_cast<int64_t>(end - column);
            score += std::max<int64_t>(
                    static_cast<int64_t>(openGapPenalty1) +
                            static_cast<int64_t>(extendGapPenalty1) * length,
                    static_cast<int64_t>(openGapPenalty2) +
                            static_cast<int64_t>(extendGapPenalty2) * length);
            column = end;
            continue;
        }
        if (queryAlignment[column] != referenceAlignment[column]) {
            score += mismatchingPenalty;
        }
        ++column;
    }
    return score;
}

int64_t ksw2ScoreOnly(const std::string &dna_q,
                      const std::string &dna_d,
                      const int32_t mismatchingPenalty,
                      const int32_t openGapPenalty1,
                      const int32_t extendGapPenalty1,
                      const int32_t openGapPenalty2,
                      const int32_t extendGapPenalty2,
                      bool &stopped) {
    if (dna_q.empty() || dna_d.empty()) {
        const int64_t length = static_cast<int64_t>(
                std::max(dna_q.size(), dna_d.size()));
        return std::max<int64_t>(
                static_cast<int64_t>(openGapPenalty1) +
                        static_cast<int64_t>(extendGapPenalty1) * length,
                static_cast<int64_t>(openGapPenalty2) +
                        static_cast<int64_t>(extendGapPenalty2) * length);
    }
    if (dna_q.size() > static_cast<std::size_t>(INT_MAX) ||
        dna_d.size() > static_cast<std::size_t>(INT_MAX)) {
        stopped = false;
        return KSW_NEG_INF;
    }
    const int8_t match = 0;
    const int8_t mismatch = static_cast<int8_t>(
            mismatchingPenalty < 0
            ? mismatchingPenalty : -mismatchingPenalty);
    int8_t matrix[25] = {
            match, mismatch, mismatch, mismatch, mismatch,
            mismatch, match, mismatch, mismatch, mismatch,
            mismatch, mismatch, match, mismatch, mismatch,
            mismatch, mismatch, mismatch, match, mismatch,
            mismatch, mismatch, mismatch, mismatch, match};
    int flags = KSW_EZ_SCORE_ONLY;
    std::vector<uint8_t> query;
    std::vector<uint8_t> target;
    encodeKswSequence(dna_q, query, flags);
    encodeKswSequence(dna_d, target, flags);
    ksw_extz_t result;
    std::memset(&result, 0, sizeof(result));
#ifdef __AVX512BW__
    ksw_extd2_avx512(
            nullptr, static_cast<int>(query.size()), query.data(),
            static_cast<int>(target.size()), target.data(),
            5, matrix,
            static_cast<int8_t>(-openGapPenalty1),
            static_cast<int8_t>(-extendGapPenalty1),
            static_cast<int8_t>(-openGapPenalty2),
            static_cast<int8_t>(-extendGapPenalty2),
            -1, -1, 0, flags, &result);
#elif __AVX2__
    ksw_extd2_avx2(
            nullptr, static_cast<int>(query.size()), query.data(),
            static_cast<int>(target.size()), target.data(),
            5, matrix,
            static_cast<int8_t>(-openGapPenalty1),
            static_cast<int8_t>(-extendGapPenalty1),
            static_cast<int8_t>(-openGapPenalty2),
            static_cast<int8_t>(-extendGapPenalty2),
            -1, -1, 0, flags, &result);
#else
    ksw_extd2_sse(
            nullptr, static_cast<int>(query.size()), query.data(),
            static_cast<int>(target.size()), target.data(),
            5, matrix,
            static_cast<int8_t>(-openGapPenalty1),
            static_cast<int8_t>(-extendGapPenalty1),
            static_cast<int8_t>(-openGapPenalty2),
            static_cast<int8_t>(-extendGapPenalty2),
            -1, -1, 0, flags, &result);
#endif
    stopped = result.stopped != 0;
    std::free(result.cigar);
    return result.score;
}

}  // namespace

int64_t alignSlidingWindow_minimap2(const std::string &dna_q, const std::string &dna_d, int64_t _length_of_q, int64_t _length_of_d,
                                    std::string &_alignment_q, std::string &_alignment_d, const int64_t &slidingWindowSize,
                                    const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                                    const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    int8_t a = 0, b = mismatchingPenalty < 0 ? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, b,
                      b, a, b, b, b,
                      b, b, a, b, b,
                      b, b, b, a, b,
                      b, b, b, b, a};
    int tl = dna_d.length(), ql = dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = dna_d.c_str();
    const char *qseq = dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);

    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table

    ts = (uint8_t *) malloc(tl);
    qs = (uint8_t *) malloc(ql);

    int i;
    int kswFlags = 0;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
        if (ts[i] == 4) kswFlags |= KSW_EZ_GENERIC_SC;
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t) qseq[i]];
        if (qs[i] == 4) kswFlags |= KSW_EZ_GENERIC_SC;
    }


#ifdef __AVX512BW__
    ksw_extd2_avx512(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, kswFlags, & ez);
#elif __AVX2__
    ksw_extd2_avx2(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, kswFlags, &ez);
#else
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, kswFlags, &ez);
#endif

    _alignment_q.clear();
    _alignment_d.clear();
    _alignment_q.reserve(dna_q.size() + dna_d.size());
    _alignment_d.reserve(_alignment_q.capacity());
    int pattern_pos = 0;
    int text_pos = 0;
    for (i = 0; i < ez.n_cigar; ++i) {
        const int cLen = static_cast<int>(ez.cigar[i] >> 4);
        const int operation = static_cast<int>(ez.cigar[i] & 0xf);
        if (operation == 0 || operation == 7 || operation == 8) {
            _alignment_q.append(dna_q, text_pos, cLen);
            _alignment_d.append(dna_d, pattern_pos, cLen);
            pattern_pos += cLen;
            text_pos += cLen;
        } else if (operation == 1) {
            _alignment_q.append(dna_q, text_pos, cLen);
            _alignment_d.append(static_cast<std::size_t>(cLen), '-');
            text_pos += cLen;
        } else if (operation == 2) {
            _alignment_q.append(static_cast<std::size_t>(cLen), '-');
            _alignment_d.append(dna_d, pattern_pos, cLen);
            pattern_pos += cLen;
        }
    }
    assert(text_pos == ql);
    assert(pattern_pos == tl);

    // With a restrictive band, KSW2 may report a frontier score that does not
    // equal the score of the emitted end-to-end CIGAR.  Selection and
    // score-certified comparison must use the actual alignment objective, not
    // that internal frontier value.  Re-scoring is linear and negligible next
    // to the DP work.
    const int64_t totalScore = scoreGappedTwoPieceAffineAlignment(
            _alignment_q, _alignment_d,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2);
    free(ez.cigar);
    free(ts);
    free(qs);

    return totalScore;
}

static int64_t alignKsw2SingletrackImpl(
        const std::string &dna_q, const std::string &dna_d,
        std::string &alignment_q, std::string &alignment_d,
        int64_t bandWidth,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2,
        bool *stopped) {
    if (stopped != nullptr) *stopped = false;
    if (dna_q.empty() || dna_d.empty()) {
        alignment_q = dna_q;
        alignment_d = dna_d;
        const std::size_t length = std::max(dna_q.size(), dna_d.size());
        if (dna_q.empty()) alignment_q.assign(length, '-');
        if (dna_d.empty()) alignment_d.assign(length, '-');
        return std::max<int64_t>(
                static_cast<int64_t>(openGapPenalty1) +
                        static_cast<int64_t>(extendGapPenalty1) * length,
                static_cast<int64_t>(openGapPenalty2) +
                        static_cast<int64_t>(extendGapPenalty2) * length);
    }
    if (dna_q.size() > static_cast<std::size_t>(INT_MAX) ||
        dna_d.size() > static_cast<std::size_t>(INT_MAX) ||
        dna_q.size() + dna_d.size() - 1 >
                static_cast<std::size_t>(INT_MAX)) {
        throw std::length_error("KSW2-Singletrack input exceeds INT_MAX");
    }
    if (bandWidth < -1 || bandWidth > INT_MAX) {
        throw std::invalid_argument(
                "KSW2-Singletrack band width must be -1..INT_MAX");
    }
    const int8_t match = 0;
    const int8_t mismatch = static_cast<int8_t>(
            mismatchingPenalty < 0
            ? mismatchingPenalty : -mismatchingPenalty);
    int8_t matrix[25] = {
            match, mismatch, mismatch, mismatch, mismatch,
            mismatch, match, mismatch, mismatch, mismatch,
            mismatch, mismatch, match, mismatch, mismatch,
            mismatch, mismatch, mismatch, match, mismatch,
            mismatch, mismatch, mismatch, mismatch, match};
    int flags = 0;
    std::vector<uint8_t> query;
    std::vector<uint8_t> target;
    encodeKswSequence(dna_q, query, flags);
    encodeKswSequence(dna_d, target, flags);
    ksw_extz_t result;
    std::memset(&result, 0, sizeof(result));
#ifdef __AVX512BW__
    ksw_extd2_singletrack_avx512(
#elif defined(__AVX2__)
    ksw_extd2_singletrack_avx2(
#else
    ksw_extd2_singletrack_sse(
#endif
            nullptr, static_cast<int>(query.size()), query.data(),
            static_cast<int>(target.size()), target.data(),
            5, matrix,
            static_cast<int8_t>(-openGapPenalty1),
            static_cast<int8_t>(-extendGapPenalty1),
            static_cast<int8_t>(-openGapPenalty2),
            static_cast<int8_t>(-extendGapPenalty2),
            static_cast<int>(bandWidth), -1, 0, flags, &result);
    if (result.stopped || result.zdropped) {
        if (stopped != nullptr) *stopped = result.stopped != 0;
        std::free(result.cigar);
        alignment_q.clear();
        alignment_d.clear();
        return KSW_NEG_INF;
    }
    alignment_q.clear();
    alignment_d.clear();
    alignment_q.reserve(dna_q.size() + dna_d.size());
    alignment_d.reserve(alignment_q.capacity());
    std::size_t queryPosition = 0;
    std::size_t targetPosition = 0;
    for (int i = 0; i < result.n_cigar; ++i) {
        const std::size_t length = result.cigar[i] >> 4;
        const int operation = result.cigar[i] & 0xf;
        if (operation == 0 || operation == 7 || operation == 8) {
            alignment_q.append(dna_q, queryPosition, length);
            alignment_d.append(dna_d, targetPosition, length);
            queryPosition += length;
            targetPosition += length;
        } else if (operation == 1) {
            alignment_q.append(dna_q, queryPosition, length);
            alignment_d.append(length, '-');
            queryPosition += length;
        } else if (operation == 2) {
            alignment_q.append(length, '-');
            alignment_d.append(dna_d, targetPosition, length);
            targetPosition += length;
        } else {
            std::free(result.cigar);
            alignment_q.clear();
            alignment_d.clear();
            return KSW_NEG_INF;
        }
    }
    std::free(result.cigar);
    // Bounds/endpoint checks are part of safe CIGAR decoding, not an
    // engine-comparison validator. They do not trigger a special KSW2-full
    // rescue; any structural/runtime failure follows the ordinary exact-tier
    // selection loop.
    if (queryPosition != dna_q.size() || targetPosition != dna_d.size()) {
        alignment_q.clear();
        alignment_d.clear();
        return KSW_NEG_INF;
    }
    // The Singletrack traceback stores only score differences. Its internal
    // KSW2 endpoint score is not a reliable public score after reconstructing
    // from those tracks (this is observable even when the emitted CIGAR is
    // byte-identical to KSW2-full). As in the established KSW2 wrapper, always
    // report the objective value of the emitted end-to-end alignment. This is
    // linear in the CIGAR length and negligible beside the DP pass.
    return scoreGappedTwoPieceAffineAlignment(
            alignment_q, alignment_d, mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2);
}

int64_t alignKsw2Singletrack(
        const std::string &dna_q, const std::string &dna_d,
        std::string &alignment_q, std::string &alignment_d,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2,
        bool *stopped) {
    return alignKsw2SingletrackImpl(
            dna_q, dna_d, alignment_q, alignment_d, -1,
            mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            stopped);
}

int64_t alignKsw2SingletrackBanded(
        const std::string &dna_q, const std::string &dna_d,
        std::string &alignment_q, std::string &alignment_d,
        int64_t bandWidth,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2,
        bool *stopped) {
    return alignKsw2SingletrackImpl(
            dna_q, dna_d, alignment_q, alignment_d, bandWidth,
            mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2,
            stopped);
}

namespace {

enum class ChainedTracebackStatus {
    Unavailable,
    Completed,
    TimeLimit
};

ChainedTracebackStatus tryChainedKsw2Traceback(
        const std::string &dna_q,
        const std::string &dna_d,
        const anchorwave::AlignmentSelectionPlan &selectionPlan,
        uint64_t maximumBand,
        const ScopedKswDeadline &deadline,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        std::string &queryAlignment,
        std::string &referenceAlignment,
        int64_t &score) {
    if (selectionPlan.certificationAnchors.size() < 2 ||
        maximumBand == 0) {
        return ChainedTracebackStatus::Unavailable;
    }
    queryAlignment.clear();
    referenceAlignment.clear();
    queryAlignment.reserve(dna_q.size() + dna_d.size());
    referenceAlignment.reserve(dna_q.size() + dna_d.size());

    const auto appendSegment = [&](std::size_t queryStart,
                                   std::size_t queryEnd,
                                   std::size_t referenceStart,
                                   std::size_t referenceEnd) {
        const std::size_t queryLength = queryEnd - queryStart;
        const std::size_t referenceLength = referenceEnd - referenceStart;
        if (queryLength == 0 || referenceLength == 0) {
            if (queryLength != 0) {
                queryAlignment.append(dna_q, queryStart, queryLength);
                referenceAlignment.append(queryLength, '-');
            } else if (referenceLength != 0) {
                queryAlignment.append(referenceLength, '-');
                referenceAlignment.append(
                        dna_d, referenceStart, referenceLength);
            }
            return true;
        }
        const uint64_t difference = queryLength > referenceLength
                ? queryLength - referenceLength
                : referenceLength - queryLength;
        if (difference > maximumBand ||
            queryLength > static_cast<std::size_t>(INT_MAX) ||
            referenceLength > static_cast<std::size_t>(INT_MAX)) {
            return false;
        }
        const uint64_t localBand = std::min<uint64_t>(
                std::max(queryLength, referenceLength), maximumBand);
        std::string localQuery;
        std::string localReference;
        const std::string querySegment = dna_q.substr(queryStart, queryLength);
        const std::string referenceSegment = dna_d.substr(
                referenceStart, referenceLength);
        alignSlidingWindow_minimap2(
                querySegment, referenceSegment,
                static_cast<int64_t>(queryLength),
                static_cast<int64_t>(referenceLength),
                localQuery, localReference,
                static_cast<int64_t>(localBand),
                mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2);
        if (deadline.stopped() ||
            localQuery.size() != localReference.size() ||
            !gappedAlignmentReconstructs(localQuery, querySegment) ||
            !gappedAlignmentReconstructs(localReference, referenceSegment)) {
            return false;
        }
        queryAlignment += localQuery;
        referenceAlignment += localReference;
        return true;
    };

    std::size_t queryPosition = 0;
    std::size_t referencePosition = 0;
    std::size_t anchorsUsed = 0;
    for (const anchorwave::AlignmentChainAnchor &anchor :
         selectionPlan.certificationAnchors) {
        const std::size_t anchorQuery = anchor.queryPosition;
        const std::size_t anchorReference = anchor.referencePosition;
        const std::size_t anchorLength = anchor.length;
        if (anchorLength == 0 || anchorQuery < queryPosition ||
            anchorReference < referencePosition ||
            anchorQuery + anchorLength > dna_q.size() ||
            anchorReference + anchorLength > dna_d.size()) {
            continue;
        }
        bool matchingAnchor = true;
        for (std::size_t offset = 0; offset < anchorLength; ++offset) {
            if (std::toupper(static_cast<unsigned char>(
                        dna_q[anchorQuery + offset])) !=
                std::toupper(static_cast<unsigned char>(
                        dna_d[anchorReference + offset]))) {
                matchingAnchor = false;
                break;
            }
        }
        if (!matchingAnchor) {
            continue;
        }
        if (!appendSegment(queryPosition, anchorQuery,
                           referencePosition, anchorReference)) {
            queryAlignment.clear();
            referenceAlignment.clear();
            return deadline.stopped()
                   ? ChainedTracebackStatus::TimeLimit
                   : ChainedTracebackStatus::Unavailable;
        }
        queryAlignment.append(dna_q, anchorQuery, anchorLength);
        referenceAlignment.append(dna_d, anchorReference, anchorLength);
        queryPosition = anchorQuery + anchorLength;
        referencePosition = anchorReference + anchorLength;
        ++anchorsUsed;
    }
    if (anchorsUsed < 2 ||
        !appendSegment(queryPosition, dna_q.size(),
                       referencePosition, dna_d.size())) {
        queryAlignment.clear();
        referenceAlignment.clear();
        return deadline.stopped()
               ? ChainedTracebackStatus::TimeLimit
               : ChainedTracebackStatus::Unavailable;
    }
    if (!gappedAlignmentReconstructs(queryAlignment, dna_q) ||
        !gappedAlignmentReconstructs(referenceAlignment, dna_d)) {
        queryAlignment.clear();
        referenceAlignment.clear();
        return ChainedTracebackStatus::Unavailable;
    }
    score = scoreGappedTwoPieceAffineAlignment(
            queryAlignment, referenceAlignment,
            mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2);
    return ChainedTracebackStatus::Completed;
}

}  // namespace

anchorwave::Ksw2CertifiedResult alignScoreCertifiedKsw2(
        const std::string &dna_q, const std::string &dna_d,
        std::string &alignment_q, std::string &alignment_d,
        int64_t initialBandWidth, int64_t maximumBandWidth,
        double maximumRuntimeSeconds,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2,
        const anchorwave::AlignmentSelectionPlan *selectionPlan,
        const int64_t *precomputedOptimalScore) {
    anchorwave::Ksw2CertifiedResult result;
    alignment_q.clear();
    alignment_d.clear();
    const uint64_t minimumLength = std::min(dna_q.size(), dna_d.size());
    const uint64_t lengthDifference = dna_q.size() > dna_d.size()
            ? dna_q.size() - dna_d.size()
            : dna_d.size() - dna_q.size();
    if (minimumLength == 0) {
        alignment_q = dna_q.empty()
                ? std::string(dna_d.size(), '-') : dna_q;
        alignment_d = dna_d.empty()
                ? std::string(dna_q.size(), '-') : dna_d;
        bool stopped = false;
        result.score = ksw2ScoreOnly(
                dna_q, dna_d, mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2, stopped);
        result.optimalScore = result.score;
        result.finalBandWidth = 0;
        result.tracebackAttempts = 1;
        result.status = anchorwave::Ksw2CertifiedStatus::Completed;
        return result;
    }
    if (initialBandWidth < 0 || maximumBandWidth < 0 ||
        static_cast<uint64_t>(maximumBandWidth) < lengthDifference) {
        result.status = anchorwave::Ksw2CertifiedStatus::NotCertified;
        return result;
    }

    ScopedKswDeadline deadline(maximumRuntimeSeconds);
    bool scoreStopped = false;
    if (precomputedOptimalScore != nullptr) {
        result.optimalScore = *precomputedOptimalScore;
    } else {
        result.optimalScore = ksw2ScoreOnly(
                dna_q, dna_d, mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2, scoreStopped);
    }
    if (scoreStopped || deadline.stopped()) {
        result.status = anchorwave::Ksw2CertifiedStatus::TimeLimit;
        return result;
    }
    if (result.optimalScore == KSW_NEG_INF) {
        result.status = anchorwave::Ksw2CertifiedStatus::Failed;
        return result;
    }

    const uint64_t maximumBand = std::min<uint64_t>(
            std::max(dna_q.size(), dna_d.size()),
            static_cast<uint64_t>(maximumBandWidth));
    if (selectionPlan != nullptr) {
        std::string chainedQuery;
        std::string chainedReference;
        int64_t chainedScore = 0;
        ++result.tracebackAttempts;
        const ChainedTracebackStatus chainedStatus = tryChainedKsw2Traceback(
                dna_q, dna_d, *selectionPlan, maximumBand, deadline,
                mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2,
                chainedQuery, chainedReference, chainedScore);
        if (chainedStatus == ChainedTracebackStatus::TimeLimit) {
            result.status = anchorwave::Ksw2CertifiedStatus::TimeLimit;
            return result;
        }
        if (chainedStatus == ChainedTracebackStatus::Completed) {
            result.score = chainedScore;
            result.bestEffortQueryAlignment = chainedQuery;
            result.bestEffortReferenceAlignment = chainedReference;
            if (chainedScore == result.optimalScore) {
                result.finalBandWidth = -2;
                alignment_q = std::move(chainedQuery);
                alignment_d = std::move(chainedReference);
                result.bestEffortQueryAlignment.clear();
                result.bestEffortReferenceAlignment.clear();
                result.status = anchorwave::Ksw2CertifiedStatus::Completed;
                return result;
            }
        }
    }
    // In the 100-Mb B73/Mo17 trace, none of the 31 rescue intervals uniquely
    // certified at the smallest planned band. Five monolithic successes were
    // reached at 2x, 4x or the cap, while a direct jump to the cap lost two
    // certifications because this KSW2 traceback implementation is not
    // path-monotone under all equal-score ties. Skip only the empirically
    // unproductive first band, then retain geometric widening. This removes
    // one complete DP pass without weakening the observed exact coverage.
    const uint64_t requestedInitial = std::max<uint64_t>(
            lengthDifference, std::max<int64_t>(1, initialBandWidth));
    const uint64_t doubledInitial = requestedInitial >
            std::numeric_limits<uint64_t>::max() / 2
            ? maximumBand : requestedInitial * 2;
    uint64_t band = std::min(maximumBand,
                             std::max(requestedInitial, doubledInitial));
    while (band >= lengthDifference && band <= maximumBand) {
        std::string candidateQuery;
        std::string candidateReference;
        ++result.tracebackAttempts;
        const int64_t bandedScore = alignSlidingWindow_minimap2(
                dna_q, dna_d,
                static_cast<int64_t>(dna_q.size()),
                static_cast<int64_t>(dna_d.size()),
                candidateQuery, candidateReference,
                static_cast<int64_t>(band),
                mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2);
        if (deadline.stopped()) {
            result.status = anchorwave::Ksw2CertifiedStatus::TimeLimit;
            return result;
        }
        const bool validReconstruction =
                candidateQuery.size() == candidateReference.size() &&
                gappedAlignmentReconstructs(candidateQuery, dna_q) &&
                gappedAlignmentReconstructs(candidateReference, dna_d);
        if (bandedScore == result.optimalScore && validReconstruction) {
            result.score = bandedScore;
            result.finalBandWidth = static_cast<int64_t>(band);
            alignment_q = std::move(candidateQuery);
            alignment_d = std::move(candidateReference);
            result.status = anchorwave::Ksw2CertifiedStatus::Completed;
            return result;
        }
        if (validReconstruction &&
            (result.bestEffortQueryAlignment.empty() ||
             bandedScore > result.score)) {
            result.score = bandedScore;
            result.bestEffortQueryAlignment = std::move(candidateQuery);
            result.bestEffortReferenceAlignment = std::move(
                    candidateReference);
        }
        if (band == maximumBand) {
            break;
        }
        const uint64_t doubled = band >
                std::numeric_limits<uint64_t>::max() / 2
                ? maximumBand : band * 2;
        const uint64_t nextBand = std::min(
                maximumBand, std::max(band + 1, doubled));
        if (nextBand <= band) {
            break;
        }
        band = nextBand;
    }
    // A failed certification still supplies the approximate tier with the
    // highest-scoring already-computed banded alignment. Evaluate the skipped
    // smallest band only on this failure path so the retained fallback sees
    // the same candidate-band set as before this optimization. Successful
    // exact rescues avoid this DP entirely.
    if (requestedInitial < std::min(maximumBand, doubledInitial)) {
        std::string candidateQuery;
        std::string candidateReference;
        ++result.tracebackAttempts;
        const int64_t bandedScore = alignSlidingWindow_minimap2(
                dna_q, dna_d,
                static_cast<int64_t>(dna_q.size()),
                static_cast<int64_t>(dna_d.size()),
                candidateQuery, candidateReference,
                static_cast<int64_t>(requestedInitial),
                mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2);
        if (deadline.stopped()) {
            result.status = anchorwave::Ksw2CertifiedStatus::TimeLimit;
            return result;
        }
        const bool validReconstruction =
                candidateQuery.size() == candidateReference.size() &&
                gappedAlignmentReconstructs(candidateQuery, dna_q) &&
                gappedAlignmentReconstructs(candidateReference, dna_d);
        if (bandedScore == result.optimalScore && validReconstruction) {
            result.score = bandedScore;
            result.finalBandWidth = static_cast<int64_t>(requestedInitial);
            alignment_q = std::move(candidateQuery);
            alignment_d = std::move(candidateReference);
            result.status = anchorwave::Ksw2CertifiedStatus::Completed;
            return result;
        }
        if (validReconstruction &&
            (result.bestEffortQueryAlignment.empty() ||
             bandedScore > result.score)) {
            result.score = bandedScore;
            result.bestEffortQueryAlignment = std::move(candidateQuery);
            result.bestEffortReferenceAlignment = std::move(
                    candidateReference);
        }
    }
    result.finalBandWidth = static_cast<int64_t>(band);
    result.status = anchorwave::Ksw2CertifiedStatus::NotCertified;
    return result;
}

int64_t alignSlidingWindow_minimap2_or_NW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d, std::string &_alignMethod,
                                          const int64_t &slidingWindowSize, const int32_t &matchingScore,
                                          const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {

    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns end
    int32_t longerSeqLength = std::max(_length_of_d, _length_of_q);

    if ((_length_of_d * 1.0 / slidingWindowSize) * (_length_of_q * 1.0 / slidingWindowSize) <= 1) { //  _length_of_d*_length_of_q <= (slidingWindowSize*slidingWindowSize) this calculated via RAM cost
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = -1;
        totalScore = alignSlidingWindow_minimap2(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                 mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        _alignMethod = "MINIMAP2";
    } else if (longerSeqLength * 0.5 < slidingWindowSize) { // this calculated with RAM cost longerSeqLength*slidingWindowSize*0.5 <= (slidingWindowSize*slidingWindowSize
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = (slidingWindowSize * 1.0 / longerSeqLength) * slidingWindowSize;
        //std::cout << "adjustBandWidth:" << adjustBandWidth << std::endl;
        totalScore = alignSlidingWindow_minimap2(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                 mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        _alignMethod = "BANDED_MINIMAP2";
    } else {
        totalScore = alignSlidingWindow(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                        mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        _alignMethod = "SLIDING_WINDOW";
    }

    return totalScore;
}

std::shared_ptr<const anchorwave::AlignmentSelectionPlan>
prepareAlignmentSelectionPlan(
        const std::string &dna_q, const std::string &dna_d,
        const int64_t &slidingWindowSize,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2) {
    const bool memorySchedulingEnabled =
            anchorwave::alignmentMemorySchedulingEnabled();
    return std::make_shared<const anchorwave::AlignmentSelectionPlan>(
            anchorwave::makeAlignmentSelectionPlan(
                    dna_q, dna_d, slidingWindowSize,
                    mismatchingPenalty,
                    openGapPenalty1, extendGapPenalty1,
                    openGapPenalty2, extendGapPenalty2,
                    anchorwave::configuredExactAlignmentMaximumEstimatedMinutes(),
                    0, 0,
                    memorySchedulingEnabled));
}

int64_t alignSlidingWindow_local_wfa2_v2(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d, std::string &_alignMethod,
                                         const int64_t &slidingWindowSize, const int32_t &matchingScore,
                                         const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                         const std::shared_ptr<const anchorwave::AlignmentSelectionPlan> &precomputedSelection) {
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";

    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns begin
    bool flag_q_all_N = std::all_of(dna_q.begin(), dna_q.end(), [](char c) { return c == 'N'; });
    bool flag_d_all_N = std::all_of(dna_d.begin(), dna_d.end(), [](char c) { return c == 'N'; });

    if (flag_q_all_N || flag_d_all_N) {
        _alignment_q = dna_q;
        _alignment_d = dna_d;

        int32_t count_ = abs(_length_of_q - _length_of_d);
        if (_length_of_q < _length_of_d) {
            _alignment_q += std::string(count_, '-');
        }

        if (_length_of_d < _length_of_q) {
            _alignment_d += std::string(count_, '-');
        }
        _alignMethod = "FILLING";
        return totalScore;
    }

    const uint64_t wfaMemoryBudget =
            anchorwave::wfaMemoryBudgetBytes(slidingWindowSize);
    const bool memorySchedulingEnabled =
            anchorwave::alignmentMemorySchedulingEnabled();
    const uint64_t transientMemoryBytes =
            anchorwave::alignmentTaskTransientMemoryBytes(
                    static_cast<uint64_t>(dna_q.size()),
                    static_cast<uint64_t>(dna_d.size()));
    // -w^2 is the normal single-task algorithm budget.  The tier-wide
    // qlen*rlen <= w^2 quality exception is expanded later from each exact
    // engine's own prediction; process-wide -M admission remains the hard
    // boundary.
    const uint64_t elasticHighWfaBudget = 0;
    const uint64_t elasticFullKsw2Budget = 0;
    std::shared_ptr<PreparedAlignmentSelection> prepared;
    if (anchorwave::currentAnchorTaskCanDefer()) {
        const std::shared_ptr<void> retained =
                anchorwave::currentAnchorTaskRetainedState();
        if (retained) {
            prepared = std::static_pointer_cast<PreparedAlignmentSelection>(
                    retained);
            // Retained state belongs to this one deferrable closure. Production
            // gap retries reuse the same cached string objects, so the plan was
            // already content/provenance-validated on the first attempt. Avoid
            // rehashing both long sequences on every memory retry. A different
            // object or changed length takes the full validation path below.
            if (prepared->queryObject != &dna_q ||
                prepared->referenceObject != &dna_d ||
                prepared->queryLength != dna_q.size() ||
                prepared->referenceLength != dna_d.size()) {
                prepared.reset();
            }
        }
    }
    if (!prepared) {
        prepared = std::make_shared<PreparedAlignmentSelection>();
        prepared->queryObject = &dna_q;
        prepared->referenceObject = &dna_d;
        prepared->queryLength = dna_q.size();
        prepared->referenceLength = dna_d.size();
        prepared->intervalId = anchorwave::nextAlignmentTraceIntervalId();
        if (precomputedSelection &&
            anchorwave::alignmentSelectionPlanMatches(
                    *precomputedSelection, dna_q, dna_d,
                    slidingWindowSize, mismatchingPenalty,
                    openGapPenalty1, extendGapPenalty1,
                    openGapPenalty2, extendGapPenalty2,
                    anchorwave::configuredExactAlignmentMaximumEstimatedMinutes(),
                    elasticHighWfaBudget, elasticFullKsw2Budget,
                    memorySchedulingEnabled)) {
            prepared->plan = *precomputedSelection;
        } else {
            prepared->plan = *prepareAlignmentSelectionPlan(
                    dna_q, dna_d, slidingWindowSize,
                    mismatchingPenalty,
                    openGapPenalty1, extendGapPenalty1,
                    openGapPenalty2, extendGapPenalty2);
        }
        if (anchorwave::currentAnchorTaskCanDefer()) {
            anchorwave::setCurrentAnchorTaskRetainedState(prepared);
        }
    }
    const anchorwave::AlignmentSelectionPlan &selection = prepared->plan;
    const bool forceExactTier = selection.exactCellEnvelopeException;
    const bool balancedTimeLimitEnabled =
            selection.provenance.exactAlignmentMaximumEstimatedMinutes > 0.0 &&
            !forceExactTier;
    // Tier 1: exact dynamic programming. Every admitted exact candidate is
    // exhausted before the executor may enter a lower-quality tier. With
    // -bt 0 the plan is exact-first; a positive -bt uses the balanced
    // P50/P90 time policy while the plan is built. Temporary memory pressure
    // never changes the quality tier.
    const auto recordSelectionOnce = [&]() {
        if (!prepared->selectionRecorded) {
            anchorwave::recordAlignmentSelectionPlan(selection);
            prepared->selectionRecorded = true;
        }
    };
    const auto recordAttempt = [&, prepared](
            anchorwave::AlignmentCandidate candidate,
            const anchorwave::AlgorithmCostEstimate *estimate,
            const CandidateMemoryEnvelope &memoryEnvelope,
            double predictedWaitMinutes,
            double actualSeconds,
            uint64_t actualMemoryBytes,
            const std::string &status) {
        if (!anchorwave::alignmentTraceEnabled()) {
            return;
        }
        anchorwave::AlignmentAttemptTrace trace;
        trace.intervalId = prepared->intervalId;
        trace.attempt = ++prepared->attemptCount;
        trace.profile = selection.profile;
        trace.provenance = selection.provenance;
        trace.approximateDecision = selection.approximateDecision;
        trace.approximateKsw2BandWidth = selection.ksw2BandWidth;
        trace.approximateSlidingWindowWidth = selection.slidingWindowWidth;
        trace.candidate = candidate;
        if (estimate != nullptr) {
            trace.predictedMinutesP50 = estimate->estimatedMinutes;
            trace.predictedMinutesP90 = estimate->estimatedMinutesP90;
            trace.predictedMemoryBytes = estimate->memoryBytes;
        }
        trace.reservedMemoryBytes = memoryEnvelope.reservationBytes;
        trace.predictedWaitMinutes = predictedWaitMinutes;
        trace.actualSeconds = actualSeconds;
        trace.configuredExactLimitMinutes = balancedTimeLimitEnabled
                ? selection.provenance.exactAlignmentMaximumEstimatedMinutes
                : 0.0;
        trace.exactRuntimeSpentSeconds =
                prepared->exactRuntimeSpentSeconds;
        trace.exactRuntimeRemainingSeconds =
                trace.configuredExactLimitMinutes > 0.0
                ? std::max(0.0,
                           trace.configuredExactLimitMinutes * 60.0 -
                                   trace.exactRuntimeSpentSeconds)
                : 0.0;
        trace.actualScore = status == "completed" ? totalScore : 0;
        if (candidate == anchorwave::AlignmentCandidate::Ksw2ScoreCertified) {
            trace.certifiedOptimalScore = prepared->certifiedOptimalScore;
            trace.certifiedInitialBandWidth =
                    selection.ksw2CertifiedInitialBandWidth;
            trace.certifiedMaximumBandWidth =
                    selection.ksw2CertifiedMaximumBandWidth;
            trace.certifiedFinalBandWidth =
                    prepared->certifiedFinalBandWidth;
            trace.certifiedTracebackAttempts =
                    prepared->certifiedTracebackAttempts;
        }
        trace.actualMemoryBytes = actualMemoryBytes;
        trace.processResidentBytes =
                anchorwave::currentProcessResidentBytes();
        const anchorwave::AlignmentResourceSnapshot resources =
                anchorwave::alignmentResourceSnapshot();
        trace.processMemoryLimitBytes = resources.maxProcessMemoryBytes;
        trace.activeReservedBytes = resources.activeReservedBytes;
        const anchorwave::AnchorTaskExecutor::LoadSnapshot load =
                anchorwave::currentAnchorTaskLoadSnapshot();
        trace.workerCount = load.workerCount;
        trace.activeTasks = load.activeTasks;
        trace.readyTasks = load.readyTasks;
        trace.deferredTasks = load.deferredTasks;
        trace.plannedTasks = load.plannedTasks;
        trace.schedulableTasks = load.schedulableTasks;
        trace.globalFutureTasks = load.globalFutureTasks;
        trace.blockedOrderedHeads = load.blockedOrderedHeads;
        trace.outstandingTasks = load.outstandingTasks;
        trace.outstandingEstimatedCost = load.outstandingEstimatedCost;
        trace.criticalEstimatedCost = load.criticalEstimatedCost;
        trace.globalTailPhase = load.globalTailPhase;
        trace.admissionTailPhase = load.admissionTailPhase;
        trace.tailPhase = load.tailPhase;
        if (status == "completed") {
            trace.resultMethod =
                    anchorwave::alignmentMethodBedLabel(_alignMethod);
        }
        trace.status = status;
        anchorwave::recordAlignmentAttemptTrace(trace);
    };
    if (!prepared->planTraced && anchorwave::alignmentTraceEnabled()) {
        for (const anchorwave::AlgorithmCostEstimate &estimate :
                selection.estimates) {
            const CandidateMemoryEnvelope memoryEnvelope =
                    candidateMemoryEnvelope(
                            selection, estimate.candidate,
                            wfaMemoryBudget, transientMemoryBytes);
            const double wait = estimate.memoryFeasible
                    ? anchorwave::estimatedAlignmentMemoryWaitMinutes(
                              memoryEnvelope.reservationBytes)
                    : std::numeric_limits<double>::infinity();
            const std::string status = !estimate.memoryFeasible
                    ? "planned_memory_infeasible"
                    : !estimate.timeFeasible
                      ? "planned_time_infeasible"
                      : "planned_feasible";
            recordAttempt(estimate.candidate, &estimate, memoryEnvelope,
                          wait, 0.0, 0, status);
        }
        prepared->planTraced = true;
    }
    const auto remainingExactRuntimeSeconds = [&, prepared]() {
        const double limitMinutes =
                selection.provenance.exactAlignmentMaximumEstimatedMinutes;
        if (!balancedTimeLimitEnabled) {
            return 0.0;
        }
        return std::max(0.0, limitMinutes * 60.0 -
                              prepared->exactRuntimeSpentSeconds);
    };
    const auto attemptWfaCandidate = [&, memorySchedulingEnabled](
            anchorwave::AlignmentCandidate candidate,
            const CandidateMemoryEnvelope &memoryEnvelope,
            double maximumRuntimeSeconds,
            uint64_t &actualMemoryBytes,
            std::string &statusText) {
        anchorwave::WfaAlignmentResult wfaResult;
        anchorwave::WfaAlignmentStatus status =
                anchorwave::WfaAlignmentStatus::Failed;
        anchorwave::WfaExecutionOptions executionOptions;
        if (maximumRuntimeSeconds > 0.0) {
            executionOptions.maximumRuntimeMilliseconds =
                    static_cast<uint64_t>(std::max(
                            1.0, std::ceil(
                                    maximumRuntimeSeconds * 1000.0)));
        }
        switch (candidate) {
            case anchorwave::AlignmentCandidate::SingletrackWfa:
                status = anchorwave::alignWithSingletrackWfa(
                        dna_q, dna_d,
                        memorySchedulingEnabled
                        ? memoryEnvelope.algorithmBudgetBytes
                        : wfaMemoryBudget,
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2, wfaResult,
                        executionOptions);
                break;
            case anchorwave::AlignmentCandidate::StandardWfa:
                status = anchorwave::alignWithStandardWfa(
                        dna_q, dna_d,
                        memorySchedulingEnabled
                        ? memoryEnvelope.algorithmBudgetBytes
                        : anchorwave::standardWfaTrialMemoryBudgetBytes(
                                  wfaMemoryBudget),
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2, wfaResult,
                        executionOptions);
                break;
            case anchorwave::AlignmentCandidate::MediumWfa:
                status = anchorwave::alignWithMediumWfa(
                        dna_q, dna_d,
                        memorySchedulingEnabled
                        ? memoryEnvelope.algorithmBudgetBytes
                        : wfaMemoryBudget,
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2, wfaResult,
                        executionOptions);
                break;
            case anchorwave::AlignmentCandidate::LowWfa:
                status = anchorwave::alignWithLowWfa(
                        dna_q, dna_d,
                        memorySchedulingEnabled
                        ? memoryEnvelope.algorithmBudgetBytes
                        : wfaMemoryBudget,
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2, wfaResult,
                        executionOptions);
                break;
            default:
                break;
        }
        // The final live allocation can be much smaller than the peak that
        // triggered a cooperative memory abort (especially medium/low WFA).
        // Scheduling calibration must use the peak, not the post-reap value.
        actualMemoryBytes = std::max(wfaResult.memoryUsedBytes,
                                     wfaResult.memoryPeakBytes);
        if (status == anchorwave::WfaAlignmentStatus::Completed) {
            totalScore = wfaResult.score;
            _alignment_q = std::move(wfaResult.queryAlignment);
            _alignment_d = std::move(wfaResult.referenceAlignment);
            _alignMethod = anchorwave::alignmentCandidateName(candidate);
            statusText = "completed";
            return true;
        }
        statusText = status == anchorwave::WfaAlignmentStatus::MemoryLimit
                     ? "memory_limit"
                     : status == anchorwave::WfaAlignmentStatus::TimeLimit
                       ? "time_limit" : "failed";
        anchorwave::recordExactAlignmentRuntimeFailure(
                status == anchorwave::WfaAlignmentStatus::MemoryLimit);
        return false;
    };

    const auto executeExactCandidate = [&, prepared](
            const RuntimeCandidate &runtime,
            uint64_t reservedBytes) {
        prepared->lastExactAttemptMemoryLimited = false;
        prepared->lastExactAttemptMemoryBytes = 0;
        if (runtime.candidate ==
                    anchorwave::AlignmentCandidate::Ksw2ScoreCertified &&
            !prepared->certifiedOptimalScoreReady) {
            prepared->certifiedOptimalScore = 0;
            prepared->certifiedFinalBandWidth = -1;
            prepared->certifiedTracebackAttempts = 0;
            prepared->retainedBandedScore =
                    std::numeric_limits<int64_t>::min();
            prepared->retainedBandedBandWidth = -1;
            prepared->retainedBandedQueryAlignment.clear();
            prepared->retainedBandedReferenceAlignment.clear();
        }
        uint64_t actualMemoryBytes = 0;
        std::string statusText = "failed";
        const uint64_t residentBefore =
                anchorwave::currentProcessResidentBytes();
        const anchorwave::AnchorTaskExecutor::LoadSnapshot calibrationLoad =
                anchorwave::currentAnchorTaskLoadSnapshot();
        const auto started = std::chrono::steady_clock::now();
        CandidateMemoryEnvelope startedEnvelope = runtime.memoryEnvelope;
        startedEnvelope.reservationBytes = reservedBytes;
        recordAttempt(runtime.candidate, runtime.estimate, startedEnvelope,
                      runtime.predictedWaitMinutes, 0.0, 0, "started");
        bool completed = false;
        const double maximumRuntimeSeconds =
                remainingExactRuntimeSeconds();
        if (balancedTimeLimitEnabled &&
            maximumRuntimeSeconds <= 0.0) {
            prepared->exactTimeLimitReached = true;
            recordAttempt(runtime.candidate, runtime.estimate,
                          startedEnvelope,
                          runtime.predictedWaitMinutes, 0.0, 0,
                          "time_limit");
            return false;
        }
        try {
            if (isWfaCandidate(runtime.candidate)) {
                completed = attemptWfaCandidate(
                        runtime.candidate, runtime.memoryEnvelope,
                        maximumRuntimeSeconds,
                        actualMemoryBytes, statusText);
            } else if (runtime.candidate ==
                       anchorwave::AlignmentCandidate::
                               Ksw2ScoreCertified) {
                if (!prepared->certifiedOptimalScoreReady) {
                    ScopedKswDeadline deadline(maximumRuntimeSeconds);
                    bool scoreStopped = false;
                    const int64_t optimalScore = ksw2ScoreOnly(
                            dna_q, dna_d, mismatchingPenalty,
                            openGapPenalty1, extendGapPenalty1,
                            openGapPenalty2, extendGapPenalty2,
                            scoreStopped);
                    if (scoreStopped || deadline.stopped()) {
                        statusText = "time_limit";
                    } else if (optimalScore == KSW_NEG_INF) {
                        statusText = "failed";
                        anchorwave::recordExactAlignmentRuntimeFailure(false);
                    } else {
                        prepared->certifiedOptimalScore = optimalScore;
                        prepared->certifiedOptimalScoreReady = true;
                        statusText = "score_ready";
                    }
                } else {
                    anchorwave::Ksw2CertifiedResult certified =
                            alignScoreCertifiedKsw2(
                                    dna_q, dna_d,
                                    _alignment_q, _alignment_d,
                                    selection.ksw2CertifiedInitialBandWidth,
                                    selection.ksw2CertifiedMaximumBandWidth,
                                    maximumRuntimeSeconds,
                                    mismatchingPenalty,
                                    openGapPenalty1, extendGapPenalty1,
                                    openGapPenalty2, extendGapPenalty2,
                                    &selection,
                                    &prepared->certifiedOptimalScore);
                prepared->certifiedOptimalScore = certified.optimalScore;
                prepared->certifiedFinalBandWidth =
                        certified.finalBandWidth;
                prepared->certifiedTracebackAttempts +=
                        certified.tracebackAttempts;
                if (certified.status ==
                    anchorwave::Ksw2CertifiedStatus::Completed) {
                    totalScore = certified.score;
                    _alignMethod = anchorwave::alignmentCandidateName(
                            runtime.candidate);
                    completed = true;
                    statusText = "completed";
                } else if (certified.status ==
                           anchorwave::Ksw2CertifiedStatus::TimeLimit) {
                    statusText = "time_limit";
                } else if (certified.status ==
                           anchorwave::Ksw2CertifiedStatus::NotCertified) {
                    if (!certified.bestEffortQueryAlignment.empty() ||
                        !certified.bestEffortReferenceAlignment.empty()) {
                        prepared->retainedBandedScore = certified.score;
                        prepared->retainedBandedBandWidth =
                                certified.finalBandWidth;
                        prepared->retainedBandedQueryAlignment = std::move(
                                certified.bestEffortQueryAlignment);
                        prepared->retainedBandedReferenceAlignment = std::move(
                                certified.bestEffortReferenceAlignment);
                    }
                    statusText = "not_certified";
                    anchorwave::recordExactAlignmentRuntimeFailure(false);
                } else {
                    statusText = "failed";
                    anchorwave::recordExactAlignmentRuntimeFailure(false);
                }
                }
            } else if (runtime.candidate ==
                       anchorwave::AlignmentCandidate::Ksw2Singletrack) {
                ScopedKswDeadline deadline(maximumRuntimeSeconds);
                bool stopped = false;
                totalScore = alignKsw2Singletrack(
                        dna_q, dna_d, _alignment_q, _alignment_d,
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2,
                        &stopped);
                if (stopped || deadline.stopped()) {
                    _alignment_q.clear();
                    _alignment_d.clear();
                    statusText = "time_limit";
                } else if (totalScore == KSW_NEG_INF) {
                    statusText = "failed";
                    anchorwave::recordExactAlignmentRuntimeFailure(false);
                } else {
                    _alignMethod = anchorwave::alignmentCandidateName(
                            runtime.candidate);
                    completed = true;
                    statusText = "completed";
                }
            } else if (runtime.candidate ==
                       anchorwave::AlignmentCandidate::Ksw2Full) {
                ScopedKswDeadline deadline(maximumRuntimeSeconds);
                totalScore = alignSlidingWindow_minimap2(
                        dna_q, dna_d, _length_of_q, _length_of_d,
                        _alignment_q, _alignment_d, -1,
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2);
                if (deadline.stopped()) {
                    _alignment_q.clear();
                    _alignment_d.clear();
                    statusText = "time_limit";
                } else {
                    _alignMethod = anchorwave::alignmentCandidateName(
                            runtime.candidate);
                    completed = true;
                    statusText = "completed";
                }
            }
        } catch (const std::bad_alloc &) {
            statusText = "memory_limit";
            anchorwave::recordExactAlignmentRuntimeFailure(true);
        } catch (...) {
            const double seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - started).count();
            recordAttempt(runtime.candidate, runtime.estimate,
                          runtime.memoryEnvelope,
                          runtime.predictedWaitMinutes, seconds,
                          actualMemoryBytes, "exception");
            throw;
        }
        const uint64_t residentAfter =
                anchorwave::currentProcessResidentBytes();
        if (actualMemoryBytes == 0 && residentAfter > residentBefore) {
            actualMemoryBytes = residentAfter - residentBefore;
        }
        const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
        prepared->exactRuntimeSpentSeconds += seconds;
        if (statusText == "score_ready" &&
            runtime.candidate == anchorwave::AlignmentCandidate::
                                         Ksw2ScoreCertified) {
            prepared->certifiedScoreRuntimeSeconds = seconds;
        }
        if (statusText == "time_limit" ||
            (balancedTimeLimitEnabled &&
             prepared->exactRuntimeSpentSeconds >=
                     selection.provenance.
                             exactAlignmentMaximumEstimatedMinutes * 60.0)) {
            prepared->exactTimeLimitReached = true;
        }
        if (completed && runtime.candidate !=
                                 anchorwave::AlignmentCandidate::
                                         Ksw2ScoreCertified) {
            anchorwave::recordAlignmentRuntimeObservation(
                    runtime.candidate,
                    runtime.estimate->estimatedMinutes,
                    runtime.estimate->estimatedMinutesP90,
                    seconds,
                    calibrationLoad.activeTasks,
                    calibrationLoad.workerCount);
        }
        CandidateMemoryEnvelope tracedEnvelope = runtime.memoryEnvelope;
        tracedEnvelope.reservationBytes = reservedBytes;
        recordAttempt(runtime.candidate, runtime.estimate, tracedEnvelope,
                      runtime.predictedWaitMinutes, seconds,
                      actualMemoryBytes, statusText);
        prepared->lastExactAttemptMemoryLimited =
                statusText == "memory_limit";
        prepared->lastExactAttemptMemoryBytes = actualMemoryBytes;
        return completed;
    };

    const auto scheduleExpandedSingletrackRetry = [&, prepared](
            const RuntimeCandidate &runtime) {
        if (runtime.candidate !=
                    anchorwave::AlignmentCandidate::SingletrackWfa ||
            !prepared->lastExactAttemptMemoryLimited) {
            return false;
        }
        prepared->singletrackMemoryModelInvalidated = true;
        if (!memorySchedulingEnabled ||
            prepared->exactTimeLimitReached) {
            return false;
        }
        const uint64_t retryBudget =
                anchorwave::expandedSingletrackRetryBudgetBytes(
                        runtime.memoryEnvelope.algorithmBudgetBytes,
                        prepared->lastExactAttemptMemoryBytes,
                        wfaMemoryBudget,
                        prepared->singletrackExpandedRetryUsed);
        if (retryBudget == 0) {
            return false;
        }
        prepared->singletrackExpandedRetryUsed = true;
        prepared->singletrackExpandedRetryPending = true;
        prepared->attempted[candidateIndex(runtime.candidate)] = false;
        CandidateMemoryEnvelope retryEnvelope;
        retryEnvelope.algorithmBudgetBytes = retryBudget;
        retryEnvelope.reservationBytes = saturatedAdd(
                addFractionalMargin(retryBudget, 115, 100),
                transientMemoryBytes);
        recordAttempt(runtime.candidate, runtime.estimate, retryEnvelope,
                      0.0, 0.0, prepared->lastExactAttemptMemoryBytes,
                      "memory_retry_planned");
        return true;
    };

    // Tier 1 is a hard quality boundary.  Within it, rank every exact engine
    // by expected finish plus a dynamic memory-time opportunity cost derived
    // from the user-provided -t/-M resources and current backlog.  At the true
    // critical tail the shadow price becomes zero, so the fastest exact engine
    // may wait for a large memory reservation instead of starting a slow mode.
    // A preferred request can still be blocked by persistent post-task RSS
    // after tracked reservations drain. Keep that distinct from structural
    // infeasibility and compare the remaining exact modes before parking or
    // retrying the task.
    std::array<bool, anchorwave::kAlignmentCandidateCount>
            runtimePressureBlocked{};
    std::array<bool, anchorwave::kAlignmentCandidateCount>
            runtimeFloorBlocked{};
    bool finalAvailabilityRecheckPerformed = false;
    bool sawTemporaryExactPressure = false;
    RuntimeCandidate temporarilyBlockedExactCandidate;
    while (true) {
        if (prepared->exactTimeLimitReached) {
            break;
        }
        const anchorwave::AlignmentResourceSnapshot resources =
                anchorwave::alignmentResourceSnapshot();
        const anchorwave::AnchorTaskExecutor::LoadSnapshot load =
                anchorwave::currentAnchorTaskLoadSnapshot();
        const double remainingRuntimeSeconds =
                remainingExactRuntimeSeconds();
        if (balancedTimeLimitEnabled &&
            remainingRuntimeSeconds <= 0.0) {
            prepared->exactTimeLimitReached = true;
            break;
        }
        std::vector<RuntimeCandidate> ranked;
        ranked.reserve(selection.exactCandidates.size());
        for (const anchorwave::AlignmentCandidate candidate :
                selection.exactCandidates) {
            if (prepared->singletrackMemoryModelInvalidated &&
                anchorwave::sharesSingletrackMemoryUnderpredictionRisk(
                        candidate)) {
                continue;
            }
            if (prepared->attempted[candidateIndex(candidate)]) {
                continue;
            }
            if (runtimePressureBlocked[candidateIndex(candidate)]) {
                continue;
            }
            if (runtimeFloorBlocked[candidateIndex(candidate)]) {
                continue;
            }
            const anchorwave::AlgorithmCostEstimate *estimate =
                    candidateEstimate(selection, candidate);
            if (estimate == nullptr) {
                prepared->attempted[candidateIndex(candidate)] = true;
                continue;
            }
            RuntimeCandidate runtime;
            runtime.candidate = candidate;
            runtime.estimate = estimate;
            runtime.effectiveEstimate =
                    anchorwave::calibratedAlignmentCostEstimate(
                            *estimate, load.activeTasks, load.workerCount);
            if (candidate == anchorwave::AlignmentCandidate::
                                     Ksw2ScoreCertified &&
                prepared->certifiedOptimalScoreReady) {
                const double spentMinutes =
                        prepared->certifiedScoreRuntimeSeconds / 60.0;
                runtime.effectiveEstimate.estimatedMinutes = std::max(
                        0.0,
                        runtime.effectiveEstimate.estimatedMinutes -
                                spentMinutes);
                runtime.effectiveEstimate.estimatedMinutesP90 = std::max(
                        runtime.effectiveEstimate.estimatedMinutes,
                        runtime.effectiveEstimate.estimatedMinutesP90 -
                                spentMinutes);
            }
            const bool scoreOnlyStage =
                    candidate == anchorwave::AlignmentCandidate::
                                         Ksw2ScoreCertified &&
                    !prepared->certifiedOptimalScoreReady;
            const CandidateMemoryEnvelope completionMemoryEnvelope =
                    candidateMemoryEnvelope(
                            selection, candidate, wfaMemoryBudget,
                            transientMemoryBytes, 0);
            runtime.memoryEnvelope = candidateMemoryEnvelope(
                    selection, candidate, wfaMemoryBudget,
                    transientMemoryBytes,
                    scoreOnlyStage
                    ? selection.ksw2CertifiedScoreOnlyMemoryBytes
                    : 0);
            runtime.expandedSingletrackRetry =
                    candidate ==
                            anchorwave::AlignmentCandidate::SingletrackWfa &&
                    prepared->singletrackExpandedRetryPending;
            if (runtime.expandedSingletrackRetry) {
                runtime.memoryEnvelope.algorithmBudgetBytes = wfaMemoryBudget;
                runtime.memoryEnvelope.reservationBytes = saturatedAdd(
                        addFractionalMargin(wfaMemoryBudget, 115, 100),
                        transientMemoryBytes);
            }
            if (scoreOnlyStage &&
                anchorwave::estimateAlignmentMemoryAvailability(
                        completionMemoryEnvelope.reservationBytes, nullptr).
                                availability ==
                        anchorwave::MemoryAvailability::StaticInfeasible) {
                prepared->attempted[candidateIndex(candidate)] = true;
                anchorwave::recordExactAlignmentRuntimeFailure(true);
                continue;
            }
            // A previous exact failure may have consumed most of this
            // interval's balanced-mode budget. Do not start another expensive
            // exact engine when its admission quantiles no longer fit the
            // remaining wall time merely to have the cooperative deadline
            // abort it. Exact-first (-bt 0) never enters this branch.
            if (balancedTimeLimitEnabled &&
                !anchorwave::exactCandidateWithinTimeLimit(
                        runtime.effectiveEstimate.estimatedMinutes,
                        runtime.effectiveEstimate.estimatedMinutesP90,
                        remainingRuntimeSeconds / 60.0)) {
                prepared->attempted[candidateIndex(candidate)] = true;
                recordAttempt(candidate, estimate, runtime.memoryEnvelope,
                              0.0, 0.0, 0,
                              "remaining_time_infeasible");
                continue;
            }
            runtime.memoryAvailability =
                    anchorwave::estimateAlignmentMemoryAvailability(
                            runtime.memoryEnvelope.reservationBytes,
                            &prepared->memoryPressure[
                                    candidateIndex(candidate)]);
            runtime.predictedWaitMinutes =
                    runtime.memoryAvailability.waitMinutes;
            if (runtime.memoryAvailability.availability ==
                anchorwave::MemoryAvailability::StableRuntimeFloor) {
                // A completed WFA attempt may leave a per-worker allocator
                // cache in RSS after its scheduler token is released. Reclaim
                // this worker's caches once before classifying the floor.
                anchorwave::releaseUnusedAlignmentMemoryToSystem();
                runtime.memoryAvailability =
                        anchorwave::estimateAlignmentMemoryAvailability(
                                runtime.memoryEnvelope.reservationBytes,
                                &prepared->memoryPressure[
                                        candidateIndex(candidate)]);
                runtime.predictedWaitMinutes =
                        runtime.memoryAvailability.waitMinutes;
            }
            if (runtime.memoryAvailability.availability ==
                anchorwave::MemoryAvailability::StaticInfeasible) {
                prepared->attempted[candidateIndex(candidate)] = true;
                anchorwave::recordExactAlignmentRuntimeFailure(true);
                continue;
            }
            if (runtime.memoryAvailability.availability ==
                anchorwave::MemoryAvailability::StableRuntimeFloor) {
                runtimeFloorBlocked[candidateIndex(candidate)] = true;
                if (!prepared->runtimeFloorRecorded[
                            candidateIndex(candidate)]) {
                    recordAttempt(candidate, estimate,
                                  runtime.memoryEnvelope,
                                  runtime.predictedWaitMinutes,
                                  0.0, 0,
                                  "stable_runtime_rss_floor");
                    prepared->runtimeFloorRecorded[
                            candidateIndex(candidate)] = true;
                }
                continue;
            }
            ranked.push_back(runtime);
        }
        // On the active AVX512 implementation, the page-aware full-KSW2
        // traceback uses no more resident memory than KSW2-Singletrack, and
        // isolated B73/Mo17 pairs show lower P50 and P90 time for every tested
        // stratum.  Remove the epsilon-dominated candidate while full KSW2 is
        // still available. If full KSW2 fails at runtime, the next ranking
        // pass restores KSW2-Singletrack as an independent Tier-1 fallback.
        const auto fullKsw = std::find_if(
                ranked.begin(), ranked.end(),
                [](const RuntimeCandidate &runtime) {
                    return runtime.candidate ==
                            anchorwave::AlignmentCandidate::Ksw2Full;
                });
        const auto singletrackKsw = std::find_if(
                ranked.begin(), ranked.end(),
                [](const RuntimeCandidate &runtime) {
                    return runtime.candidate ==
                            anchorwave::AlignmentCandidate::Ksw2Singletrack;
                });
        constexpr double kKswImplementationReliabilityTolerance = 0.002;
        if (fullKsw != ranked.end() &&
            singletrackKsw != ranked.end() &&
            anchorwave::exactCandidateEpsilonDominates(
                    *fullKsw->estimate, *singletrackKsw->estimate,
                    kKswImplementationReliabilityTolerance)) {
            ranked.erase(singletrackKsw);
        }
        const bool fullKsw2RuntimeAvailable = std::any_of(
                ranked.begin(), ranked.end(),
                [](const RuntimeCandidate &runtime) {
                    return runtime.candidate ==
                            anchorwave::AlignmentCandidate::Ksw2Full;
                });
        if (anchorwave::deferScoreCertifiedKsw2BehindFullKsw2(
                    prepared->singletrackMemoryModelInvalidated,
                    fullKsw2RuntimeAvailable)) {
            ranked.erase(
                    std::remove_if(
                            ranked.begin(), ranked.end(),
                            [](const RuntimeCandidate &runtime) {
                                return runtime.candidate ==
                                        anchorwave::AlignmentCandidate::
                                                Ksw2ScoreCertified;
                            }),
                    ranked.end());
        }
        const bool independentKsw2RuntimeAvailable = std::any_of(
                ranked.begin(), ranked.end(),
                [](const RuntimeCandidate &runtime) {
                    return anchorwave::
                            isIndependentExactAfterSingletrackMemoryFailure(
                                    runtime.candidate);
                });
        if (independentKsw2RuntimeAvailable &&
            anchorwave::preferIndependentKsw2ForAmbiguousWfaProfile(
                    selection.profile)) {
            ranked.erase(
                    std::remove_if(
                            ranked.begin(), ranked.end(),
                            [](const RuntimeCandidate &runtime) {
                                return isWfaCandidate(runtime.candidate);
                            }),
                    ranked.end());
        }
        if (prepared->singletrackExpandedRetryPending) {
            const bool independentExactAvailable = std::any_of(
                    ranked.begin(), ranked.end(),
                    [](const RuntimeCandidate &runtime) {
                        return anchorwave::
                                isIndependentExactAfterSingletrackMemoryFailure(
                                        runtime.candidate);
                    });
            if (independentExactAvailable) {
                ranked.erase(
                        std::remove_if(
                                ranked.begin(), ranked.end(),
                                [](const RuntimeCandidate &runtime) {
                                    return runtime.expandedSingletrackRetry;
                                }),
                        ranked.end());
            }
        }
        if (ranked.empty()) {
            // StableRuntimeFloor is a bounded runtime conclusion, not a
            // startup-size rejection.  Probe every such exact candidate once
            // more at the tier boundary: an RSS drop that happened during
            // ranking must revive exact DP instead of opening banded/sliding.
            if (!finalAvailabilityRecheckPerformed) {
                finalAvailabilityRecheckPerformed = true;
                bool revivedExactCandidate = false;
                for (const anchorwave::AlignmentCandidate candidate :
                        selection.exactCandidates) {
                    const std::size_t index = candidateIndex(candidate);
                    if (!runtimeFloorBlocked[index] ||
                        prepared->attempted[index]) {
                        continue;
                    }
                    CandidateMemoryEnvelope memoryEnvelope =
                            candidateMemoryEnvelope(
                                    selection, candidate, wfaMemoryBudget,
                                    transientMemoryBytes,
                                    candidate ==
                                                    anchorwave::AlignmentCandidate::
                                                            Ksw2ScoreCertified &&
                                            !prepared->certifiedOptimalScoreReady
                                    ? selection.
                                              ksw2CertifiedScoreOnlyMemoryBytes
                                    : 0);
                    const anchorwave::MemoryAvailabilityEstimate availability =
                            anchorwave::estimateAlignmentMemoryAvailability(
                                    memoryEnvelope.reservationBytes,
                                    &prepared->memoryPressure[index]);
                    if (availability.availability !=
                                anchorwave::MemoryAvailability::
                                        StableRuntimeFloor &&
                        availability.availability !=
                                anchorwave::MemoryAvailability::
                                        StaticInfeasible) {
                        runtimeFloorBlocked[index] = false;
                        prepared->runtimeFloorRecorded[index] = false;
                        revivedExactCandidate = true;
                    }
                }
                if (revivedExactCandidate) {
                    continue;
                }
            }
            break;
        }
        // Solve the complete exact-attempt chain as a small stochastic
        // shortest-path problem. There are at most seven candidates, so the
        // Bellman subset DP is exact (O(n*2^n)). One version minimizes the
        // interval's expected finish; the resource-aware version also prices
        // memory displaced from the global backlog. This avoids ranking the
        // first engine with a resource-blind fallback chain.
        std::vector<anchorwave::ExactAttemptSystemCost> completionAttempts;
        std::vector<anchorwave::ExactAttemptSystemCost> systemAttempts;
        completionAttempts.reserve(ranked.size());
        systemAttempts.reserve(ranked.size());
        for (RuntimeCandidate &runtime : ranked) {
            runtime.memoryOpportunityMinutes = dynamicMemoryShadowMinutes(
                    runtime.effectiveEstimate.estimatedMinutesP90,
                    runtime.memoryEnvelope.reservationBytes,
                    resources, load);
            const double attemptMinutes = runtime.predictedWaitMinutes +
                    anchorwave::alignmentRiskAdjustedMinutes(
                            runtime.effectiveEstimate);
            completionAttempts.push_back(
                    anchorwave::ExactAttemptSystemCost{
                            runtime.candidate, attemptMinutes,
                            runtime.effectiveEstimate.successProbability});
            systemAttempts.push_back(anchorwave::ExactAttemptSystemCost{
                    runtime.candidate,
                    attemptMinutes + runtime.memoryOpportunityMinutes,
                    runtime.effectiveEstimate.successProbability});
        }
        for (std::size_t i = 0; i < ranked.size(); ++i) {
            ranked[i].expectedCompletionMinutes =
                    anchorwave::exactAttemptExpectedSystemMinutes(
                            completionAttempts, i);
            ranked[i].systemCost =
                    anchorwave::exactAttemptExpectedSystemMinutes(
                            systemAttempts, i);
        }
        double fastestExactExpectedMinutes =
                std::numeric_limits<double>::infinity();
        for (const RuntimeCandidate &runtime : ranked) {
            fastestExactExpectedMinutes = std::min(
                    fastestExactExpectedMinutes,
                    runtime.expectedCompletionMinutes -
                            runtime.predictedWaitMinutes);
        }
        for (RuntimeCandidate &runtime : ranked) {
            // Whichever exact engine is fastest for this interval enters the
            // same immediate-start lane. No algorithm name or memory mode is
            // privileged inside Tier 1.
            if (anchorwave::exactCandidateFastLaneEligible(
                        runtime.predictedWaitMinutes,
                        runtime.expectedCompletionMinutes -
                                runtime.predictedWaitMinutes,
                        fastestExactExpectedMinutes,
                        runtime.memoryOpportunityMinutes)) {
                // The fast lane waives only this immediately runnable
                // attempt's memory opportunity cost. Failure still enters the
                // Bellman-optimal remaining chain with its own resource costs.
                runtime.systemCost = std::max(
                        0.0, runtime.systemCost -
                                     runtime.memoryOpportunityMinutes);
            }
        }

        // Compare all waiting and runnable Tier-1 alternatives on the same
        // conservative bounds. A short wait is worthwhile only when the
        // waiting candidate's P90 still beats the ready candidate's P50.
        const RuntimeCandidate *bestWaiting = nullptr;
        const RuntimeCandidate *bestReady = nullptr;
        double bestWaitingBound = std::numeric_limits<double>::infinity();
        double bestReadyBound = std::numeric_limits<double>::infinity();
        double bestWaitingShadow = 0.0;
        double bestReadyShadow = 0.0;
        for (const RuntimeCandidate &runtime : ranked) {
            if (runtime.memoryAvailability.availability !=
                        anchorwave::MemoryAvailability::Ready &&
                std::isfinite(runtime.predictedWaitMinutes)) {
                const double p90 = std::max(
                        runtime.effectiveEstimate.estimatedMinutes,
                        runtime.effectiveEstimate.estimatedMinutesP90);
                const double shadow = dynamicMemoryShadowMinutes(
                        p90, runtime.memoryEnvelope.reservationBytes,
                        resources, load);
                const double bound = runtime.predictedWaitMinutes + p90 +
                                     shadow;
                if (bound < bestWaitingBound) {
                    bestWaiting = &runtime;
                    bestWaitingBound = bound;
                    bestWaitingShadow = shadow;
                }
            } else if (runtime.memoryAvailability.availability ==
                               anchorwave::MemoryAvailability::Ready) {
                const double p50 = std::max(
                        0.0, runtime.effectiveEstimate.estimatedMinutes);
                const double shadow = dynamicMemoryShadowMinutes(
                        p50, runtime.memoryEnvelope.reservationBytes,
                        resources, load);
                const double bound = p50 + shadow;
                if (bound < bestReadyBound) {
                    bestReady = &runtime;
                    bestReadyBound = bound;
                    bestReadyShadow = shadow;
                }
            }
        }
        anchorwave::AlignmentCandidate dominantWaitingCandidate =
                anchorwave::AlignmentCandidate::SlidingWindow;
        const bool holdForWaitingCandidate = bestWaiting != nullptr &&
                bestReady != nullptr &&
                anchorwave::waitingExactDominatesRunnableExact(
                        bestWaiting->effectiveEstimate,
                        bestWaiting->predictedWaitMinutes,
                        bestReady->effectiveEstimate,
                        bestWaitingShadow, bestReadyShadow);
        if (holdForWaitingCandidate) {
            dominantWaitingCandidate = bestWaiting->candidate;
        }
        std::stable_sort(
                ranked.begin(), ranked.end(),
                [](const RuntimeCandidate &left,
                   const RuntimeCandidate &right) {
                    if (left.systemCost != right.systemCost) {
                        return left.systemCost < right.systemCost;
                    }
                    return candidateIndex(left.candidate) <
                           candidateIndex(right.candidate);
                });
        if (holdForWaitingCandidate &&
            ranked.front().candidate != dominantWaitingCandidate) {
            const auto dominant = std::find_if(
                    ranked.begin(), ranked.end(),
                    [dominantWaitingCandidate](
                            const RuntimeCandidate &runtime) {
                        return runtime.candidate == dominantWaitingCandidate;
                    });
            if (dominant != ranked.end()) {
                std::rotate(ranked.begin(), dominant, dominant + 1);
            }
        }
        bool restartRanking = false;
        bool waiting = false;
        RuntimeCandidate waitingCandidate;
        for (const RuntimeCandidate &runtime : ranked) {
            anchorwave::AlignmentMemoryTryResult admission =
                    anchorwave::tryAcquireAlignmentMemory(
                            runtime.memoryEnvelope.reservationBytes,
                            runtime.effectiveEstimate.estimatedMinutesP90);
            if (admission.admission ==
                anchorwave::AlignmentMemoryAdmission::PermanentlyInfeasible) {
                prepared->attempted[candidateIndex(runtime.candidate)] = true;
                anchorwave::recordExactAlignmentRuntimeFailure(true);
                restartRanking = true;
                continue;
            }
            if (!admission) {
                waiting = true;
                waitingCandidate = runtime;
                break;
            }
            recordSelectionOnce();
            const bool scoreOnlyStage =
                    runtime.candidate == anchorwave::AlignmentCandidate::
                                                 Ksw2ScoreCertified &&
                    !prepared->certifiedOptimalScoreReady;
            prepared->attempted[candidateIndex(runtime.candidate)] = true;
            if (runtime.expandedSingletrackRetry) {
                prepared->singletrackExpandedRetryPending = false;
            }
            const uint64_t reservedBytes =
                    admission.reservation.reservedBytes();
            if (executeExactCandidate(runtime, reservedBytes)) {
                return totalScore;
            }
            scheduleExpandedSingletrackRetry(runtime);
            if (scoreOnlyStage && prepared->certifiedOptimalScoreReady &&
                !prepared->exactTimeLimitReached) {
                // The low-memory score phase is complete. Release its token
                // at the end of this scope, then re-enter admission for the
                // traceback stage without recomputing the exact score.
                prepared->attempted[candidateIndex(runtime.candidate)] = false;
            }
            restartRanking = true;
            break;
        }
        if (waiting) {
            const bool drainMemory =
                    anchorwave::currentAnchorTaskShouldDrainMemory();
            if (anchorwave::currentAnchorTaskCanDefer() && !drainMemory) {
                const uint64_t deferrals =
                        anchorwave::currentAnchorTaskDeferralCount();
                if (deferrals == 0 || (deferrals & (deferrals - 1)) == 0) {
                    recordAttempt(waitingCandidate.candidate,
                                  waitingCandidate.estimate,
                                  waitingCandidate.memoryEnvelope,
                                  waitingCandidate.predictedWaitMinutes,
                                  0.0, 0, "deferred_memory");
                }
                const double waitMilliseconds = std::isfinite(
                        waitingCandidate.predictedWaitMinutes)
                        ? waitingCandidate.predictedWaitMinutes * 30000.0
                        : 2000.0;
                throw anchorwave::AlignmentTaskDeferred(
                        std::chrono::milliseconds(
                                static_cast<int64_t>(std::max(
                                        250.0, std::min(
                                                5000.0,
                                                waitMilliseconds)))));
            }
            anchorwave::AlignmentMemoryTryResult preferred =
                    anchorwave::acquirePreferredAlignmentMemory(
                            waitingCandidate.memoryEnvelope.reservationBytes,
                            waitingCandidate.effectiveEstimate.
                                    estimatedMinutesP90);
            if (preferred) {
                recordSelectionOnce();
                const bool scoreOnlyStage =
                        waitingCandidate.candidate ==
                                anchorwave::AlignmentCandidate::
                                        Ksw2ScoreCertified &&
                        !prepared->certifiedOptimalScoreReady;
                prepared->attempted[candidateIndex(
                        waitingCandidate.candidate)] = true;
                if (waitingCandidate.expandedSingletrackRetry) {
                    prepared->singletrackExpandedRetryPending = false;
                }
                if (executeExactCandidate(
                            waitingCandidate,
                            preferred.reservation.reservedBytes())) {
                    return totalScore;
                }
                scheduleExpandedSingletrackRetry(waitingCandidate);
                if (scoreOnlyStage &&
                    prepared->certifiedOptimalScoreReady &&
                    !prepared->exactTimeLimitReached) {
                    prepared->attempted[candidateIndex(
                            waitingCandidate.candidate)] = false;
                }
                restartRanking = true;
            } else if (preferred.admission ==
                       anchorwave::AlignmentMemoryAdmission::
                               PermanentlyInfeasible) {
                prepared->attempted[candidateIndex(
                        waitingCandidate.candidate)] = true;
                anchorwave::recordExactAlignmentRuntimeFailure(true);
                restartRanking = true;
            } else {
                runtimePressureBlocked[candidateIndex(
                        waitingCandidate.candidate)] = true;
                sawTemporaryExactPressure = true;
                temporarilyBlockedExactCandidate = waitingCandidate;
                recordAttempt(waitingCandidate.candidate,
                              waitingCandidate.estimate,
                              waitingCandidate.memoryEnvelope,
                              waitingCandidate.predictedWaitMinutes,
                              0.0, 0, "persistent_rss_unavailable");
                restartRanking = true;
            }
        }
        if (!restartRanking && !waiting) {
            break;
        }
    }

    if (sawTemporaryExactPressure) {
        // Current RSS/queue occupancy is not evidence that an exact engine is
        // structurally impossible.  Never cross the tier-1 boundary for that
        // reason.  Deferrable genome tasks retry after consumers/caches have
        // released memory; legacy synchronous callers receive an explicit
        // error instead of silently returning a lower-quality alignment.
        if (anchorwave::currentAnchorTaskCanDefer()) {
            recordAttempt(temporarilyBlockedExactCandidate.candidate,
                          temporarilyBlockedExactCandidate.estimate,
                          temporarilyBlockedExactCandidate.memoryEnvelope,
                          temporarilyBlockedExactCandidate.predictedWaitMinutes,
                          0.0, 0, "deferred_persistent_rss");
            throw anchorwave::AlignmentTaskDeferred(
                    std::chrono::milliseconds(1000));
        }
        throw std::runtime_error(
                "exact alignment is temporarily blocked by current process "
                "RSS inside the configured -M envelope");
    }

    // Tier two contains both approximate strategies. The selector predicts
    // which method will produce the larger score and places it first. Run only
    // that candidate when it succeeds; the next candidate is a failure rescue,
    // not a mandatory second alignment.
    std::vector<anchorwave::AlignmentCandidate> fallbackCandidates =
            selection.approximateCandidates;
    if (fallbackCandidates.empty()) {
        fallbackCandidates = selection.bandedCandidates;
        fallbackCandidates.insert(
                fallbackCandidates.end(),
                selection.lastResortCandidates.begin(),
                selection.lastResortCandidates.end());
    }
    for (anchorwave::AlignmentCandidate candidate : fallbackCandidates) {
        if (candidate != anchorwave::AlignmentCandidate::Ksw2Banded &&
            candidate != anchorwave::AlignmentCandidate::SlidingWindow) {
            continue;
        }
        const CandidateMemoryEnvelope memoryEnvelope =
                candidateMemoryEnvelope(selection, candidate,
                                        wfaMemoryBudget,
                                        transientMemoryBytes);
        const anchorwave::AlgorithmCostEstimate *const estimate =
                candidateEstimate(selection, candidate);
        const double predictedMinutes = estimate != nullptr
                ? anchorwave::alignmentRiskAdjustedMinutes(*estimate) : 0.0;
        anchorwave::AlignmentMemoryTryResult admission =
                anchorwave::tryAcquireAlignmentMemory(
                        memoryEnvelope.reservationBytes, predictedMinutes);
        const std::size_t approximateIndex = candidateIndex(candidate);
        if (!admission &&
            !prepared->allocatorReclaimAttempted[approximateIndex]) {
            const anchorwave::AlignmentResourceSnapshot resources =
                    anchorwave::alignmentResourceSnapshot();
            if (resources.enabled && resources.activeReservations == 0) {
                prepared->allocatorReclaimAttempted[approximateIndex] = true;
                anchorwave::releaseUnusedAlignmentMemoryToSystem();
                admission = anchorwave::tryAcquireAlignmentMemory(
                        memoryEnvelope.reservationBytes, predictedMinutes);
                recordAttempt(
                        candidate, estimate, memoryEnvelope, 0.0, 0.0, 0,
                        admission ? "allocator_reclaim_admitted"
                                  : "allocator_reclaim_incomplete");
            }
        }
        if (admission.admission ==
                anchorwave::AlignmentMemoryAdmission::
                        PermanentlyInfeasible) {
            continue;
        }
        if (!admission) {
            recordAttempt(candidate, estimate, memoryEnvelope, 0.0,
                          0.0, 0, "deferred_approximate_memory");
            if (anchorwave::currentAnchorTaskCanDefer()) {
                throw anchorwave::AlignmentTaskDeferred(
                        std::chrono::milliseconds(1000));
            }
            admission = anchorwave::acquirePreferredAlignmentMemory(
                    memoryEnvelope.reservationBytes, predictedMinutes);
            if (!admission) {
                throw std::runtime_error(
                        "approximate alignment is temporarily blocked by "
                        "current process RSS inside the configured -M "
                        "envelope");
            }
        }
        anchorwave::AlignmentMemoryReservation memoryReservation =
                std::move(admission.reservation);
        recordSelectionOnce();
        const uint64_t residentBefore =
                anchorwave::currentProcessResidentBytes();
        const auto started = std::chrono::steady_clock::now();
        if (candidate == anchorwave::AlignmentCandidate::Ksw2Banded) {
            anchorwave::recordBandedFallbackExecution();
            if (prepared->retainedBandedBandWidth ==
                        selection.ksw2BandWidth &&
                (!prepared->retainedBandedQueryAlignment.empty() ||
                 !prepared->retainedBandedReferenceAlignment.empty())) {
                totalScore = prepared->retainedBandedScore;
                _alignment_q = std::move(
                        prepared->retainedBandedQueryAlignment);
                _alignment_d = std::move(
                        prepared->retainedBandedReferenceAlignment);
            } else {
                totalScore = alignSlidingWindow_minimap2(
                        dna_q, dna_d, _length_of_q, _length_of_d,
                        _alignment_q, _alignment_d, selection.ksw2BandWidth,
                        mismatchingPenalty,
                        openGapPenalty1, extendGapPenalty1,
                        openGapPenalty2, extendGapPenalty2);
            }
            _alignMethod = anchorwave::alignmentCandidateName(candidate);
        } else {
            anchorwave::recordSlidingFallbackExecution();
            totalScore = alignSlidingWindow(
                    dna_q, dna_d, _length_of_q, _length_of_d,
                    _alignment_q, _alignment_d,
                    selection.slidingWindowWidth,
                    matchingScore, mismatchingPenalty,
                    openGapPenalty1, extendGapPenalty1,
                    openGapPenalty2, extendGapPenalty2);
            _alignMethod = "SLIDING_WINDOW";
        }
        // A failed score-certification attempt may already have produced a
        // valid banded result. Comparing that already-paid-for score does not
        // add a second alignment, so retain the larger measured score.
        if (candidate == anchorwave::AlignmentCandidate::SlidingWindow &&
            prepared->retainedBandedScore > totalScore &&
            (!prepared->retainedBandedQueryAlignment.empty() ||
             !prepared->retainedBandedReferenceAlignment.empty())) {
            totalScore = prepared->retainedBandedScore;
            _alignment_q = std::move(
                    prepared->retainedBandedQueryAlignment);
            _alignment_d = std::move(
                    prepared->retainedBandedReferenceAlignment);
            _alignMethod = anchorwave::alignmentCandidateName(
                    anchorwave::AlignmentCandidate::Ksw2Banded);
        }
        const uint64_t residentAfter =
                anchorwave::currentProcessResidentBytes();
        const double actualSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
        if (!gappedAlignmentReconstructs(_alignment_q, dna_q) ||
            !gappedAlignmentReconstructs(_alignment_d, dna_d)) {
            recordAttempt(candidate, estimate, memoryEnvelope, 0.0,
                          actualSeconds,
                          residentAfter > residentBefore
                          ? residentAfter - residentBefore : 0,
                          "invalid_reconstruction");
            _alignment_q.clear();
            _alignment_d.clear();
            continue;
        }
        recordAttempt(candidate, estimate, memoryEnvelope, 0.0,
                      actualSeconds,
                      residentAfter > residentBefore
                      ? residentAfter - residentBefore : 0,
                      "completed");
        return totalScore;
    }

    // makeAlignmentSelectionPlan always supplies the approximate tier. Retain
    // a defensive fallback in case a future selector violates that contract.
    const CandidateMemoryEnvelope memoryEnvelope = candidateMemoryEnvelope(
            selection, anchorwave::AlignmentCandidate::SlidingWindow,
            wfaMemoryBudget, transientMemoryBytes);
    anchorwave::AlignmentMemoryTryResult admission =
            anchorwave::tryAcquireAlignmentMemory(
                    memoryEnvelope.reservationBytes);
    if (!admission && admission.admission ==
            anchorwave::AlignmentMemoryAdmission::TemporarilyUnavailable) {
        if (anchorwave::currentAnchorTaskCanDefer()) {
            throw anchorwave::AlignmentTaskDeferred(
                    std::chrono::milliseconds(1000));
        }
        admission = anchorwave::acquirePreferredAlignmentMemory(
                memoryEnvelope.reservationBytes);
    }
    if (!admission) {
        throw std::runtime_error(
                admission.admission ==
                        anchorwave::AlignmentMemoryAdmission::
                                PermanentlyInfeasible
                ? "no fallback alignment fits inside the configured -M "
                  "memory envelope"
                : "fallback alignment is temporarily blocked by current "
                  "process RSS inside the configured -M memory envelope");
    }
    anchorwave::AlignmentMemoryReservation memoryReservation =
            std::move(admission.reservation);
    recordSelectionOnce();
    anchorwave::recordSlidingFallbackExecution();
    const uint64_t residentBefore = anchorwave::currentProcessResidentBytes();
    const auto started = std::chrono::steady_clock::now();
    totalScore = alignSlidingWindow(
            dna_q, dna_d, _length_of_q, _length_of_d,
            _alignment_q, _alignment_d,
            selection.slidingWindowWidth,
            matchingScore, mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2);
    _alignMethod = "SLIDING_WINDOW";
    const uint64_t residentAfter = anchorwave::currentProcessResidentBytes();
    recordAttempt(anchorwave::AlignmentCandidate::SlidingWindow,
                  candidateEstimate(
                          selection,
                          anchorwave::AlignmentCandidate::SlidingWindow),
                  memoryEnvelope, 0.0,
                  std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - started).count(),
                  residentAfter > residentBefore
                  ? residentAfter - residentBefore : 0,
                  "completed");
    return totalScore;
}

int64_t alignSlidingWindow(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d, std::string &_alignMethod,
                           const int64_t &slidingWindowSize,  const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    return alignSlidingWindow_local_wfa2_v2(dna_q, dna_d, _alignment_q, _alignment_d, _alignMethod, slidingWindowSize, matchingScore,
                                            mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                                            std::shared_ptr<const anchorwave::AlignmentSelectionPlan>());
}

int64_t alignSlidingWindow(
        std::string &dna_q, std::string &dna_d,
        std::string &alignment_q, std::string &alignment_d,
        std::string &alignMethod,
        const int64_t &slidingWindowSize, const int32_t &matchingScore,
        const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2,
        const std::shared_ptr<const anchorwave::AlignmentSelectionPlan>
                &precomputedSelection) {
    return alignSlidingWindow_local_wfa2_v2(
            dna_q, dna_d, alignment_q, alignment_d, alignMethod,
            slidingWindowSize, matchingScore,
            mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
            openGapPenalty2, extendGapPenalty2, precomputedSelection);
}

int64_t alignSlidingWindowNW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                             const int64_t &slidingWindowSize, const int32_t &matchingScore,
                             const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    // This entry point directly executes the sliding implementation from the
    // shared approximate tier.
    anchorwave::recordSlidingFallbackExecution();
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns begin
    bool flag_q_all_N = std::all_of(dna_q.begin(), dna_q.end(), [](char c) { return c == 'N'; });
    bool flag_d_all_N = std::all_of(dna_d.begin(), dna_d.end(), [](char c) { return c == 'N'; });

    if (flag_q_all_N || flag_d_all_N) {
        _alignment_q = dna_q;
        _alignment_d = dna_d;

        int32_t count_ = abs(_length_of_q - _length_of_d);
        if (_length_of_q < _length_of_d) {
            _alignment_q += std::string(count_, '-');
        }

        if (_length_of_d < _length_of_q) {
            _alignment_d += std::string(count_, '-');
        }

        return totalScore;
    }

    totalScore = alignSlidingWindow(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                    mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

    return totalScore;
}
