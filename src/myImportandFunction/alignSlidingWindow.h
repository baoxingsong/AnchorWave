//
// Created by song on 8/5/18.
//

#pragma once

#include "../cpu_arch.h"

#ifdef __SSE2NEON__
#include "../../sse2neon.h"
#else
#include <immintrin.h>
#endif // __SSE2NEON__

#include "../impl/getSequencesFromGff.h"
#include "../model//Score.h"
#include "../util/myutil.h"

#include "../../minimap2/ksw2.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>

namespace anchorwave {
struct AlignmentSelectionPlan;

enum class Ksw2CertifiedStatus {
    Completed,
    TimeLimit,
    NotCertified,
    Failed
};

struct Ksw2CertifiedResult {
    Ksw2CertifiedStatus status = Ksw2CertifiedStatus::Failed;
    int64_t score = 0;
    int64_t optimalScore = 0;
    int64_t finalBandWidth = -1;
    uint64_t tracebackAttempts = 0;
    // When certification exhausts its band budget, retain the final valid
    // banded traceback as a best-effort fallback. Public output arguments
    // remain empty unless status is Completed.
    std::string bestEffortQueryAlignment;
    std::string bestEffortReferenceAlignment;
};
}

// Exact, end-to-end 2-piece affine KSW2 alignment. This lower-level entry
// point is exposed so the exact-algorithm consistency tests can verify that
// KSW2 and every WFA memory mode optimize the same score.
int64_t alignSlidingWindow_minimap2(
        const std::string &dna_q, const std::string &dna_d,
        int64_t length_of_q, int64_t length_of_d,
        std::string &alignment_q, std::string &alignment_d,
        const int64_t &bandWidth,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2);

// Compute the exact unbanded KSW2 score without traceback storage, then grow
// a banded traceback geometrically until its score certifies that the emitted
// alignment is globally optimal. maximumRuntimeSeconds==0 disables the
// cooperative deadline. A non-completed result never writes an uncertified
// alignment to the output arguments; its result object may retain the final
// valid banded traceback for the production fallback selector.
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
        const anchorwave::AlignmentSelectionPlan *selectionPlan = nullptr,
        const int64_t *precomputedOptimalScore = nullptr);

int64_t alignSlidingWindow(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d, std::string &_alignMethod,
                           const int64_t &slidingWindowSize, const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2);

// Precompute the selector plan without retaining either sequence. The gap
// scheduler uses the resulting compact plan to rank long work globally, and
// the alignment attempt reuses it so profiling is not repeated.
std::shared_ptr<const anchorwave::AlignmentSelectionPlan>
prepareAlignmentSelectionPlan(
        const std::string &dna_q, const std::string &dna_d,
        const int64_t &slidingWindowSize,
        const int32_t &mismatchingPenalty,
        const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1,
        const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2);

int64_t alignSlidingWindow(
        std::string &dna_q, std::string &dna_d,
        std::string &alignment_q, std::string &alignment_d,
        std::string &alignMethod,
        const int64_t &slidingWindowSize, const int32_t &matchingScore,
        const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
        const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2,
        const int32_t &extendGapPenalty2,
        const std::shared_ptr<const anchorwave::AlignmentSelectionPlan>
                &precomputedSelection);

int64_t alignSlidingWindowNW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                             const int64_t &slidingWindowSize, const int32_t &matchingScore,
                             const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2);
