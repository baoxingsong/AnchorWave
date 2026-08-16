#pragma once

// KSW2 Singletrack integration adapted from LorienLV/singletrack commit
// 7185195cb4049fd5290875ab4fc503384c4891dd (KSW2 base commit
// 289609bd9e5381a13b16239d0a7703f1ff03f9ca). See
// KSW2_SINGLETRACK_NOTICE for the upstream KSW2 MIT license.

#include "ksw2.h"

#ifdef __cplusplus
extern "C" {
#endif

// Exact two-piece affine global DP. The optimized score recurrence is KSW2;
// Singletrack stores the horizontal and vertical score differences needed
// to reconstruct an optimal path instead of KSW2's packed transition matrix.
void ksw_extd2_singletrack_sse(
        void *km, int qlen, const uint8_t *query,
        int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
        int8_t q, int8_t e, int8_t q2, int8_t e2,
        int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

// Shared width-independent traceback used by the SSE, AVX2 and AVX512
// recurrence kernels after they store the two Singletrack difference tracks.
void ksw2_singletrack_backtrace_affine2p(
        void *km, int qlen, const uint8_t *query,
        int tlen, const uint8_t *target, int8_t alphabetSize,
        const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
        int n_col, const int *off, const int *offEnd,
        const int8_t *p, const int8_t *p2, ksw_extz_t *ez);

#ifdef __AVX2__
void ksw_extd2_singletrack_avx2(
        void *km, int qlen, const uint8_t *query,
        int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
        int8_t q, int8_t e, int8_t q2, int8_t e2,
        int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif

#ifdef __AVX512BW__
void ksw_extd2_singletrack_avx512(
        void *km, int qlen, const uint8_t *query,
        int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
        int8_t q, int8_t e, int8_t q2, int8_t e2,
        int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif

#ifdef __cplusplus
}
#endif
