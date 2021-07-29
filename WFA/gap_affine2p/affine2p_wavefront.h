/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WaveFront-Alignment for pairwise gap-affine 2-pieces alignment
 */

#ifndef AFFINE2P_WAVEFRONT_H_
#define AFFINE2P_WAVEFRONT_H_

#include "../utils/commons.h"
#include "../system/profiler_counter.h"
#include "../system/profiler_timer.h"
#include "../system/mm_allocator.h"
#include "../system/mm_stack.h"
#include "../alignment/cigar.h"
#include "affine2p_penalties.h"

/*
 * Constants
 */
#define AFFINE2P_WAVEFRONT_OFFSET_NULL (INT32_MIN)

/*
 * Translate k and offset to coordinates h,v
 */
#define AFFINE2P_WAVEFRONT_V(k,offset) ((offset)-(k))
#define AFFINE2P_WAVEFRONT_H(k,offset) (offset)

#define AFFINE2P_WAVEFRONT_DIAGONAL(h,v) ((h)-(v))
#define AFFINE2P_WAVEFRONT_OFFSET(h,v)   (h)

/*
 * Offset size
 */
typedef int32_t awf2p_offset_t;

/*
 * Wavefront
 */
typedef struct {
  // Range
  bool null;                  // Is null interval?
  int lo;                     // Effective lowest diagonal (inclusive)
  int hi;                     // Effective highest diagonal (inclusive)
  int lo_base;                // Lowest diagonal before reduction (inclusive)
  int hi_base;                // Highest diagonal before reduction (inclusive)
  // Offsets
  awf2p_offset_t* offsets;    // Offsets
} affine2p_wavefront_t;

/*
 * Wavefront Set
 */
typedef struct {
  /* In Wavefronts*/
  affine2p_wavefront_t* in_mwavefront_sub;
  affine2p_wavefront_t* in_mwavefront_gap1;
  affine2p_wavefront_t* in_mwavefront_gap2;
  affine2p_wavefront_t* in_i1wavefront_ext;
  affine2p_wavefront_t* in_i2wavefront_ext;
  affine2p_wavefront_t* in_d1wavefront_ext;
  affine2p_wavefront_t* in_d2wavefront_ext;
  /* Out Wavefronts */
  affine2p_wavefront_t* out_mwavefront;
  affine2p_wavefront_t* out_i1wavefront;
  affine2p_wavefront_t* out_i2wavefront;
  affine2p_wavefront_t* out_d1wavefront;
  affine2p_wavefront_t* out_d2wavefront;
} affine2p_wavefront_set_t;

/*
 * Gap-Affine 2-pieces Wavefronts
 */
typedef struct {
  // Dimensions
  int pattern_length;                          // Pattern length
  int text_length;                             // Text length
  int num_wavefronts;                          // Total number of allocatable wavefronts
  int max_allocated_wavefront;                 // Maximum index/score of allocated wavefront
  // Wavefronts
  affine2p_wavefront_t** mwavefronts;          // M-wavefronts
  affine2p_wavefront_t** i1wavefronts;         // I1-wavefronts
  affine2p_wavefront_t** i2wavefronts;         // I2-wavefronts
  affine2p_wavefront_t** d1wavefronts;         // D1-wavefronts
  affine2p_wavefront_t** d2wavefronts;         // D2-wavefronts
  affine2p_wavefront_t wavefront_null;         // Null wavefront (orthogonal reading)
  affine2p_wavefront_t wavefront_victim;       // Null wavefront (orthogonal writing)
//  // Reduction
//  affine2p_wavefronts_reduction_t reduction;   // Reduction parameters // TODO
  // Penalties
  affine2p_penalties_t penalties;              // Penalties parameters
  // CIGAR
  cigar_t cigar;                               // Alignment CIGAR
  // MM
  bool mm_allocator_own;                       // Ownership of MM-Allocator
  mm_allocator_t* mm_allocator;                // MM-Allocator (General memory allocator)
  mm_stack_t* mm_stack;                        // MM-Stack (Specific fast malloc/free wavefronts' memory)
  affine2p_wavefront_t* wavefronts_mem;        // MM-Slab  (Specific fast malloc/free wavefronts-ptr => base)
  affine2p_wavefront_t* wavefronts_current;    // MM-Slab  (Specific fast malloc/free wavefronts-ptr => next)
} affine2p_wavefronts_t;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Setup
 */
affine2p_wavefronts_t *affine2p_wavefronts_new_complete(
        const int pattern_length,
        const int text_length,
        affine2p_penalties_t *const penalties,
        mm_allocator_t *const mm_allocator);

#ifdef __cplusplus
}
#endif

affine2p_wavefronts_t* affine2p_wavefronts_new_adaptive(
    const int pattern_length,
    const int text_length,
    affine2p_penalties_t* const penalties,
    const int min_wavefront_length,
    const int max_distance_threshold,
    mm_allocator_t* const mm_allocator);
void affine2p_wavefronts_clear(
    affine2p_wavefronts_t* const affine2p_wavefronts);

#ifdef __cplusplus
extern "C" {
#endif
void affine2p_wavefronts_delete(
    affine2p_wavefronts_t* const affine2p_wavefronts);
#ifdef __cplusplus
}
#endif
/*
 * Fetch & allocate wavefronts
 */
void affine2p_wavefront_set_fetch_input(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    affine2p_wavefront_set_t* const wavefront_set,
    const int score);
void affine2p_wavefront_set_allocate_output(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    affine2p_wavefront_set_t* const wavefront_set,
    const int score,
    const int lo_effective,
    const int hi_effective);

/*
 * Initial Condition
 */
void affine2p_wavefront_initialize(
    affine2p_wavefronts_t* const affine2p_wavefronts);

#endif /* AFFINE2P_WAVEFRONT_H_ */
