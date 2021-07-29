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

#include "affine2p_penalties.h"
#include "affine2p_wavefront.h"

/*
 * Setup
 */
void affine2p_wavefronts_allocate_null(
    affine2p_wavefronts_t* const affine2p_wavefronts) {
  // Parameters
  const int wavefront_length = 2*affine2p_wavefronts->num_wavefronts + 1;
  // Allocate null wavefront
  awf2p_offset_t* const offsets_null = mm_allocator_calloc(
      affine2p_wavefronts->mm_allocator,wavefront_length,awf2p_offset_t,false);
  // Initialize
  affine2p_wavefronts->wavefront_null.null = true;
  affine2p_wavefronts->wavefront_null.lo =  1;
  affine2p_wavefronts->wavefront_null.hi = -1;
  affine2p_wavefronts->wavefront_null.lo_base =  1;
  affine2p_wavefronts->wavefront_null.hi_base = -1;
  affine2p_wavefronts->wavefront_null.offsets = offsets_null + affine2p_wavefronts->num_wavefronts; // Center at k=0
  int i;
  for (i=0;i<wavefront_length;++i) {
    offsets_null[i] = AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
void affine2p_wavefronts_allocate_victim(
    affine2p_wavefronts_t* const affine2p_wavefronts) {
  // Parameters
  const int wavefront_length = 2*affine2p_wavefronts->num_wavefronts + 1;
  // Allocate null wavefront
  awf2p_offset_t* const offsets_victim = mm_allocator_calloc(
      affine2p_wavefronts->mm_allocator,wavefront_length,awf2p_offset_t,false);
  // Initialize
  affine2p_wavefronts->wavefront_victim.null = true;
  affine2p_wavefronts->wavefront_victim.lo =  1;
  affine2p_wavefronts->wavefront_victim.hi = -1;
  affine2p_wavefronts->wavefront_victim.lo_base =  1;
  affine2p_wavefronts->wavefront_victim.hi_base = -1;
  affine2p_wavefronts->wavefront_victim.offsets = offsets_victim + affine2p_wavefronts->num_wavefronts; // Center at k=0
}
void affine2p_wavefronts_allocate_components(
    affine2p_wavefronts_t* const affine2p_wavefronts) {
  // Parameters
  mm_allocator_t* const mm_allocator = affine2p_wavefronts->mm_allocator;
  // Initialize wavefronts
  const int num_wavefronts = affine2p_wavefronts->num_wavefronts;
  affine2p_wavefronts->mwavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine2p_wavefront_t*,true);
  affine2p_wavefronts->i1wavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine2p_wavefront_t*,true);
  affine2p_wavefronts->i2wavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine2p_wavefront_t*,true);
  affine2p_wavefronts->d1wavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine2p_wavefront_t*,true);
  affine2p_wavefronts->d2wavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine2p_wavefront_t*,true);
  // Allocate bulk-memory (for all wavefronts)
  affine2p_wavefront_t* const wavefronts_mem =
      mm_allocator_calloc(mm_allocator,5*num_wavefronts,affine2p_wavefront_t,false);
  affine2p_wavefronts->wavefronts_mem = wavefronts_mem;
  affine2p_wavefronts->wavefronts_current = wavefronts_mem;
}
affine2p_wavefronts_t* affine2p_wavefronts_new(
    const int pattern_length,
    const int text_length,
    affine2p_penalties_t* const penalties,
    mm_allocator_t* mm_allocator) {
  // MM
  bool mm_allocator_own = false;
  if (mm_allocator == NULL) {
    mm_allocator = mm_allocator_new(BUFFER_SIZE_4M); // For handlers, ptr-vectors, etc
    mm_allocator_own = true;
  }
  mm_stack_t* const mm_stack = mm_stack_new(BUFFER_SIZE_4M); // For WF-Offsets
  // Handler
  affine2p_wavefronts_t* const affine2p_wavefronts = mm_allocator_alloc(mm_allocator,affine2p_wavefronts_t);
  affine2p_wavefronts->mm_allocator = mm_allocator;
  affine2p_wavefronts->mm_allocator_own = mm_allocator_own;
  affine2p_wavefronts->mm_stack = mm_stack;
  // Dimensions
  const int max_seq_length = ABS(pattern_length-text_length);
  const int max_score_misms = MIN(pattern_length,text_length) * penalties->mismatch;
  const int max_score_indel1 = penalties->gap_opening1 + max_seq_length * penalties->gap_extension1;
  const int max_score_indel2 = penalties->gap_opening2 + max_seq_length * penalties->gap_extension2;
  const int max_score_indel = MAX(max_score_indel1,max_score_indel2);
  const int num_wavefronts = max_score_misms + max_score_indel;
  affine2p_wavefronts->pattern_length = pattern_length;
  affine2p_wavefronts->text_length = text_length;
  affine2p_wavefronts->num_wavefronts = num_wavefronts;
  affine2p_wavefronts->max_allocated_wavefront = 0;
  // Penalties
  affine2p_wavefronts->penalties = *penalties;
  affine2p_wavefronts->penalties.match = 0; // TODO
  // Allocate wavefronts
  affine2p_wavefronts_allocate_components(affine2p_wavefronts);
  affine2p_wavefronts_allocate_null(affine2p_wavefronts);
  affine2p_wavefronts_allocate_victim(affine2p_wavefronts);
  // CIGAR
  cigar_allocate(&affine2p_wavefronts->cigar,pattern_length+text_length,mm_allocator);
  // Return
  return affine2p_wavefronts;
}

#ifdef __cplusplus
extern "C" {
#endif
affine2p_wavefronts_t* affine2p_wavefronts_new_complete(
    const int pattern_length,
    const int text_length,
    affine2p_penalties_t* const penalties,
    mm_allocator_t* const mm_allocator) {
  // Create new
  affine2p_wavefronts_t* const affine2p_wavefronts =
      affine2p_wavefronts_new(pattern_length,text_length,penalties,mm_allocator);
//  // Reduction
//  affine2p_wavefronts_reduction_set_none(&affine2p_wavefronts->reduction); // TODO
  // Return
  return affine2p_wavefronts;
}
#ifdef __cplusplus
}
#endif

affine2p_wavefronts_t* affine2p_wavefronts_new_adaptive(
    const int pattern_length,
    const int text_length,
    affine2p_penalties_t* const penalties,
    const int min_wavefront_length,
    const int max_distance_threshold,
    mm_allocator_t* const mm_allocator) {
  // Create new
  affine2p_wavefronts_t* const affine2p_wavefronts =
      affine2p_wavefronts_new(pattern_length,text_length,penalties,mm_allocator);
//  // Reduction
//  affine2p_wavefronts_reduction_set_adaptive(
//      &affine2p_wavefronts->reduction,min_wavefront_length,max_distance_threshold); // TODO
  // Return
  return affine2p_wavefronts;
}
void affine2p_wavefronts_clear(
    affine2p_wavefronts_t* const affine2p_wavefronts) {
  // Clear wavefronts memory
  const int num_wavefronts = MIN(affine2p_wavefronts->max_allocated_wavefront,affine2p_wavefronts->num_wavefronts);
  memset(affine2p_wavefronts->mwavefronts,0,num_wavefronts*sizeof(affine2p_wavefront_t*));
  memset(affine2p_wavefronts->i1wavefronts,0,num_wavefronts*sizeof(affine2p_wavefront_t*));
  memset(affine2p_wavefronts->i2wavefronts,0,num_wavefronts*sizeof(affine2p_wavefront_t*));
  memset(affine2p_wavefronts->d1wavefronts,0,num_wavefronts*sizeof(affine2p_wavefront_t*));
  memset(affine2p_wavefronts->d2wavefronts,0,num_wavefronts*sizeof(affine2p_wavefront_t*));
  mm_stack_clear(affine2p_wavefronts->mm_stack);
  // Clear CIGAR
  cigar_clear(&affine2p_wavefronts->cigar);
  // Clear wavefronts-ptr
  affine2p_wavefronts->wavefronts_current = affine2p_wavefronts->wavefronts_mem;
}
void affine2p_wavefronts_delete(
    affine2p_wavefronts_t* const affine2p_wavefronts) {
  // Parameters
  mm_allocator_t* const mm_allocator = affine2p_wavefronts->mm_allocator;
  // Free MID-Wavefronts
  mm_allocator_free(mm_allocator,affine2p_wavefronts->mwavefronts);
  mm_allocator_free(mm_allocator,affine2p_wavefronts->i1wavefronts);
  mm_allocator_free(mm_allocator,affine2p_wavefronts->i2wavefronts);
  mm_allocator_free(mm_allocator,affine2p_wavefronts->d1wavefronts);
  mm_allocator_free(mm_allocator,affine2p_wavefronts->d2wavefronts);
  mm_allocator_free(mm_allocator,affine2p_wavefronts->wavefront_null.offsets - affine2p_wavefronts->num_wavefronts);
  mm_allocator_free(mm_allocator,affine2p_wavefronts->wavefront_victim.offsets - affine2p_wavefronts->num_wavefronts);
  // Free bulk memory
  mm_allocator_free(mm_allocator,affine2p_wavefronts->wavefronts_mem);
  // CIGAR
  cigar_free(&affine2p_wavefronts->cigar);
  // DEBUG
#ifdef AFFINE2P_WAVEFRONT_DEBUG
  affine_matrix_free(&affine2p_wavefronts->affine_matrix,mm_allocator);
#endif
  // MM
  mm_stack_delete(affine2p_wavefronts->mm_stack);
  if (affine2p_wavefronts->mm_allocator_own) {
    mm_allocator_delete(affine2p_wavefronts->mm_allocator);
  }
  // Handler
  mm_allocator_free(mm_allocator,affine2p_wavefronts);
}
/*
 * Accessors
 */
affine2p_wavefront_t* affine2p_wavefronts_get_source_mwavefront(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score) {
  return (score < 0 || affine2p_wavefronts->mwavefronts[score] == NULL) ?
      &affine2p_wavefronts->wavefront_null : affine2p_wavefronts->mwavefronts[score];
}
affine2p_wavefront_t* affine2p_wavefronts_get_source_i1wavefront(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score) {
  return (score < 0 || affine2p_wavefronts->i1wavefronts[score] == NULL) ?
      &affine2p_wavefronts->wavefront_null : affine2p_wavefronts->i1wavefronts[score];
}
affine2p_wavefront_t* affine2p_wavefronts_get_source_i2wavefront(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score) {
  return (score < 0 || affine2p_wavefronts->i2wavefronts[score] == NULL) ?
      &affine2p_wavefronts->wavefront_null : affine2p_wavefronts->i2wavefronts[score];
}
affine2p_wavefront_t* affine2p_wavefronts_get_source_d1wavefront(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score) {
  return (score < 0 || affine2p_wavefronts->d1wavefronts[score] == NULL) ?
      &affine2p_wavefronts->wavefront_null : affine2p_wavefronts->d1wavefronts[score];
}
affine2p_wavefront_t* affine2p_wavefronts_get_source_d2wavefront(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score) {
  return (score < 0 || affine2p_wavefronts->d2wavefronts[score] == NULL) ?
      &affine2p_wavefronts->wavefront_null : affine2p_wavefronts->d2wavefronts[score];
}
/*
 * Fetch & allocate wavefronts
 */
affine2p_wavefront_t* affine2p_wavefront_allocate(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int lo_base,
    const int hi_base) {
  // Compute limits
  const int wavefront_length = hi_base - lo_base + 2; // (+1) for k=0
  // Allocate wavefront
  affine2p_wavefront_t* const wavefront = affine2p_wavefronts->wavefronts_current;
  ++(affine2p_wavefronts->wavefronts_current); // Next
  // Configure offsets
  wavefront->null = false;
  wavefront->lo = lo_base;
  wavefront->hi = hi_base;
  wavefront->lo_base = lo_base;
  wavefront->hi_base = hi_base;
  // Allocate offsets
  awf2p_offset_t* const offsets_mem = mm_stack_calloc(
      affine2p_wavefronts->mm_stack,wavefront_length,awf2p_offset_t,false);
  awf2p_offset_t* const offsets = offsets_mem - lo_base; // Center at k=0
  wavefront->offsets = offsets;
  // Return
  return wavefront;
}
void affine2p_wavefront_set_fetch_input(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    affine2p_wavefront_set_t* const wavefront_set,
    const int score) {
  // Compute scores
  const affine2p_penalties_t* const penalties = &(affine2p_wavefronts->penalties);
  const int mismatch = score - penalties->mismatch;
  const int gap_open1 = score - penalties->gap_opening1 - penalties->gap_extension1;
  const int gap_extend1 = score - penalties->gap_extension1;
  const int gap_open2 = score - penalties->gap_opening2 - penalties->gap_extension2;
  const int gap_extend2 = score - penalties->gap_extension2;
  // Fetch wavefronts
  wavefront_set->in_mwavefront_sub = affine2p_wavefronts_get_source_mwavefront(affine2p_wavefronts,mismatch);
  wavefront_set->in_mwavefront_gap1 = affine2p_wavefronts_get_source_mwavefront(affine2p_wavefronts,gap_open1);
  wavefront_set->in_mwavefront_gap2 = affine2p_wavefronts_get_source_mwavefront(affine2p_wavefronts,gap_open2);
  wavefront_set->in_i1wavefront_ext = affine2p_wavefronts_get_source_i1wavefront(affine2p_wavefronts,gap_extend1);
  wavefront_set->in_i2wavefront_ext = affine2p_wavefronts_get_source_i2wavefront(affine2p_wavefronts,gap_extend2);
  wavefront_set->in_d1wavefront_ext = affine2p_wavefronts_get_source_d1wavefront(affine2p_wavefronts,gap_extend1);
  wavefront_set->in_d2wavefront_ext = affine2p_wavefronts_get_source_d2wavefront(affine2p_wavefronts,gap_extend2);
}
void affine2p_wavefront_set_allocate_output(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    affine2p_wavefront_set_t* const wavefront_set,
    const int score,
    const int lo_effective,
    const int hi_effective) {
  // Allocate M-Wavefront
  wavefront_set->out_mwavefront = affine2p_wavefront_allocate(affine2p_wavefronts,lo_effective,hi_effective);
  affine2p_wavefronts->mwavefronts[score] = wavefront_set->out_mwavefront;
  // Allocate I-Wavefront
  if (!wavefront_set->in_mwavefront_gap1->null || !wavefront_set->in_i1wavefront_ext->null) {
    wavefront_set->out_i1wavefront = affine2p_wavefront_allocate(affine2p_wavefronts,lo_effective,hi_effective);
    affine2p_wavefronts->i1wavefronts[score] = wavefront_set->out_i1wavefront;
  } else {
    wavefront_set->out_i1wavefront = &affine2p_wavefronts->wavefront_victim;
  }
  if (!wavefront_set->in_mwavefront_gap2->null || !wavefront_set->in_i2wavefront_ext->null) {
    wavefront_set->out_i2wavefront = affine2p_wavefront_allocate(affine2p_wavefronts,lo_effective,hi_effective);
    affine2p_wavefronts->i2wavefronts[score] = wavefront_set->out_i2wavefront;
  } else {
    wavefront_set->out_i2wavefront = &affine2p_wavefronts->wavefront_victim;
  }
  // Allocate D-Wavefront
  if (!wavefront_set->in_mwavefront_gap1->null || !wavefront_set->in_d1wavefront_ext->null) {
    wavefront_set->out_d1wavefront = affine2p_wavefront_allocate(affine2p_wavefronts,lo_effective,hi_effective);
    affine2p_wavefronts->d1wavefronts[score] = wavefront_set->out_d1wavefront;
  } else {
    wavefront_set->out_d1wavefront = &affine2p_wavefronts->wavefront_victim;
  }
  if (!wavefront_set->in_mwavefront_gap2->null || !wavefront_set->in_d2wavefront_ext->null) {
    wavefront_set->out_d2wavefront = affine2p_wavefront_allocate(affine2p_wavefronts,lo_effective,hi_effective);
    affine2p_wavefronts->d2wavefronts[score] = wavefront_set->out_d2wavefront;
  } else {
    wavefront_set->out_d2wavefront = &affine2p_wavefronts->wavefront_victim;
  }
  // Increase max-wavefront
  affine2p_wavefronts->max_allocated_wavefront = MAX(affine2p_wavefronts->max_allocated_wavefront,score);
}
/*
 * Initial Conditions and finalization
 */
void affine2p_wavefront_initialize(
    affine2p_wavefronts_t* const affine2p_wavefronts) {
  affine2p_wavefronts->mwavefronts[0] = affine2p_wavefront_allocate(affine2p_wavefronts,0,0);
  affine2p_wavefronts->mwavefronts[0]->offsets[0] = 0;
}



