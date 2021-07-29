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

#include "../utils/string_padded.h"
#include "affine2p_penalties.h"
#include "affine2p_wavefront_align.h"
#include "affine2p_wavefront_extend.h"
#include "affine2p_wavefront_backtrace.h"
#include "affine2p_wavefront_display.h"

/*
 * Compute wavefront-set limits
 */
void affine2p_wavefronts_compute_limits(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const affine2p_wavefront_set_t* const wavefront_set,
    const int score,
    int* const lo_effective,
    int* const hi_effective) {
  // Set limits (min_lo)
  int lo = wavefront_set->in_mwavefront_sub->lo;
  if (lo > wavefront_set->in_mwavefront_gap1->lo) lo = wavefront_set->in_mwavefront_gap1->lo;
  if (lo > wavefront_set->in_mwavefront_gap2->lo) lo = wavefront_set->in_mwavefront_gap2->lo;
  if (lo > wavefront_set->in_i1wavefront_ext->lo) lo = wavefront_set->in_i1wavefront_ext->lo;
  if (lo > wavefront_set->in_i2wavefront_ext->lo) lo = wavefront_set->in_i2wavefront_ext->lo;
  if (lo > wavefront_set->in_d1wavefront_ext->lo) lo = wavefront_set->in_d1wavefront_ext->lo;
  if (lo > wavefront_set->in_d2wavefront_ext->lo) lo = wavefront_set->in_d2wavefront_ext->lo;
  --lo;
  // Set limits (max_hi)
  int hi = wavefront_set->in_mwavefront_sub->hi;
  if (hi < wavefront_set->in_mwavefront_gap1->hi) hi = wavefront_set->in_mwavefront_gap1->hi;
  if (hi < wavefront_set->in_mwavefront_gap2->hi) hi = wavefront_set->in_mwavefront_gap2->hi;
  if (hi < wavefront_set->in_i1wavefront_ext->hi) hi = wavefront_set->in_i1wavefront_ext->hi;
  if (hi < wavefront_set->in_i2wavefront_ext->hi) hi = wavefront_set->in_i2wavefront_ext->hi;
  if (hi < wavefront_set->in_d1wavefront_ext->hi) hi = wavefront_set->in_d1wavefront_ext->hi;
  if (hi < wavefront_set->in_d2wavefront_ext->hi) hi = wavefront_set->in_d2wavefront_ext->hi;
  ++hi;
  // Set effective limits values
  *hi_effective = hi;
  *lo_effective = lo;
}
/*
 * Compute wavefront offsets
 */
#define AFFINE2P_WAVEFRONT_DECLARE(wavefront,prefix) \
  const awf2p_offset_t* const prefix ## _offsets = wavefront->offsets; \
  const int prefix ## _hi = wavefront->hi; \
  const int prefix ## _lo = wavefront->lo
#define AFFINE2P_WAVEFRONT_COND_FETCH(prefix,index,value) \
  (prefix ## _lo <= (index) && (index) <= prefix ## _hi) ? (value) : AFFINE2P_WAVEFRONT_OFFSET_NULL
/*
 * Compute next-wavefront (and kernel specializations)
 */
void affine2p_wavefronts_compute_next_idm(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const affine2p_wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap1,m_gap1);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap2,m_gap2);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_i1wavefront_ext,i1_ext);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_i2wavefront_ext,i2_ext);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_d1wavefront_ext,d1_ext);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_d2wavefront_ext,d2_ext);
  awf2p_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  awf2p_offset_t* const out_i1offsets = wavefront_set->out_i1wavefront->offsets;
  awf2p_offset_t* const out_i2offsets = wavefront_set->out_i2wavefront->offsets;
  awf2p_offset_t* const out_d1offsets = wavefront_set->out_d1wavefront->offsets;
  awf2p_offset_t* const out_d2offsets = wavefront_set->out_d2wavefront->offsets;
  // Compute loop peeling offset (min_hi)
  int min_hi = wavefront_set->in_mwavefront_sub->hi;
  if (!wavefront_set->in_mwavefront_gap1->null && min_hi > wavefront_set->in_mwavefront_gap1->hi-1) min_hi = wavefront_set->in_mwavefront_gap1->hi-1;
  if (!wavefront_set->in_mwavefront_gap2->null && min_hi > wavefront_set->in_mwavefront_gap2->hi-1) min_hi = wavefront_set->in_mwavefront_gap2->hi-1;
  if (!wavefront_set->in_i1wavefront_ext->null && min_hi > wavefront_set->in_i1wavefront_ext->hi+1) min_hi = wavefront_set->in_i1wavefront_ext->hi+1;
  if (!wavefront_set->in_i2wavefront_ext->null && min_hi > wavefront_set->in_i2wavefront_ext->hi+1) min_hi = wavefront_set->in_i2wavefront_ext->hi+1;
  if (!wavefront_set->in_d1wavefront_ext->null && min_hi > wavefront_set->in_d1wavefront_ext->hi-1) min_hi = wavefront_set->in_d1wavefront_ext->hi-1;
  if (!wavefront_set->in_d2wavefront_ext->null && min_hi > wavefront_set->in_d2wavefront_ext->hi-1) min_hi = wavefront_set->in_d2wavefront_ext->hi-1;
  // Compute loop peeling offset (max_lo)
  int max_lo = wavefront_set->in_mwavefront_sub->lo;
  if (!wavefront_set->in_mwavefront_gap1->null && max_lo < wavefront_set->in_mwavefront_gap1->lo+1) max_lo = wavefront_set->in_mwavefront_gap1->lo+1;
  if (!wavefront_set->in_mwavefront_gap2->null && max_lo < wavefront_set->in_mwavefront_gap2->lo+1) max_lo = wavefront_set->in_mwavefront_gap2->lo+1;
  if (!wavefront_set->in_i1wavefront_ext->null && max_lo < wavefront_set->in_i1wavefront_ext->lo+1) max_lo = wavefront_set->in_i1wavefront_ext->lo+1;
  if (!wavefront_set->in_i2wavefront_ext->null && max_lo < wavefront_set->in_i2wavefront_ext->lo+1) max_lo = wavefront_set->in_i2wavefront_ext->lo+1;
  if (!wavefront_set->in_d1wavefront_ext->null && max_lo < wavefront_set->in_d1wavefront_ext->lo-1) max_lo = wavefront_set->in_d1wavefront_ext->lo-1;
  if (!wavefront_set->in_d2wavefront_ext->null && max_lo < wavefront_set->in_d2wavefront_ext->lo-1) max_lo = wavefront_set->in_d2wavefront_ext->lo-1;
  // Compute score wavefronts (prologue)
  int k;
  for (k=lo;k<max_lo;++k) {
    // Update I1
    const awf2p_offset_t ins1_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap1,k-1,m_gap1_offsets[k-1]);
    const awf2p_offset_t ins1_i = AFFINE2P_WAVEFRONT_COND_FETCH(i1_ext,k-1,i1_ext_offsets[k-1]);
    const awf2p_offset_t ins1 = MAX(ins1_g,ins1_i) + 1;
    out_i1offsets[k] = ins1;
    // Update I2
    const awf2p_offset_t ins2_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap2,k-1,m_gap2_offsets[k-1]);
    const awf2p_offset_t ins2_i = AFFINE2P_WAVEFRONT_COND_FETCH(i2_ext,k-1,i2_ext_offsets[k-1]);
    const awf2p_offset_t ins2 = MAX(ins2_g,ins2_i) + 1;
    out_i2offsets[k] = ins2;
    // Update I
    const awf2p_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const awf2p_offset_t del1_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap1,k+1,m_gap1_offsets[k+1]);
    const awf2p_offset_t del1_d = AFFINE2P_WAVEFRONT_COND_FETCH(d1_ext,k+1,d1_ext_offsets[k+1]);
    const awf2p_offset_t del1 = MAX(del1_g,del1_d);
    out_d1offsets[k] = del1;
    // Update D2
    const awf2p_offset_t del2_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap2,k+1,m_gap2_offsets[k+1]);
    const awf2p_offset_t del2_d = AFFINE2P_WAVEFRONT_COND_FETCH(d2_ext,k+1,d2_ext_offsets[k+1]);
    const awf2p_offset_t del2 = MAX(del2_g,del2_d);
    out_d2offsets[k] = del2;
    // Update D
    const awf2p_offset_t del = MAX(del1,del2);
    // Update M
    const awf2p_offset_t sub = AFFINE2P_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(del,MAX(sub,ins));
  }
  // Compute score wavefronts (core)
#if defined(__clang__)
  #pragma clang loop vectorize(enable)
#elif defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=max_lo;k<=min_hi;++k) {
    // Update I1
    const awf2p_offset_t ins1_g = m_gap1_offsets[k-1];
    const awf2p_offset_t ins1_i = i1_ext_offsets[k-1];
    const awf2p_offset_t ins1 = MAX(ins1_g,ins1_i) + 1;
    out_i1offsets[k] = ins1;
    // Update I2
    const awf2p_offset_t ins2_g = m_gap2_offsets[k-1];
    const awf2p_offset_t ins2_i = i2_ext_offsets[k-1];
    const awf2p_offset_t ins2 = MAX(ins2_g,ins2_i) + 1;
    out_i2offsets[k] = ins2;
    // Update I
    const awf2p_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const awf2p_offset_t del1_g = m_gap1_offsets[k+1];
    const awf2p_offset_t del1_d = d1_ext_offsets[k+1];
    const awf2p_offset_t del1 = MAX(del1_g,del1_d);
    out_d1offsets[k] = del1;
    // Update D2
    const awf2p_offset_t del2_g = m_gap2_offsets[k+1];
    const awf2p_offset_t del2_d = d2_ext_offsets[k+1];
    const awf2p_offset_t del2 = MAX(del2_g,del2_d);
    out_d2offsets[k] = del2;
    // Update D
    const awf2p_offset_t del = MAX(del1,del2);
    // Update M
    const awf2p_offset_t sub = m_sub_offsets[k] + 1;
    out_moffsets[k] = MAX(del,MAX(sub,ins));
  }
  // Compute score wavefronts (epilogue)
  for (k=min_hi+1;k<=hi;++k) {
    // Update I1
    const awf2p_offset_t ins1_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap1,k-1,m_gap1_offsets[k-1]);
    const awf2p_offset_t ins1_i = AFFINE2P_WAVEFRONT_COND_FETCH(i1_ext,k-1,i1_ext_offsets[k-1]);
    const awf2p_offset_t ins1 = MAX(ins1_g,ins1_i) + 1;
    out_i1offsets[k] = ins1;
    // Update I2
    const awf2p_offset_t ins2_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap2,k-1,m_gap2_offsets[k-1]);
    const awf2p_offset_t ins2_i = AFFINE2P_WAVEFRONT_COND_FETCH(i2_ext,k-1,i2_ext_offsets[k-1]);
    const awf2p_offset_t ins2 = MAX(ins2_g,ins2_i) + 1;
    out_i2offsets[k] = ins2;
    // Update I
    const awf2p_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const awf2p_offset_t del1_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap1,k+1,m_gap1_offsets[k+1]);
    const awf2p_offset_t del1_d = AFFINE2P_WAVEFRONT_COND_FETCH(d1_ext,k+1,d1_ext_offsets[k+1]);
    const awf2p_offset_t del1 = MAX(del1_g,del1_d);
    out_d1offsets[k] = del1;
    // Update D2
    const awf2p_offset_t del2_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap2,k+1,m_gap2_offsets[k+1]);
    const awf2p_offset_t del2_d = AFFINE2P_WAVEFRONT_COND_FETCH(d2_ext,k+1,d2_ext_offsets[k+1]);
    const awf2p_offset_t del2 = MAX(del2_g,del2_d);
    out_d2offsets[k] = del2;
    // Update D
    const awf2p_offset_t del = MAX(del1,del2);
    // Update M
    const awf2p_offset_t sub = AFFINE2P_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(del,MAX(sub,ins));
  }
}
void affine2p_wavefronts_compute_next_im(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const affine2p_wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap1,m_gap1);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap2,m_gap2);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_i1wavefront_ext,i1_ext);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_i2wavefront_ext,i2_ext);
  awf2p_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  awf2p_offset_t* const out_i1offsets = wavefront_set->out_i1wavefront->offsets;
  awf2p_offset_t* const out_i2offsets = wavefront_set->out_i2wavefront->offsets;
  // Compute score wavefronts
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=lo;k<=hi;++k) {
    // Update I1
    const awf2p_offset_t ins1_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap1,k-1,m_gap1_offsets[k-1]);
    const awf2p_offset_t ins1_i = AFFINE2P_WAVEFRONT_COND_FETCH(i1_ext,k-1,i1_ext_offsets[k-1]);
    const awf2p_offset_t ins1 = MAX(ins1_g,ins1_i) + 1;
    out_i1offsets[k] = ins1;
    // Update I2
    const awf2p_offset_t ins2_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap2,k-1,m_gap2_offsets[k-1]);
    const awf2p_offset_t ins2_i = AFFINE2P_WAVEFRONT_COND_FETCH(i2_ext,k-1,i2_ext_offsets[k-1]);
    const awf2p_offset_t ins2 = MAX(ins2_g,ins2_i) + 1;
    out_i2offsets[k] = ins2;
    // Update I
    const awf2p_offset_t ins = MAX(ins1,ins2);
    // Update M
    const awf2p_offset_t sub = AFFINE2P_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(ins,sub);
  }
}
void affine2p_wavefronts_compute_next_dm(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const affine2p_wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap1,m_gap1);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap2,m_gap2);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_d1wavefront_ext,d1_ext);
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_d2wavefront_ext,d2_ext);
  awf2p_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  awf2p_offset_t* const out_d1offsets = wavefront_set->out_d1wavefront->offsets;
  awf2p_offset_t* const out_d2offsets = wavefront_set->out_d2wavefront->offsets;
  // Compute score wavefronts
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=lo;k<=hi;++k) {
    // Update D1
    const awf2p_offset_t del1_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap1,k+1,m_gap1_offsets[k+1]);
    const awf2p_offset_t del1_d = AFFINE2P_WAVEFRONT_COND_FETCH(d1_ext,k+1,d1_ext_offsets[k+1]);
    const awf2p_offset_t del1 = MAX(del1_g,del1_d);
    out_d1offsets[k] = del1;
    // Update D2
    const awf2p_offset_t del2_g = AFFINE2P_WAVEFRONT_COND_FETCH(m_gap2,k+1,m_gap2_offsets[k+1]);
    const awf2p_offset_t del2_d = AFFINE2P_WAVEFRONT_COND_FETCH(d2_ext,k+1,d2_ext_offsets[k+1]);
    const awf2p_offset_t del2 = MAX(del2_g,del2_d);
    out_d2offsets[k] = del2;
    // Update D
    const awf2p_offset_t del = MAX(del1,del2);
    // Update M
    const awf2p_offset_t sub = AFFINE2P_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(del,sub);
  }
}
void affine2p_wavefronts_compute_next_m(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const affine2p_wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE2P_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  awf2p_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  // Compute score wavefronts
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=lo;k<=hi;++k) {
    // Update M
    out_moffsets[k] = AFFINE2P_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
  }
}
/*
 * Compute wavefront
 */
void affine2p_wavefronts_compute_next(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Select wavefronts
  affine2p_wavefront_set_t wavefront_set;
  affine2p_wavefront_set_fetch_input(affine2p_wavefronts,&wavefront_set,score);
  // Check null wavefronts
  if (wavefront_set.in_mwavefront_sub->null &&
      wavefront_set.in_mwavefront_gap1->null &&
      wavefront_set.in_mwavefront_gap2->null &&
      wavefront_set.in_i1wavefront_ext->null &&
      wavefront_set.in_i2wavefront_ext->null &&
      wavefront_set.in_d1wavefront_ext->null &&
      wavefront_set.in_d2wavefront_ext->null) {
    return;
  }
  // Set limits
  int hi, lo;
  affine2p_wavefronts_compute_limits(affine2p_wavefronts,&wavefront_set,score,&lo,&hi);
  // Allocate score-wavefronts
  affine2p_wavefront_set_allocate_output(affine2p_wavefronts,&wavefront_set,score,lo,hi);
  // Compute WF
  affine2p_wavefronts_compute_next_idm(affine2p_wavefronts,&wavefront_set,lo,hi);
}
/*
 * Alignment end reached
 */
bool affine2p_wavefront_end_reached(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int pattern_length,
    const int text_length,
    const int score) {
  // Parameters
  const int alignment_k = AFFINE2P_WAVEFRONT_DIAGONAL(text_length,pattern_length);
  const int alignment_offset = AFFINE2P_WAVEFRONT_OFFSET(text_length,pattern_length);
  // Fetch wavefront and check termination
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront!=NULL) {
    awf2p_offset_t* const offsets = mwavefront->offsets;
    if (mwavefront->lo <= alignment_k &&
        alignment_k <= mwavefront->hi &&
        offsets[alignment_k] >= alignment_offset) {
      return true;
    }
  }
  return false;
}
/*
 * Computation using Wavefronts
 */
void affine2p_wavefronts_align(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Init padded strings
  strings_padded_t* const strings_padded =
      strings_padded_new_rhomb(
          pattern,pattern_length,text,text_length,
          AFFINE2P_WAVEFRONT_PADDING,affine2p_wavefronts->mm_allocator);
  // Initialize wavefront
  affine2p_wavefront_initialize(affine2p_wavefronts);
  // Compute wavefronts for increasing score
  int score = 0;
  while (true) {
    // Exact extend s-wavefront
    affine2p_wavefronts_extend_wavefront_packed(
        affine2p_wavefronts,strings_padded->pattern_padded,pattern_length,
        strings_padded->text_padded,text_length,score);
    // Exit condition
    if (affine2p_wavefront_end_reached(affine2p_wavefronts,pattern_length,text_length,score)) {
      // Backtrace & check alignment reached
      affine2p_wavefronts_backtrace(
          affine2p_wavefronts,strings_padded->pattern_padded,pattern_length,
          strings_padded->text_padded,text_length,score);
      break;
    }
    // Generate next wavefronts
    ++score; // Increase score
    affine2p_wavefronts_compute_next(
        affine2p_wavefronts,strings_padded->pattern_padded,pattern_length,
        strings_padded->text_padded,text_length,score);
  }
  // Free
  strings_padded_delete(strings_padded);
}


/*
 * Computation using Wavefronts, function customed by Baoxing
 */

int64_t affine2p_wavefronts_align_anchorwavenumberthreshold_version(
        affine2p_wavefronts_t* const affine2p_wavefronts,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length, const int64_t anchorwavenumberthreshold) {
//    printf("WFA line 419 \n");
    // Init padded strings
    strings_padded_t* const strings_padded =
            strings_padded_new_rhomb(
                    pattern,pattern_length,text,text_length,
                    AFFINE2P_WAVEFRONT_PADDING,affine2p_wavefronts->mm_allocator);
    // Initialize wavefront
    affine2p_wavefront_initialize(affine2p_wavefronts);
    // Compute wavefronts for increasing score
    int score = 0;
    int64_t numberFronts = 0;
    while (true) {

        {
            affine2p_wavefront_t *const mwavefront = affine2p_wavefronts->mwavefronts[score];
            if (mwavefront != NULL){
                numberFronts = numberFronts + mwavefront->hi - mwavefront->lo + 1;
                if (numberFronts > anchorwavenumberthreshold) {
                    strings_padded_delete(strings_padded);
                    return numberFronts;
                }
            }
        }

        // Exact extend s-wavefront
        affine2p_wavefronts_extend_wavefront_packed(
                affine2p_wavefronts,strings_padded->pattern_padded,pattern_length,
                strings_padded->text_padded,text_length,score);
        // Exit condition
        if (affine2p_wavefront_end_reached(affine2p_wavefronts,pattern_length,text_length,score)) {
            // Backtrace & check alignment reached
            affine2p_wavefronts_backtrace(
                    affine2p_wavefronts,strings_padded->pattern_padded,pattern_length,
                    strings_padded->text_padded,text_length,score);
            break;
        }
        // Generate next wavefronts
        ++score; // Increase score


        //The following 10 lines of source code, could be further optimized.
        // But, I that need to changed some other functions, at this stage, I do not want to change any function implemented by Santiago.
        // Select wavefronts
        affine2p_wavefront_set_t wavefront_set;
        affine2p_wavefront_set_fetch_input(affine2p_wavefronts,&wavefront_set,score); // this function did not change the value of affine2p_wavefronts
        int hi, lo;
        affine2p_wavefronts_compute_limits(affine2p_wavefronts,&wavefront_set,score,&lo,&hi);
        int count = hi-lo+1;
        numberFronts = numberFronts + count*5;
        if ( numberFronts > anchorwavenumberthreshold ){
            strings_padded_delete(strings_padded);
            return numberFronts;
        }

        affine2p_wavefronts_compute_next(
                affine2p_wavefronts,strings_padded->pattern_padded,pattern_length,
                strings_padded->text_padded,text_length,score);
    }
    // Free
    strings_padded_delete(strings_padded);
    return 0;
}
