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

#include "affine2p_wavefront_display.h"

/*
 * Display
 */
#define AFFINE2P_WAVEFRONTS_PRINT_ELEMENT(wavefront,k) \
  /* Check limits */ \
  if (wavefront!=NULL && wavefront->lo <= k && k <= wavefront->hi) { \
    if (wavefront->offsets[k] >= 0) { \
      fprintf(stream,"[%2d]",(int)wavefront->offsets[k]); \
    } else { \
      fprintf(stream,"[  ]"); \
    } \
  } else { \
    fprintf(stream,"    "); \
  }
void affine2p_wavefronts_print_wavefronts_block(
    FILE* const stream,
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score_begin,
    const int score_end) {
  // Compute min/max k
  int s, max_k=0, min_k=0;
  for (s=score_begin;s<=score_end;++s) {
    affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[s];
    if (mwavefront != NULL) {
      max_k = MAX(max_k,mwavefront->hi);
      min_k = MIN(min_k,mwavefront->lo);
    }
    affine2p_wavefront_t* const i1wavefront = affine2p_wavefronts->i1wavefronts[s];
    if (i1wavefront != NULL) {
      max_k = MAX(max_k,i1wavefront->hi);
      min_k = MIN(min_k,i1wavefront->lo);
    }
    affine2p_wavefront_t* const i2wavefront = affine2p_wavefronts->i2wavefronts[s];
    if (i2wavefront != NULL) {
      max_k = MAX(max_k,i2wavefront->hi);
      min_k = MIN(min_k,i2wavefront->lo);
    }
    affine2p_wavefront_t* const d1wavefront = affine2p_wavefronts->d1wavefronts[s];
    if (d1wavefront != NULL) {
      max_k = MAX(max_k,d1wavefront->hi);
      min_k = MIN(min_k,d1wavefront->lo);
    }
    affine2p_wavefront_t* const d2wavefront = affine2p_wavefronts->d2wavefronts[s];
    if (d2wavefront != NULL) {
      max_k = MAX(max_k,d2wavefront->hi);
      min_k = MIN(min_k,d2wavefront->lo);
    }
  }
  // Headers
  fprintf(stream,">[SCORE %3d-%3d]\n",score_begin,score_end);
  fprintf(stream,"         ");
  for (s=score_begin;s<=score_end;++s) {
    fprintf(stream,"[ M][I1][I2][D1][D2] ");
  }
  fprintf(stream,"\n");
  fprintf(stream,"        ");
  for (s=score_begin;s<=score_end;++s) {
    fprintf(stream,"+--------------------");
  }
  fprintf(stream,"\n");
  // Traverse all diagonals
  int k;
  for (k=max_k;k>=min_k;k--) {
    fprintf(stream,"[k=%3d] ",k);
    // Traverse all scores
    for (s=score_begin;s<=score_end;++s) {
      fprintf(stream,"|");
      // Fetch wavefront
      affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[s];
      affine2p_wavefront_t* const i1wavefront = affine2p_wavefronts->i1wavefronts[s];
      affine2p_wavefront_t* const i2wavefront = affine2p_wavefronts->i2wavefronts[s];
      affine2p_wavefront_t* const d1wavefront = affine2p_wavefronts->d1wavefronts[s];
      affine2p_wavefront_t* const d2wavefront = affine2p_wavefronts->d2wavefronts[s];
      AFFINE2P_WAVEFRONTS_PRINT_ELEMENT(mwavefront,k);
      AFFINE2P_WAVEFRONTS_PRINT_ELEMENT(i1wavefront,k);
      AFFINE2P_WAVEFRONTS_PRINT_ELEMENT(i2wavefront,k);
      AFFINE2P_WAVEFRONTS_PRINT_ELEMENT(d1wavefront,k);
      AFFINE2P_WAVEFRONTS_PRINT_ELEMENT(d2wavefront,k);
    }
    fprintf(stream,"|\n");
  }
  // Headers
  fprintf(stream,"        ");
  for (s=score_begin;s<=score_end;++s) {
    fprintf(stream,"+--------------------");
  }
  fprintf(stream,"\n");
  fprintf(stream,"SCORE   ");
  for (s=score_begin;s<=score_end;++s) {
    fprintf(stream,"          %2d         ",s);
  }
  fprintf(stream,"\n");
}
void affine2p_wavefronts_print_wavefronts(
    FILE* const stream,
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int current_score) {
  // Print wavefronts by chunks
  const int display_block_wf = 5;
  int s;
  for (s=0;s<=current_score;s+=display_block_wf) {
    affine2p_wavefronts_print_wavefronts_block(stream,
        affine2p_wavefronts,s,MIN(s+display_block_wf,current_score));
  }
}


