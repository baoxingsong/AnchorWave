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
#include "affine2p_wavefront_extend.h"

/*
 * Wavefront offset extension comparing characters
 */
void affine2p_wavefronts_extend_mwavefront_compute_packed(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Fetch m-wavefront
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront==NULL) return;
  // Extend diagonally each wavefront point
  awf2p_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    // Fetch offset & positions
    const awf2p_offset_t offset = offsets[k];
    const int v = AFFINE2P_WAVEFRONT_V(k,offset);
    const int h = AFFINE2P_WAVEFRONT_H(k,offset);
    if (v < 0 || v >= pattern_length) continue;
    if (h < 0 || h >= text_length) continue;
    // Fetch pattern/text blocks
    uint64_t* pattern_blocks = (uint64_t*)(pattern+v);
    uint64_t* text_blocks = (uint64_t*)(text+h);
    uint64_t pattern_block = *pattern_blocks;
    uint64_t text_block = *text_blocks;
    // Compare 64-bits blocks
    uint64_t cmp = pattern_block ^ text_block;
    while (__builtin_expect(!cmp,0)) {
      // Increment offset (full block)
      offsets[k] += 8;
      // Next blocks
      ++pattern_blocks;
      ++text_blocks;
      // Fetch
      pattern_block = *pattern_blocks;
      text_block = *text_blocks;
      // Compare
      cmp = pattern_block ^ text_block;
    }
    // Count equal characters
    const int equal_right_bits = __builtin_ctzl(cmp);
    const int equal_chars = DIV_FLOOR(equal_right_bits,8);
    // Increment offset
    offsets[k] += equal_chars;
  }
}
void affine2p_wavefronts_extend_mwavefront_compute(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Fetch m-wavefront
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront==NULL) return;
  // Extend diagonally each wavefront point
  awf2p_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    // Exact extend
    const awf2p_offset_t offset = offsets[k];
    int v = AFFINE2P_WAVEFRONT_V(k,offset);
    int h = AFFINE2P_WAVEFRONT_H(k,offset);
    if (v <= 0 || v >= pattern_length) continue;
    if (h <= 0 || h >= text_length) continue;
    while (pattern[v++]==text[h++]) {
      ++(offsets[k]);
    }
  }
}
/*
 * Gap-Affine Wavefront exact extension
 */
void affine2p_wavefronts_extend_wavefront_packed(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Extend wavefront
  affine2p_wavefronts_extend_mwavefront_compute_packed(
      affine2p_wavefronts,pattern,pattern_length,
      text,text_length,score);
}
