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
 * DESCRIPTION: Gap-affine 2-pieces alignment algorithms wrapper (including WFA)
 */

#include "../benchmark/benchmark_gap_affine2p.h"

#include "../gap_affine2p/affine2p_matrix.h"
#include "../gap_affine2p/affine2p_dp.h"
#include "../gap_affine2p/affine2p_wavefront.h"
#include "../gap_affine2p/affine2p_wavefront_align.h"

/*
 * Benchmark Gap-Affine 2-Pieces
 */
void benchmark_gap_affine2p_dp(
    align_input_t* const align_input,
    affine2p_penalties_t* const penalties) {
  // Allocate
  affine2p_matrix_t matrix;
  affine2p_matrix_allocate(
      &matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t cigar;
  cigar_allocate(&cigar,
      align_input->pattern_length+align_input->text_length,
      align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  affine2p_dp_compute(&matrix,penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,&cigar);
  timer_stop(&align_input->timer);
  // Free
  affine2p_matrix_free(&matrix,align_input->mm_allocator);
  cigar_free(&cigar);
}
void benchmark_gap_affine2p_wavefront(
    align_input_t* const align_input,
    affine2p_penalties_t* const penalties) {
  // Allocate
  affine2p_wavefronts_t* const affine2p_wavefronts =
      affine2p_wavefronts_new_complete(
          align_input->pattern_length,align_input->text_length,
          penalties,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  affine2p_wavefronts_align(affine2p_wavefronts,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length);
  timer_stop(&align_input->timer);
  // Free
  affine2p_wavefronts_delete(affine2p_wavefronts);
}
