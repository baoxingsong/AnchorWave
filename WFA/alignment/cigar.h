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
 * DESCRIPTION: Cigar data-structure (match/mismatch/insertion/deletion)
 */

#ifndef CIGAR_H_
#define CIGAR_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"
#include "../gap_lineal/lineal_penalties.h"
#include "../gap_affine/affine_penalties.h"
#include "../gap_affine2p/affine2p_penalties.h"

/*
 * CIGAR
 */
typedef struct {
  // Operations buffer
  char* operations;
  int max_operations;
  int begin_offset;
  int end_offset;
  // Score
  int score;
  // MM
  mm_allocator_t* mm_allocator;
} cigar_t;

/*
 * Distance metrics
 */
typedef enum {
  edit,
  gap_lineal,
  gap_affine,
  gap_affine_2piece
} distance_metric_t;

/*
 * Setup
 */
void cigar_allocate(
    cigar_t* const cigar,
    const int max_operations,
    mm_allocator_t* const mm_allocator);
void cigar_clear(
    cigar_t* const cigar);
void cigar_free(
    cigar_t* const cigar);

/*
 * Accessors
 */
int cigar_get_matches(
    cigar_t* const cigar);
void cigar_add_mismatches(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    cigar_t* const cigar);

/*
 * Score
 */
int cigar_score_edit(
    cigar_t* const cigar);
int cigar_score_gap_lineal(
    cigar_t* const cigar,
    lineal_penalties_t* const penalties);
int cigar_score_gap_affine(
    cigar_t* const cigar,
    affine_penalties_t* const penalties);
#ifdef __cplusplus
extern "C" {
#endif
int cigar_score_gap_affine2p(
    cigar_t* const cigar,
    affine2p_penalties_t* const penalties);
#ifdef __cplusplus
}
#endif
/*
 * Utils
 */
int cigar_cmp(
    cigar_t* const cigar_a,
    cigar_t* const cigar_b);
void cigar_copy(
    cigar_t* const cigar_dst,
    cigar_t* const cigar_src);
bool cigar_check_alignment(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar,
    const bool verbose);

/*
 * Display
 */
void cigar_print(
    FILE* const stream,
    cigar_t* const cigar);
void cigar_print_pretty(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar,
    mm_allocator_t* const mm_allocator);

#endif /* CIGAR_H_ */
