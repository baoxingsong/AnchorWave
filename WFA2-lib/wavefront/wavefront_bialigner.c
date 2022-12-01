/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
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
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "wavefront_bialigner.h"
#include "wavefront_aligner.h"
#include "wavefront_attributes.h"
#include "wavefront_heuristic.h"

/*
 * Setup
 */
wavefront_bialigner_t* wavefront_bialigner_new(
    wavefront_aligner_attr_t* const attributes,
    wavefront_plot_t* const plot) {
  // Allocate
  wavefront_bialigner_t* const wf_bialigner = malloc(sizeof(wavefront_bialigner_t));
  // Configure subsidiary aligners
  wavefront_aligner_attr_t subsidiary_attr = wavefront_aligner_attr_default;
  // Inherit attributes from master aligner
  subsidiary_attr.distance_metric = attributes->distance_metric;
  subsidiary_attr.linear_penalties = attributes->linear_penalties;
  subsidiary_attr.affine_penalties = attributes->affine_penalties;
  subsidiary_attr.affine2p_penalties = attributes->affine2p_penalties;
  subsidiary_attr.match_funct = attributes->match_funct;
  subsidiary_attr.match_funct_arguments = attributes->match_funct_arguments;
  // Set specifics for subsidiary aligners
  subsidiary_attr.heuristic = attributes->heuristic; // Inherit same heuristic
  subsidiary_attr.memory_mode = wavefront_memory_high; // Classic WFA
  subsidiary_attr.alignment_scope = compute_score;
  // Set other parameter for subsidiary aligners
  subsidiary_attr.system = attributes->system;
  // Allocate forward/reverse aligners
  wf_bialigner->alg_forward = wavefront_aligner_new(&subsidiary_attr);
  wf_bialigner->alg_forward->align_mode = wf_align_biwfa_breakpoint_forward;
  wf_bialigner->alg_forward->plot = plot;
  wf_bialigner->alg_reverse = wavefront_aligner_new(&subsidiary_attr);
  wf_bialigner->alg_reverse->align_mode = wf_align_biwfa_breakpoint_reverse;
  wf_bialigner->alg_reverse->plot = plot;
  // Allocate subsidiary aligner
  subsidiary_attr.alignment_scope = compute_alignment;
  wf_bialigner->alg_subsidiary = wavefront_aligner_new(&subsidiary_attr);
  wf_bialigner->alg_subsidiary->align_mode = wf_align_biwfa_subsidiary;
  wf_bialigner->alg_subsidiary->plot = plot;
  // Return
  return wf_bialigner;
}
void wavefront_bialigner_reap(
    wavefront_bialigner_t* const wf_bialigner) {
  wavefront_aligner_reap(wf_bialigner->alg_forward);
  wavefront_aligner_reap(wf_bialigner->alg_reverse);
  wavefront_aligner_reap(wf_bialigner->alg_subsidiary);
}
void wavefront_bialigner_delete(
    wavefront_bialigner_t* const wf_bialigner) {
  wavefront_aligner_delete(wf_bialigner->alg_forward);
  wavefront_aligner_delete(wf_bialigner->alg_reverse);
  wavefront_aligner_delete(wf_bialigner->alg_subsidiary);
  free(wf_bialigner);
}
/*
 * Accessors
 */
uint64_t wavefront_bialigner_get_size(
    wavefront_bialigner_t* const wf_bialigner) {
  return wavefront_aligner_get_size(wf_bialigner->alg_forward) +
      wavefront_aligner_get_size(wf_bialigner->alg_reverse) +
      wavefront_aligner_get_size(wf_bialigner->alg_subsidiary);
}
void wavefront_bialigner_set_heuristic(
    wavefront_bialigner_t* const wf_bialigner,
    wavefront_heuristic_t* const heuristic) {
  wf_bialigner->alg_forward->heuristic = *heuristic;
  wf_bialigner->alg_reverse->heuristic = *heuristic;
  wf_bialigner->alg_subsidiary->heuristic = *heuristic;
}
void wavefront_bialigner_set_match_funct(
    wavefront_bialigner_t* const wf_bialigner,
    int (*match_funct)(int,int,void*),
    void* const match_funct_arguments) {
  wf_bialigner->alg_forward->match_funct = match_funct;
  wf_bialigner->alg_forward->match_funct_arguments = match_funct_arguments;
  wf_bialigner->alg_reverse->match_funct = match_funct;
  wf_bialigner->alg_reverse->match_funct_arguments = match_funct_arguments;
  wf_bialigner->alg_subsidiary->match_funct = match_funct;
  wf_bialigner->alg_subsidiary->match_funct_arguments = match_funct_arguments;
}
void wavefront_bialigner_set_max_alignment_score(
    wavefront_bialigner_t* const wf_bialigner,
    const int max_alignment_score) {
  wf_bialigner->alg_forward->system.max_alignment_score = max_alignment_score;
  wf_bialigner->alg_reverse->system.max_alignment_score = max_alignment_score;
  wf_bialigner->alg_subsidiary->system.max_alignment_score = max_alignment_score;
}
void wavefront_bialigner_set_max_memory(
    wavefront_bialigner_t* const wf_bialigner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort) {
  wf_bialigner->alg_forward->system.max_memory_resident = max_memory_resident;
  wf_bialigner->alg_forward->system.max_memory_abort = max_memory_abort;
  wf_bialigner->alg_reverse->system.max_memory_resident = max_memory_resident;
  wf_bialigner->alg_reverse->system.max_memory_abort = max_memory_abort;
  wf_bialigner->alg_subsidiary->system.max_memory_resident = max_memory_resident;
  wf_bialigner->alg_subsidiary->system.max_memory_abort = max_memory_abort;
}
void wavefront_bialigner_set_max_num_threads(
        wavefront_bialigner_t* const wf_bialigner,
        const int max_num_threads) {
    wf_bialigner->alg_forward->system.max_num_threads = max_num_threads;
    wf_bialigner->alg_reverse->system.max_num_threads = max_num_threads;
    wf_bialigner->alg_subsidiary->system.max_num_threads = max_num_threads;
}
void wavefront_bialigner_set_min_offsets_per_thread(
        wavefront_bialigner_t* const wf_bialigner,
        const int min_offsets_per_thread) {
    wf_bialigner->alg_forward->system.min_offsets_per_thread = min_offsets_per_thread;
    wf_bialigner->alg_reverse->system.min_offsets_per_thread = min_offsets_per_thread;
    wf_bialigner->alg_subsidiary->system.min_offsets_per_thread = min_offsets_per_thread;
}
