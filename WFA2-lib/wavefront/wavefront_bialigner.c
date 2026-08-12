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

#include "utils/commons.h"
#include "wavefront_bialigner.h"
#include "wavefront_aligner.h"
#include "wavefront_attributes.h"
#include "wavefront_heuristic.h"

/* Report a high-WFA leaf's growth as aggregate BiWFA memory. */
static bool wavefront_bialigner_leaf_memory_probe(
    void* arguments,
    const uint64_t leaf_memory_used,
    const uint64_t leaf_memory_limit,
    const int score) {
  (void)leaf_memory_used;
  (void)leaf_memory_limit;
  wavefront_bialigner_t* const bialigner =
      (wavefront_bialigner_t*)arguments;
  wavefront_aligner_t* const master = bialigner->wf_master;
  bialigner->max_num_threads_used = MAX(
      bialigner->max_num_threads_used,
      bialigner->wf_base->align_status.max_num_threads_used);
  const uint64_t aggregate_memory = wavefront_bialigner_get_size(bialigner);
  master->align_status.memory_used = aggregate_memory;
  master->align_status.memory_peak =
      MAX(master->align_status.memory_peak,aggregate_memory);
  ++master->align_status.memory_probe_count;
  alignment_system_t* const system = &master->system;
  const bool callback_allowed = system->memory_probe == NULL ||
      system->memory_probe(system->memory_probe_arguments,
                           aggregate_memory,
                           system->biwfa_max_memory_abort,
                           score);
  return callback_allowed &&
      aggregate_memory <= system->biwfa_max_memory_abort;
}

/*
 * Setup
 */
wavefront_bialigner_t* wavefront_bialigner_new(
    wavefront_aligner_t* const wf_master,
    wavefront_aligner_attr_t* const attributes,
    wavefront_plot_t* const plot) {
  // Allocate
  wavefront_bialigner_t* const wf_bialigner = malloc(sizeof(wavefront_bialigner_t));
  wf_bialigner->wf_master = wf_master;
  // Configure subsidiary aligners
  wavefront_aligner_attr_t subsidiary_attr = wavefront_aligner_attr_default;
  // Inherit attributes from master aligner
  subsidiary_attr.distance_metric = attributes->distance_metric;
  subsidiary_attr.linear_penalties = attributes->linear_penalties;
  subsidiary_attr.affine_penalties = attributes->affine_penalties;
  subsidiary_attr.affine2p_penalties = attributes->affine2p_penalties;
  // Set specifics for subsidiary aligners
  subsidiary_attr.heuristic = attributes->heuristic;   // Inherit same heuristic
  subsidiary_attr.memory_mode = wavefront_memory_high; // Classic WFA
  subsidiary_attr.alignment_scope = compute_score;     // BiWFAs are score-only
  subsidiary_attr.alignment_form.extension = false;    // Deactivate extension mode
  // Set other parameter for subsidiary aligners
  subsidiary_attr.system = attributes->system;         // Inherit system configuration
  // The master BiWFA reports aggregate memory. Calling the same callback from
  // each subsidiary would double-count forward/reverse/base working sets.
  subsidiary_attr.system.memory_probe = NULL;
  subsidiary_attr.system.memory_probe_arguments = NULL;
  // Allocate forward/reverse aligners
  wf_bialigner->wf_forward = wavefront_aligner_new(&subsidiary_attr);
  wf_bialigner->wf_forward->align_mode = wf_align_biwfa_breakpoint_forward;
  wf_bialigner->wf_forward->plot = plot;
  wf_bialigner->wf_reverse = wavefront_aligner_new(&subsidiary_attr);
  wf_bialigner->wf_reverse->align_mode = wf_align_biwfa_breakpoint_reverse;
  wf_bialigner->wf_reverse->plot = plot;
  // Allocate subsidiary aligner
  subsidiary_attr.alignment_scope = compute_alignment;
  subsidiary_attr.heuristic.strategy = wf_heuristic_none; // Not inherited
  subsidiary_attr.system.memory_probe =
      &wavefront_bialigner_leaf_memory_probe;
  subsidiary_attr.system.memory_probe_arguments = wf_bialigner;
  wf_bialigner->wf_base = wavefront_aligner_new(&subsidiary_attr);
  wf_bialigner->wf_base->align_mode = wf_align_biwfa_subsidiary;
  wf_bialigner->wf_base->plot = plot;
  wf_bialigner->num_leaf_alignments = 0;
  wf_bialigner->max_leaf_score = 0;
  wf_bialigner->max_num_threads_used = 1;
  // Return
  return wf_bialigner;
}
void wavefront_bialigner_reap(
    wavefront_bialigner_t* const wf_bialigner) {
  wavefront_aligner_reap(wf_bialigner->wf_forward);
  wavefront_aligner_reap(wf_bialigner->wf_reverse);
  wavefront_aligner_reap(wf_bialigner->wf_base);
}
void wavefront_bialigner_delete(
    wavefront_bialigner_t* const wf_bialigner) {
  wavefront_aligner_delete(wf_bialigner->wf_forward);
  wavefront_aligner_delete(wf_bialigner->wf_reverse);
  wavefront_aligner_delete(wf_bialigner->wf_base);
  free(wf_bialigner);
}
/*
 * Sequences
 */
void wavefront_bialigner_set_sequences_ascii(
    wavefront_bialigner_t* const wf_bialigner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  wavefront_sequences_init_ascii(
      &wf_bialigner->wf_forward->sequences,
      pattern,pattern_length,text,text_length,false);
  wavefront_sequences_init_ascii(
      &wf_bialigner->wf_reverse->sequences,
      pattern,pattern_length,text,text_length,true);
  wavefront_sequences_init_ascii(
      &wf_bialigner->wf_base->sequences,
      pattern,pattern_length,text,text_length,false);
}
void wavefront_bialigner_set_sequences_lambda(
    wavefront_bialigner_t* const wf_bialigner,
    alignment_match_funct_t match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length) {
  wavefront_sequences_init_lambda(&wf_bialigner->wf_forward->sequences,
      match_funct,match_funct_arguments,pattern_length,text_length,false);
  wavefront_sequences_init_lambda(&wf_bialigner->wf_reverse->sequences,
      match_funct,match_funct_arguments,pattern_length,text_length,true);
  wavefront_sequences_init_lambda(&wf_bialigner->wf_base->sequences,
      match_funct,match_funct_arguments,pattern_length,text_length,false);
}
void wavefront_bialigner_set_sequences_packed2bits(
    wavefront_bialigner_t* const wf_bialigner,
    const uint8_t* const pattern,
    const int pattern_length,
    const uint8_t* const text,
    const int text_length) {
  wavefront_sequences_init_packed2bits(
      &wf_bialigner->wf_forward->sequences,
      pattern,pattern_length,text,text_length,false);
  wavefront_sequences_init_packed2bits(
      &wf_bialigner->wf_reverse->sequences,
      pattern,pattern_length,text,text_length,true);
  wavefront_sequences_init_packed2bits(
      &wf_bialigner->wf_base->sequences,
      pattern,pattern_length,text,text_length,false);
}
void wavefront_bialigner_set_sequences_bounds(
    wavefront_bialigner_t* const wf_bialigner,
    const int pattern_begin,
    const int pattern_end,
    const int text_begin,
    const int text_end) {
  wavefront_sequences_set_bounds(
      &wf_bialigner->wf_forward->sequences,
      pattern_begin,pattern_end,text_begin,text_end);
  wavefront_sequences_set_bounds(
      &wf_bialigner->wf_reverse->sequences,
      pattern_begin,pattern_end,text_begin,text_end);
  wavefront_sequences_set_bounds(
      &wf_bialigner->wf_base->sequences,
      pattern_begin,pattern_end,text_begin,text_end);
}
/*
 * Accessors
 */
uint64_t wavefront_bialigner_get_size(
    wavefront_bialigner_t* const wf_bialigner) {
  return wavefront_aligner_get_size(wf_bialigner->wf_forward) +
      wavefront_aligner_get_size(wf_bialigner->wf_reverse) +
      wavefront_aligner_get_size(wf_bialigner->wf_base);
}
void wavefront_bialigner_set_heuristic(
    wavefront_bialigner_t* const wf_bialigner,
    wavefront_heuristic_t* const heuristic) {
  wf_bialigner->wf_forward->heuristic = *heuristic;
  wf_bialigner->wf_reverse->heuristic = *heuristic;
  // Heuristics are not inherited to wf_base
}
void wavefront_bialigner_set_max_alignment_steps(
    wavefront_bialigner_t* const wf_bialigner,
    const int max_alignment_steps) {
  wf_bialigner->wf_forward->system.max_alignment_steps = max_alignment_steps;
  wf_bialigner->wf_reverse->system.max_alignment_steps = max_alignment_steps;
  wf_bialigner->wf_base->system.max_alignment_steps = max_alignment_steps;
}
void wavefront_bialigner_set_max_memory(
    wavefront_bialigner_t* const wf_bialigner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort) {
  wf_bialigner->wf_forward->system.max_memory_resident = max_memory_resident;
  wf_bialigner->wf_forward->system.max_memory_abort = max_memory_abort;
  wf_bialigner->wf_reverse->system.max_memory_resident = max_memory_resident;
  wf_bialigner->wf_reverse->system.max_memory_abort = max_memory_abort;
  wf_bialigner->wf_base->system.max_memory_resident = max_memory_resident;
  wf_bialigner->wf_base->system.max_memory_abort = max_memory_abort;
}
void wavefront_bialigner_set_memory_probe(
    wavefront_bialigner_t* const wf_bialigner,
    wavefront_memory_probe_funct_t memory_probe,
    void* memory_probe_arguments) {
  /* Aggregate probing belongs to the master BiWFA. Keep child callbacks off. */
  (void)memory_probe;
  (void)memory_probe_arguments;
  wf_bialigner->wf_forward->system.memory_probe = NULL;
  wf_bialigner->wf_forward->system.memory_probe_arguments = NULL;
  wf_bialigner->wf_reverse->system.memory_probe = NULL;
  wf_bialigner->wf_reverse->system.memory_probe_arguments = NULL;
  wf_bialigner->wf_base->system.memory_probe =
      &wavefront_bialigner_leaf_memory_probe;
  wf_bialigner->wf_base->system.memory_probe_arguments = wf_bialigner;
}
void wavefront_bialigner_set_leaf_score(
    wavefront_bialigner_t* const wf_bialigner,
    const int biwfa_leaf_score) {
  const int score = MAX(biwfa_leaf_score,0);
  wf_bialigner->wf_forward->system.biwfa_leaf_score = score;
  wf_bialigner->wf_reverse->system.biwfa_leaf_score = score;
  wf_bialigner->wf_base->system.biwfa_leaf_score = score;
}
void wavefront_bialigner_set_max_num_threads(
    wavefront_bialigner_t* const wf_bialigner,
    const int max_num_threads) {
  wf_bialigner->wf_forward->system.max_num_threads = max_num_threads;
  wf_bialigner->wf_reverse->system.max_num_threads = max_num_threads;
  wf_bialigner->wf_base->system.max_num_threads = max_num_threads;
}
void wavefront_bialigner_set_min_offsets_per_thread(
    wavefront_bialigner_t* const wf_bialigner,
    const int min_offsets_per_thread) {
  wf_bialigner->wf_forward->system.min_offsets_per_thread = min_offsets_per_thread;
  wf_bialigner->wf_reverse->system.min_offsets_per_thread = min_offsets_per_thread;
  wf_bialigner->wf_base->system.min_offsets_per_thread = min_offsets_per_thread;
}
