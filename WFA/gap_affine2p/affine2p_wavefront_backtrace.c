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
 * DESCRIPTION: WFA extend backtrace component
 */

#include "affine2p_wavefront_backtrace.h"
#include "affine2p_matrix.h"

/*
 * Wavefront type
 */
#define AFFINE2P_BACKTRACE_TYPE_BITS                   4 // 4-bits for piggyback
#define AFFINE2P_BACKTRACE_TYPE_MASK 0x000000000000000Fl // Extract mask

#define AFFINE2P_BACKTRACE_PIGGYBACK_SET(offset,backtrace_type) \
  (( ((int64_t)(offset)) << AFFINE2P_BACKTRACE_TYPE_BITS) | backtrace_type)

#define AFFINE2P_BACKTRACE_PIGGYBACK_GET_TYPE(offset) \
  ((offset) & AFFINE2P_BACKTRACE_TYPE_MASK)
#define AFFINE2P_BACKTRACE_PIGGYBACK_GET_OFFSET(offset) \
  ((offset) >> AFFINE2P_BACKTRACE_TYPE_BITS)

typedef enum {
  affine2p_backtrace_M       = 1,
  affine2p_backtrace_I1_open = 2,
  affine2p_backtrace_I1_ext  = 3,
  affine2p_backtrace_I2_open = 4,
  affine2p_backtrace_I2_ext  = 5,
  affine2p_backtrace_D1_open = 6,
  affine2p_backtrace_D1_ext  = 7,
  affine2p_backtrace_D2_open = 8,
  affine2p_backtrace_D2_ext  = 9
} affine2p_backtrace_type;

/*
 * Backtrace Detect Limits
 */
bool affine2p_wavefronts_valid_location(
    const int k,
    const awf2p_offset_t offset,
    const int pattern_length,
    const int text_length) {
  // Locate offset (remember that backtrace is always +1 offset ahead)
  const int v = AFFINE2P_WAVEFRONT_V(k,offset);
  const int h = AFFINE2P_WAVEFRONT_H(k,offset);
  return (v > 0 && v <= pattern_length &&
          h > 0 && h <= text_length);
}
void affine2p_wavefronts_add_trailing_gap(
    cigar_t* const cigar,
    const int k,
    const int alignment_k) {
  // Parameters
  char* const operations = cigar->operations;
  int op_sentinel = cigar->begin_offset;
  // Add trailing gap
  int i;
  if (k < alignment_k) {
    for (i=k;i<alignment_k;++i) operations[op_sentinel--] = 'I';
  } else if (k > alignment_k) {
    for (i=alignment_k;i<k;++i) operations[op_sentinel--] = 'D';
  }
  cigar->begin_offset = op_sentinel;
}
/*
 * Backtrace Trace Patch Match/Mismsmatch
 */
int64_t affine2p_wavefronts_trace_misms(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k &&
      k <= mwavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k]+1,affine2p_backtrace_M);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
void affine2p_wavefronts_trace_matches(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const char* const text,
    const int k,
    awf2p_offset_t offset,
    const bool valid_location,
    const int num_matches,
    cigar_t* const cigar) {
  int i;
  for (i=0;i<num_matches;++i) {
    // DEBUG
#ifdef AFFINE2P_WAVEFRONT_DEBUG
    const int v = AFFINE2P_WAVEFRONT_V(k,offset);
    const int h = AFFINE2P_WAVEFRONT_H(k,offset);
    if (!valid_location) { // Check inside matrix
      fprintf(stderr,"Backtrace error: Match outside DP-Matrix\n");
      exit(1);
    } else if (pattern[v-1] != text[h-1]) { // Check match
      fprintf(stderr,"Backtrace error: Not a match traceback\n");
      exit(1);
    }
#endif
    // Set Match
    cigar->operations[(cigar->begin_offset)--] = 'M';
    // Update state
    --offset;
  }
}
/*
 * Backtrace Trace Patch Deletion
 */
int64_t affine2p_wavefronts_trace_del1_open(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k+1 &&
      k+1 <= mwavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],affine2p_backtrace_D1_open);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
int64_t affine2p_wavefronts_trace_del2_open(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k+1 &&
      k+1 <= mwavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],affine2p_backtrace_D2_open);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
int64_t affine2p_wavefronts_trace_del1_ext(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const d1wavefront = affine2p_wavefronts->d1wavefronts[score];
  if (d1wavefront != NULL &&
      d1wavefront->lo_base <= k+1 &&
      k+1 <= d1wavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(d1wavefront->offsets[k+1],affine2p_backtrace_D1_ext);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
int64_t affine2p_wavefronts_trace_del2_ext(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const d2wavefront = affine2p_wavefronts->d2wavefronts[score];
  if (d2wavefront != NULL &&
      d2wavefront->lo_base <= k+1 &&
      k+1 <= d2wavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(d2wavefront->offsets[k+1],affine2p_backtrace_D2_ext);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace Trace Patch Insertion
 */
int64_t affine2p_wavefronts_trace_ins1_open(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k-1 &&
      k-1 <= mwavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,affine2p_backtrace_I1_open);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
int64_t affine2p_wavefronts_trace_ins2_open(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const mwavefront = affine2p_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k-1 &&
      k-1 <= mwavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,affine2p_backtrace_I2_open);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
int64_t affine2p_wavefronts_trace_ins1_ext(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const i1wavefront = affine2p_wavefronts->i1wavefronts[score];
  if (i1wavefront != NULL &&
      i1wavefront->lo_base <= k-1 &&
      k-1 <= i1wavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(i1wavefront->offsets[k-1]+1,affine2p_backtrace_I1_ext);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
int64_t affine2p_wavefronts_trace_ins2_ext(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const int score,
    const int k) {
  if (score < 0) return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  affine2p_wavefront_t* const i2wavefront = affine2p_wavefronts->i2wavefronts[score];
  if (i2wavefront != NULL &&
      i2wavefront->lo_base <= k-1 &&
      k-1 <= i2wavefront->hi_base) {
    return AFFINE2P_BACKTRACE_PIGGYBACK_SET(i2wavefront->offsets[k-1]+1,affine2p_backtrace_I2_ext);
  } else {
    return AFFINE2P_WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace
 */
void affine2p_wavefronts_backtrace(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int alignment_score) {
  // Parameters
  const affine2p_penalties_t* const wavefront_penalties = &(affine2p_wavefronts->penalties);
  cigar_t* const cigar = &affine2p_wavefronts->cigar;
  const int alignment_k = AFFINE2P_WAVEFRONT_DIAGONAL(text_length,pattern_length);
  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  awf2p_offset_t offset = affine2p_wavefronts->mwavefronts[alignment_score]->offsets[k];
  bool valid_location = affine2p_wavefronts_valid_location(k,offset,pattern_length,text_length);
  // Trace the alignment back
  affine2p_matrix_type matrix_type = affine2p_matrix_M;
  int v = AFFINE2P_WAVEFRONT_V(k,offset);
  int h = AFFINE2P_WAVEFRONT_H(k,offset);
  while (v > 0 && h > 0 && score > 0) {
    // Check location
    if (!valid_location) {
      valid_location = affine2p_wavefronts_valid_location(k,offset,pattern_length,text_length);
      if (valid_location) {
        affine2p_wavefronts_add_trailing_gap(cigar,k,alignment_k);
      }
    }
    // Compute scores
    const int mismatch = score - wavefront_penalties->mismatch;
    const int gap_open1 = score - wavefront_penalties->gap_opening1 - wavefront_penalties->gap_extension1;
    const int gap_open2 = score - wavefront_penalties->gap_opening2 - wavefront_penalties->gap_extension2;
    const int gap_extend1 = score - wavefront_penalties->gap_extension1;
    const int gap_extend2 = score - wavefront_penalties->gap_extension2;
    // Compute source offsets
    int64_t max_all;
    switch (matrix_type) {
      case affine2p_matrix_M: {
        const int64_t misms = affine2p_wavefronts_trace_misms(affine2p_wavefronts,mismatch,k);
        const int64_t ins1_open = affine2p_wavefronts_trace_ins1_open(affine2p_wavefronts,gap_open1,k);
        const int64_t ins1_ext  = affine2p_wavefronts_trace_ins1_ext(affine2p_wavefronts,gap_extend1,k);
        const int64_t ins2_open = affine2p_wavefronts_trace_ins2_open(affine2p_wavefronts,gap_open2,k);
        const int64_t ins2_ext  = affine2p_wavefronts_trace_ins2_ext(affine2p_wavefronts,gap_extend2,k);
        const int64_t del1_open = affine2p_wavefronts_trace_del1_open(affine2p_wavefronts,gap_open1,k);
        const int64_t del1_ext  = affine2p_wavefronts_trace_del1_ext(affine2p_wavefronts,gap_extend1,k);
        const int64_t del2_open = affine2p_wavefronts_trace_del2_open(affine2p_wavefronts,gap_open2,k);
        const int64_t del2_ext  = affine2p_wavefronts_trace_del2_ext(affine2p_wavefronts,gap_extend2,k);
        const int64_t max_ins_o = MAX(ins1_open,ins2_open);
        const int64_t max_ins_e = MAX(ins1_ext,ins2_ext);
        const int64_t max_del_o = MAX(del1_open,del2_open);
        const int64_t max_del_e = MAX(del1_ext,del2_ext);
        const int64_t max_ins = MAX(max_ins_o,max_ins_e);
        const int64_t max_del = MAX(max_del_o,max_del_e);
        max_all = MAX(misms,MAX(max_ins,max_del));
        break;
      }
      case affine2p_matrix_I1: {
        const int64_t ins1_open = affine2p_wavefronts_trace_ins1_open(affine2p_wavefronts,gap_open1,k);
        const int64_t ins1_ext  = affine2p_wavefronts_trace_ins1_ext(affine2p_wavefronts,gap_extend1,k);
        max_all = MAX(ins1_open,ins1_ext);
        break;
      }
      case affine2p_matrix_I2: {
        const int64_t ins2_open = affine2p_wavefronts_trace_ins2_open(affine2p_wavefronts,gap_open2,k);
        const int64_t ins2_ext  = affine2p_wavefronts_trace_ins2_ext(affine2p_wavefronts,gap_extend2,k);
        max_all = MAX(ins2_open,ins2_ext);
        break;
      }
      case affine2p_matrix_D1: {
        const int64_t del1_open = affine2p_wavefronts_trace_del1_open(affine2p_wavefronts,gap_open1,k);
        const int64_t del1_ext  = affine2p_wavefronts_trace_del1_ext(affine2p_wavefronts,gap_extend1,k);
        max_all = MAX(del1_open,del1_ext);
        break;
      }
      case affine2p_matrix_D2: {
        const int64_t del2_open = affine2p_wavefronts_trace_del2_open(affine2p_wavefronts,gap_open2,k);
        const int64_t del2_ext  = affine2p_wavefronts_trace_del2_ext(affine2p_wavefronts,gap_extend2,k);
        max_all = MAX(del2_open,del2_ext);
        break;
      }
      default:
        fprintf(stderr,"Affine-2P Wavefront-Backtrace. Wrong type trace (I)\n");
        exit(1);
        break;
    }
    // Traceback Matches
    if (matrix_type == affine2p_matrix_M) {
      const int max_offset = AFFINE2P_BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
      const int num_matches = offset - max_offset;
      affine2p_wavefronts_trace_matches(affine2p_wavefronts,
          pattern,text,k,offset,valid_location,num_matches,cigar);
      offset = max_offset;
      // Update coordinates
      v = AFFINE2P_WAVEFRONT_V(k,offset);
      h = AFFINE2P_WAVEFRONT_H(k,offset);
      if (v <= 0 || h <= 0) break;
    }
    // Traceback Operation
    const affine2p_backtrace_type backtrace_type = AFFINE2P_BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case affine2p_backtrace_M:
        score = mismatch;
        matrix_type = affine2p_matrix_M;
        break;
      case affine2p_backtrace_I1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case affine2p_backtrace_I1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_I1;
        break;
      case affine2p_backtrace_I2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case affine2p_backtrace_I2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_I2;
        break;
      case affine2p_backtrace_D1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case affine2p_backtrace_D1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_D1;
        break;
      case affine2p_backtrace_D2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case affine2p_backtrace_D2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_D2;
        break;
      default:
        fprintf(stderr,"Affine-2P Wavefront-Backtrace. Wrong type trace (II)\n");
        exit(1);
        break;
    }
    switch (backtrace_type) {
      case affine2p_backtrace_M:
        if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case affine2p_backtrace_I1_open:
      case affine2p_backtrace_I1_ext:
      case affine2p_backtrace_I2_open:
      case affine2p_backtrace_I2_ext:
        if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case affine2p_backtrace_D1_open:
      case affine2p_backtrace_D1_ext:
      case affine2p_backtrace_D2_open:
      case affine2p_backtrace_D2_ext:
        if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"Affine-2P Wavefront-Backtrace. Wrong type trace (II)\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = AFFINE2P_WAVEFRONT_V(k,offset);
    h = AFFINE2P_WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (score == 0) {
    // Account for last stroke of matches
    affine2p_wavefronts_trace_matches(affine2p_wavefronts,
        pattern,text,k,offset,valid_location,offset,cigar);
  } else {
    // Account for last stroke of insertion/deletion
    while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
    while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  }
  ++(cigar->begin_offset); // Set CIGAR length
}
