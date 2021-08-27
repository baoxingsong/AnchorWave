//
// Created by song on 8/5/18.
//

#ifndef PROALI_ALIGNSLIDINGWINDOW_H
#define PROALI_ALIGNSLIDINGWINDOW_H

#include <string>
#include <stack>
#include <immintrin.h>
#include <map>
#include "../impl/impl.h"
#include "../util/nucleotideCodeSubstitutionMatrix.h"
#include "../util/parameters.h"
#include "../../WFA/gap_affine/affine_wavefront_align.h"

#include "../../WFA/utils/commons.h"
#include "../../WFA/system/profiler_timer.h"

#include "../../WFA/alignment/score_matrix.h"
#include "../../WFA/edit/edit_dp.h"
#include "../../WFA/gap_lineal/nw.h"
#include "../../WFA/gap_affine/affine_wavefront.h"
#include "../../WFA/gap_affine/swg.h"

#include "../../WFA/benchmark/benchmark_edit.h"
#include "../../WFA/benchmark/benchmark_gap_lineal.h"
#include "../../WFA/benchmark/benchmark_gap_affine.h"
#include "../../WFA/benchmark/benchmark_gap_affine2p.h"

#include "../../WFA/gap_affine2p/affine2p_penalties.h"
#include "../../WFA/gap_affine2p/affine2p_wavefront.h"
#include "../../WFA/gap_affine2p/affine2p_wavefront_align.h"

#include "../../WFA/gap_affine/affine_wavefront_align.h"
#include <stdlib.h>
#include "../../minimap2/ksw2.h"
//#include "../../minimap2/kalloc.h"

//int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d,
//                            std::string & _alignment_q, std::string & _alignment_d,
//                            affine_penalties_t* const affine_penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & wfaAdativeSize, const int32_t & matchingScore, const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
//                            const int32_t & extendGapPenalty1, const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m, bool tryAdapt );

int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m, std::map<std::string, std::string>& parameters );

int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m);

int64_t alignSlidingWindow( std::string& align_ref2, std::string& align_query, const std::string& dna_ref1,
                            std::string & align_ref1, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2);


//int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
//                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
//                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1,  const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
//                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m );
//
//int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
//                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize,  const int32_t & matchingScore,
//                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,  const int32_t & extendGapPenalty1, const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
//                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m );
//
//void inversionAlignment(  std::string& _dna_q,  std::string& _dna_d, const int32_t & matchingScore,
//                          const int32_t & mismatchingPenalty, const  int32_t & _open_gap_penalty, const int32_t & _extend_gap_penalty,
//                          const int32_t & _inversion_penalty, int32_t & maxDnaDIndex);
//
//void inversionAlignment(  std::string& _dna_q,  std::string& _dna_d, const int32_t & matchingScore,
//                          const int32_t & mismatchingPenalty, const  int32_t & _open_gap_penalty, const int32_t & _extend_gap_penalty,
//                          const int32_t & _inversion_penalty, int32_t & maxDnaDIndex, int32_t & maxDnaDIndexe);

int32_t needleAlignment(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                        const int32_t & mismatchingPenalty, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2);

int64_t alignSlidingWindow( const std::string& dna_q, const std::string& dna_d, int64_t _length_of_q, int64_t _length_of_d,
                            std::string & _alignment_q, std::string & _alignment_d, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2);

int64_t alignSlidingWindowNW(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                               affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                               const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1,  const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                               const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m );
//
//int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
//                             affine2p_penalties_t* const penalties/*, mm_allocator_t* const mm_allocatorddd*/, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
//                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
//                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m );

int64_t alignSlidingWindow_minimap2( const std::string& dna_q, const std::string& dna_d, int64_t _length_of_q, int64_t _length_of_d,
                                     std::string & _alignment_q, std::string & _alignment_d, const int64_t & slidingWindowSize,
                                     const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                                     const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2);
#endif //PROALI_ALIGNSLIDINGWINDOW_H
