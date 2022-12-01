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
#include "../../WFA2-lib/bindings/cpp/WFAligner.hpp"
#include <stdlib.h>
#include "../../minimap2/ksw2.h"

int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                              const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m, std::map<std::string, std::string>& parameters );

int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                              const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m);

int64_t alignSlidingWindow( std::string& align_ref2, std::string& align_query, const std::string& dna_ref1,
                            std::string & align_ref1, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2);

int32_t needleAlignment(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                        const int32_t & mismatchingPenalty, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2);

int64_t alignSlidingWindow( const std::string& dna_q, const std::string& dna_d, int64_t _length_of_q, int64_t _length_of_d,
                            std::string & _alignment_q, std::string & _alignment_d, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2);

int64_t alignSlidingWindowNW(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                                const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                               const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1,  const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                               const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m );

int64_t alignSlidingWindow_minimap2( const std::string& dna_q, const std::string& dna_d, int64_t _length_of_q, int64_t _length_of_d,
                                     std::string & _alignment_q, std::string & _alignment_d, const int64_t & slidingWindowSize,
                                     const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                                     const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2);

int32_t needleAlignment(const std::vector<std::string> & _dna_refs, const std::vector<std::string> & _dna_queries, std::vector<std::stack<char>> & SQs, std::vector<std::stack<char>> & SRs,
                        const int32_t & mismatchingPenalty, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2);

#endif //PROALI_ALIGNSLIDINGWINDOW_H
