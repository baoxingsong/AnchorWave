//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H
#define SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H

#include <algorithm>
#include "sequenceAlignment.h"
#include "../model/model.h"
#include <vector>
//#include "./fasta.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <math.h>
#include <iostream>
//#include "calculateLambda.h"
#include <limits>
#include <set>





std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                                                                           const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,const int32_t & matchingScore,
                                                                           const int32_t & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int32_t & w, const int32_t & xDrop);



// the second algorithm, continue from the output of first result
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a 2-piece gap cost approach
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                                                                           const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2, const int32_t & matchingScore,
                                                                           const int32_t & mismatchingPenalty, const Scorei & m,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                                                                           std::vector<Seed> & x_extend_seeds, class Matrix & T);

#endif //SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H
