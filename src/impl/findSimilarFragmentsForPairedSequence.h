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
                                                                           std::string & seq1_string, std::string & seq2_string, const int32_t & scoreThreshold,
                                                                           const int32_t & w, const int32_t & xDrop);


void link_seeds( std::set<Seed> & seeds );

//set up seeds using smither-waterman method
std::set<Seed> getSeeds( int8_t * seq1, int8_t * seq2, int8_t * seq1_rev_com, int8_t * seq2_rev_com,
                         const int32_t & length1, const int32_t & length2, const int32_t & miniReferenceSize, const int32_t & _open_gap_penalty1,
                         const int32_t & _extend_gap_penalty1,
                         const int32_t & startPosition2, const int32_t & mini_cns_score,
                         const int32_t & windowsSize, const Scorei & m, const int32_t & step_size );


//same with the above functions just give different parameters, and call the above function
void getSeeds( int8_t * seq1, int8_t * seq2, int8_t * seq1_rev_com, int8_t * seq2_rev_com, const int32_t & length1,
               const int32_t & length2, const int32_t & miniReferenceSize, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
               const int32_t & matrix_boundary_distance, const int32_t & startPosition2, const int32_t & mini_cns_score,
               const int32_t & windowsSize, int32_t & this_windowsSize_return, const Scorei & m,
               const int32_t & step_size, std::set<Seed> & seeds);

#endif //SONG_CNS_FINDSIMILARFRAGMENTSFORPAIREDSEQUENCE_H
