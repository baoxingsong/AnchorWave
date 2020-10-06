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


int64_t alignSlidingWindow( const std::string& dna_q, const std::string& dna_d,
                            std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize, const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1 );

#endif //PROALI_ALIGNSLIDINGWINDOW_H
