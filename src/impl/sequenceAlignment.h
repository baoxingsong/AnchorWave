//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_SMITHWATERMAN_H
#define SONG_CNS_SMITHWATERMAN_H

#include <string>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include "../model/model.h"
#include <immintrin.h>
#include <stdio.h>
#include <string.h>


std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                   bool reverseAlignment, bool returnCigar);

// if you only want the maximum score and the position then use this one
void SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m, bool & positions);


void ssw_int8(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int8_t &_open_gap_penalty, const int8_t &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m, bool positions);

// this is a 2-piece affine gap cost smithwaterman method which return cigar
std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                   const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                   bool reverseAlignment, bool returnCigar);

void SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                   const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m);

// not well read yet, and not tested
std::vector<uint32_t> NeedlemanWunsch(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                   const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                   const Scorei & m, bool reverseAlignment, bool returnCigar);

std::vector<uint32_t> SemiGlobal(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                   const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                   bool returnCigar, const int32_t & zdrop, const int32_t & w, Matrix & T);

std::vector<uint32_t> mapCnsToGenome(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                     const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                     const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                     int32_t &maxScore, int32_t &endPosition1, const Scorei & m, Matrix & T);

void mapCnsToGenome(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                     const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                     const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                     int32_t &maxScore, int32_t &endPosition1, const Scorei & m);

void SemiGlobal_xextend(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                   const int32_t & xdrop, const int32_t & w);


//same with above one, but returen cigar
std::vector<uint32_t> SemiGlobal_xextend(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                         const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                                         int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & mi,
                                         const int32_t & xdrop, const int32_t & w, Matrix & T, int32_t & iii, int32_t  &jjj);

std::vector<uint32_t> SemiGlobal(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int16_t * weights, Score & score,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                   bool returnCigar, const int32_t & zdrop, const int32_t & w, Matrix & T);

std::vector<uint32_t> SemiGlobal_single_gap_penalty(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int16_t * weights, Score & score,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                   bool returnCigar, const int32_t & z, const int32_t & w, Matrix & T);
#endif //SONG_CNS_SMITHWATERMAN_H
