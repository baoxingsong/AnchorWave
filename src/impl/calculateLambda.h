//
// Created by bs674 on 6/7/19.
//

#ifndef AND_CNS_CALCULATELAMBDA_H
#define AND_CNS_CALCULATELAMBDA_H

#include <string>
#include <math.h>
#include <iostream>
double calculateLambda(const int8_t * seq1,  const int8_t * seq2, const uint32_t & length1, const uint32_t & length2, const int & matchingScore, const int & mismatchingPenalty, const double & epsilon, const double & precision);


#endif //AND_CNS_CALCULATELAMBDA_H
