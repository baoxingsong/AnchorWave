//
// Created by bs674 on 6/7/19.
//

// references:
// https://math.stackexchange.com/questions/129504/solving-a-sum-of-exponentials
// https://www.ugrad.math.ubc.ca/coursedoc/math100/notes/approx/newton.html
// http://www.bioinfo.rpi.edu/bystrc/courses/biol4540/lecture7.pdf // this is the only document I could find that gives clear equation to calculated lambda values
// use a newton method to calculate the lambda value

#include "calculateLambda.h"
size_t charFrequency(const int8_t * seq, const uint32_t & length, double * freq){
    size_t l = 0;
    int i;
    for ( i=0; i<length; ++i ){
        if( seq[i] == 0 ){
            freq[0] = freq[0] + 1.0;
            ++l;
        }else if ( seq[i] == 1 ){
            freq[1] = freq[1] + 1.0;
            ++l;
        }else if ( seq[i] == 3 ){
            freq[2] = freq[2] + 1.0;
            ++l;
        }else if ( seq[i] == 4 ){
            freq[3] = freq[3] + 1.0;
            ++l;
        }
    }
    for( i=0; i<4; ++i ){
        freq[i] = freq[i] / double(l);
    }
    return l;
}

double getF( const double * freq1, const double * freq2,  const int & matchingScore, const int & mismatchingPenalty, const double & lambda){
    double fValue = 0;
    int j;
    for( int i=0; i<4; ++i ){
        for( j=0; j<4; ++j ){
            if( i==j ){
                fValue += (freq1[i]*freq1[j]*exp(double(matchingScore)*lambda));
            }else{
                fValue += (freq1[i]*freq1[j]*exp(double(mismatchingPenalty)*lambda));
            }
        }
    }
    return fValue-1.0;
}

double getFp( const double * freq1, const double * freq2,  const int & matchingScore, const int & mismatchingPenalty, const double & lambda, const double & epsilon){
    double fPValue = 0;
    double lambda2 = lambda + epsilon;
    double f1 = getF( freq1, freq2,  matchingScore, mismatchingPenalty, lambda);
    double f2 = getF( freq1, freq2,  matchingScore, mismatchingPenalty, lambda2);
    fPValue = (f2-f1)/epsilon;
    return fPValue;
}

double getF2fp( const double * freq1, const double * freq2,  const int & matchingScore, const int & mismatchingPenalty, const double & lambda, const double & epsilon ){
    return getF( freq1, freq2, matchingScore, mismatchingPenalty, lambda)/getFp( freq1, freq2, matchingScore, mismatchingPenalty, lambda, epsilon);
}

double calculateLambda(const int8_t * seq1,  const int8_t * seq2, const uint32_t & length1, const uint32_t & length2, const int & matchingScore, const int & mismatchingPenalty, const double & epsilon, const double & precision){
    double freq1[4];
    double freq2[4];
    size_t l1 = charFrequency(seq1, length1, freq1);
    size_t l2 = charFrequency(seq2, length2, freq2);

    int i = 10000;
    double difference=10000;
    double lambda = 1.0;
    while( getFp( freq1, freq2, matchingScore, mismatchingPenalty, lambda, epsilon ) <0 ){
        lambda = lambda * 2.0;
    }
    while( i>0 &&  difference>precision){
        double lambda1 = lambda - getF2fp( freq1, freq2, matchingScore, mismatchingPenalty, lambda, epsilon );
        difference = abs(lambda1-lambda);
        lambda = lambda1;
    }
    return lambda;
}
