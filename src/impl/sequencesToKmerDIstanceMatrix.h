//
// Created by bs674 on 8/1/22.
//

#ifndef ANCHORWAVE_SEQUENCESTOKMERDISTANCEMATRIX_H
#define ANCHORWAVE_SEQUENCESTOKMERDISTANCEMATRIX_H

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include "../model/model.h"
#include <sstream>

void sequencesToKmerCountMatrix(std::map<std::string, std::string> &sequences, int8_t &kmer_size, std::map<std::string/*sequence name*/, std::map<std::string/*kmer*/, int16_t> > &kmerCountMatrix);

void kmerCountMatrixToDistanceMatrix(std::map<std::string/*sequence name*/, std::map<std::string/*kmer*/, int16_t> > &kmerCountMatrix, std::vector<std::string> &seqNames, float **distanceMatrix);


#endif //ANCHORWAVE_SEQUENCESTOKMERDISTANCEMATRIX_H
