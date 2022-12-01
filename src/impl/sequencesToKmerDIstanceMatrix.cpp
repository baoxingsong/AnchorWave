//
// Created by bs674 on 8/1/22.
//

#include "sequencesToKmerDIstanceMatrix.h"

void sequencesToKmerCountMatrix(std::map<std::string, std::string> &sequences, int8_t &kmer_size, std::map<std::string/*sequence name*/, std::map<std::string/*kmer*/, int16_t> > &kmerCountMatrix) {
    for (std::map<std::string, std::string>::iterator it = sequences.begin(); it != sequences.end(); ++it) {
        kmerCountMatrix[it->first] = std::map<std::string, int16_t>();
        for (int32_t position = 0; position <= it->second.length() - kmer_size; ++position) {
            std::string currentKmer = it->second.substr(position, kmer_size);
            if (kmerCountMatrix[it->first].find(currentKmer) != kmerCountMatrix[it->first].end()) {
                kmerCountMatrix[it->first][currentKmer] = kmerCountMatrix[it->first][currentKmer] + 1;
            } else {
                kmerCountMatrix[it->first][currentKmer] = 1;
            }
        }
    }
}

void kmerCountMatrixToDistanceMatrix(std::map<std::string/*sequence name*/, std::map<std::string/*kmer*/, int16_t> > &kmerCountMatrix, std::vector<std::string> &seqNames, float **distanceMatrix) {
    std::set<std::string> allKmers;
    for (std::map<std::string/*sequence name*/, std::map<std::string/*kmer*/, int16_t> >::iterator it = kmerCountMatrix.begin(); it != kmerCountMatrix.end(); ++it) {
        seqNames.push_back(it->first);
        for (std::map<std::string/*kmer*/, int16_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            allKmers.insert(it2->first);
        }
    }

    for (size_t i = 0; i < seqNames.size(); ++i) {
        distanceMatrix[i][i] = 0;
        for (size_t j = i + 1; j < seqNames.size(); ++j) {
            int32_t distance = 0;
            for (std::string kmer: allKmers) {
                int16_t counti = 0;
                int16_t countj = 0;
                if (kmerCountMatrix[seqNames[i]].find(kmer) != kmerCountMatrix[seqNames[i]].end()) {
                    counti = kmerCountMatrix[seqNames[i]][kmer];
                }
                if (kmerCountMatrix[seqNames[j]].find(kmer) != kmerCountMatrix[seqNames[j]].end()) {
                    countj = kmerCountMatrix[seqNames[j]][kmer];
                }
                distance += (counti - countj) * (counti - countj);
            }
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
    }
}