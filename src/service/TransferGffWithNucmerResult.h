//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
#define ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H

#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include <map>
#include <regex>
#include <cstdlib>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "../../minimap2/minimap.h"
#include "../../minimap2/kseq.h"


void setupAnchorsWithSpliceAlignmentResult(const std::string &gffFilePath, const std::string &cdsSequenceFile, const std::string &samFile, std::map<std::string, std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                           double &inversion_PENALTY, double &MIN_ALIGNMENT_SCORE, bool &considerInversion, const int &minExon, const int64_t &windownWidth, const double &minimumSimilarity, const double &minimumSimilarity2,
                                           std::map<std::string, std::string> &parameters, std::map<std::string, std::string> &referenceGenome, std::map<std::string, std::string> &queryGenome,
                                           int &expectedCopies, double &maximumSimilarity, const std::string &referenceSamFilePath, const int32_t &wfaSize3, const bool &searchForNewAnchors, const bool &exonModel
        /*, const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
        const int32_t & extendGapPenalty1, const int & k, const int & w, const bool & H*/);

void setupAnchorsWithSpliceAlignmentResultQuota(const std::string &gffFilePath, const std::string &samFile, const std::string &cdsSequenceFile, std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                                double &INDEL_SCORE, double &GAP_OPEN_PENALTY, double &MIN_ALIGNMENT_SCORE, int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                                                double &calculateIndelDistance, const int &minExon, const int64_t &windownWidth, const double &minimumSimilarity, const double &minimumSimilarity2, std::map<std::string, std::string> &parameters,
                                                std::map<std::string, std::string> &referenceGenome,
                                                std::map<std::string, std::string> &queryGenome, int &expectedCopies, const int32_t &wfaSize3,
                                                double &maximumSimilarity, const std::string &referenceSamFilePath,
        /*const int32_t & matchingScore, const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int & k, const int & w, const bool & H,*/
                                                bool &searchForNewAnchors, const bool &exonModel);


#endif //ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
