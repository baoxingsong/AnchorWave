//
// Created by Baoxing Song on 2019-03-13.
//

#ifndef PROALI_LONGESTPATH_H
#define PROALI_LONGESTPATH_H

#include "../model/model.h"
#include "readGffFileWithEverything.h"
#include "readFastaFile.h"
#include "TranscriptUpdateInformation.h"
#include "checkOrfState.h"
#include <algorithm>

void longestPath (std::vector<AlignmentMatch> & pairedSimilarFragments, std::vector<AlignmentMatch> & sortedOrthologPairs, const bool & keepTandemDuplication);
void myAlignmentMatchSort( std::vector<AlignmentMatch> & pairedSimilarFragments, const double & score, const double & penalty, const double & scoreThreshold, const bool & keepTandemDuplication);


void myOrthologPairsSortQueryQuota( std::vector<OrthologPair2> & pairedSimilarFragments);
void myOrthologPairsSortQuota( std::vector<OrthologPair2> & pairedSimilarFragments);

void longestPathQuotav2 (std::vector<OrthologPair2> pairedSimilarFragments, std::vector<std::vector<OrthologPair2>> & sortedOrthologPairChains,
        double & INDEL_SCORE, double & GAP_OPEN_PENALTY,
                         double & MIN_ALIGNMENT_SCORE, const int & MAX_DIST_BETWEEN_MATCHES, int & refMaximumTimes, int & queryMaximumTimes,
                         double & calculateIndelDistance );

#endif //PROALI_LONGESTPATH_H
