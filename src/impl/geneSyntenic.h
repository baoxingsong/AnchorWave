//
// Created by Baoxing Song on 2019-03-13.
//

#pragma once

#include "readFastaFile.h"
#include "TranscriptUpdateInformation.h"
#include "../model/AlignmentMatch.h"

#include <algorithm>
#include <assert.h>
#include <cmath>

void longestPath(std::vector<AlignmentMatch> &pairedSimilarFragments, std::vector<AlignmentMatch> &sortedOrthologPairs, const bool &keepTandemDuplication, double &scoreThreshold);

void myAlignmentMatchSort(std::vector<AlignmentMatch> &pairedSimilarFragments, const double &penalty, const double &scoreThreshold, const bool &keepTandemDuplication, const bool &considerInversion);

void myOrthologPairsSortQueryQuota(std::vector<AlignmentMatch> &pairedSimilarFragments);

void myOrthologPairsSortQuota(std::vector<AlignmentMatch> &pairedSimilarFragments);

void longestPathQuotav2(std::vector<AlignmentMatch> pairedSimilarFragments, std::vector<std::vector<AlignmentMatch>> &sortedOrthologPairChains,
                        std::map<std::string, std::map<int64_t, AlignmentMatch>> &refIndexMap, std::map<std::string, std::map<int64_t, AlignmentMatch>> &queryIndexMap,
                        double &INDEL_SCORE, double &GAP_OPEN_PENALTY,
                        double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                        double &calculateIndelDistance, bool withNovelAnchors);
