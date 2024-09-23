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

void orthologPairSortMatchBin(std::vector<AlignmentMatch> &pairedSimilarFragments); // for delete tandem

void orthologPairSortQuery(std::vector<AlignmentMatch> &pairedSimilarFragments); // for delete tandem

void orthologPairSortReference(std::vector<AlignmentMatch> &pairedSimilarFragments); // for delete tandem

void orthologPairSortPosition(std::vector<AlignmentMatch> &pairedSimilarFragments);

void longestPathQuotaGeneSubAccelerate(std::vector<AlignmentMatch> pairedSimilarFragments,
                                                std::map<std::string, std::map<int, std::string>> &refIndexMap /*chr, index, refGeneName*/, std::map<std::string, std::map<int, std::string>> &queryIndexMap,
                                                double &GAP_EXTENSION_PENALTY, double &GAP_OPEN_PENALTY,
                                                double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES,
                                                int &refMaximumTimes, int &queryMaximumTimes,
                                                const int &count_style, const int &get_all_collinear_gene_pair,
                                                std::map<std::string, int64_t> &refTimes, std::map<std::string, int64_t> &queryTimes,
                                                std::set<std::string> &untouchedRefChrs, std::set<std::string> &untouchedQueryChrs,
                                                std::vector<double> &scoreArray, std::vector<int64_t> &prev);

void longestPathQuotaGeneNonStrandSubAccelerate(std::vector<AlignmentMatch> pairedSimilarFragments,
                                                std::map<std::string, std::map<int, std::string>> &refIndexMap /*chr, index, refGeneName*/, std::map<std::string, std::map<int, std::string>> &queryIndexMap,
                                                double &GAP_EXTENSION_PENALTY, double &GAP_OPEN_PENALTY,
                                                double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES,
                                                int &refMaximumTimes, int &queryMaximumTimes,
                                                const int &count_style, const int &get_all_collinear_gene_pair,
                                                std::map<std::string, int64_t> &refTimes, std::map<std::string, int64_t> &queryTimes,
                                                std::set<std::string> &untouchedRefChrs, std::set<std::string> &untouchedQueryChrs,
                                                std::vector<double> &scoreArray_positive, std::vector<double> &scoreArray_negative,
                                                std::vector<int64_t> &prev_positive, std::vector<int64_t> &prev_negative);

void longestPathQuotav2(std::vector<AlignmentMatch> pairedSimilarFragments, std::vector<std::vector<AlignmentMatch>> &sortedOrthologPairChains,
                        std::map<std::string, std::map<int64_t, AlignmentMatch>> &refIndexMap, std::map<std::string, std::map<int64_t, AlignmentMatch>> &queryIndexMap,
                        double &INDEL_SCORE, double &GAP_OPEN_PENALTY,
                        double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                        double &calculateIndelDistance, bool withNovelAnchors);

void longestPathQuotaGene(std::vector<AlignmentMatch> pairedSimilarFragments, std::vector<std::vector<AlignmentMatch>> &sortedOrthologPairChains,
                          std::map<std::string, std::map<int, std::string>> &refIndexMap /*chr, index, refGeneName*/, std::map<std::string, std::map<int, std::string>> &queryIndexMap,
                          double &GAP_EXTENSION_PENALTY, double &GAP_OPEN_PENALTY,
                          double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES,
                          int &refMaximumTimes, int &queryMaximumTimes,
                          std::vector<double> &block_score, const int &count_style, const int &get_all_collinear_gene_pair,
                          std::map<std::string, int64_t> &refTimes, std::map<std::string, int64_t> &queryTimes);

void longestPathQuotaGeneNonStrand(std::vector<AlignmentMatch> pairedSimilarFragments, std::vector<std::vector<AlignmentMatch>> &sortedOrthologPairChains,
                          std::map<std::string, std::map<int, std::string>> &refIndexMap /*chr, index, refGeneName*/, std::map<std::string, std::map<int, std::string>> &queryIndexMap,
                          double &GAP_EXTENSION_PENALTY, double &GAP_OPEN_PENALTY,
                          double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES,
                          int &refMaximumTimes, int &queryMaximumTimes,
                          std::vector<double> &block_score, const int &count_style, const int &get_all_collinear_gene_pair,
                          std::map<std::string, int64_t> &refTimes, std::map<std::string, int64_t> &queryTimes);
