//
// Created by bs674 on 6/8/21.
//

#ifndef ANCHORWAVE_READMAFFORMSA_H
#define ANCHORWAVE_READMAFFORMSA_H

#include <string>
#include <iostream>
#include "../myImportandFunction/alignSlidingWindow.h"
#include <map>
#include "../model/model.h"
#include "../util/util.h"
#include "GetReverseComplementary.h"
#include "getSubsequence.h"

void readMafForMsa(const std::string &filePath, std::vector<BlocksForMsa> &blocksForMsas);

void ancestorInversion(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas);



//
// Created by bs674 on 6/8/21.
//

#include "readMafForMsa.h"

void readMafForMsa(const std::string &filePath, std::vector<BlocksForMsa> &blocksForMsas);

void ancestorInversion(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas);

void ancestorLink(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas);

void generateMsa(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas,
                 std::map<std::string, std::string> &reference1Genome, std::map<std::string, std::string> &reference2Genome);

void generateMsa(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas,
                 std::map<std::string, std::string> &reference1Genome, std::map<std::string, std::string> &reference2Genome,
                 int32_t &matchingScore, int32_t &mismatchingPenalty, int32_t &open_gap_penalty1, int32_t &extend_gap_penalty1,
                 int32_t &open_gap_penalty2, int32_t &extend_gap_penalty2, int64_t &slidingWindowSize, int32_t &wfaSize,
                 int32_t &min_wavefront_length, int32_t &max_distance_threshold, int &miniInsertionSize, int &maxDistance);

void outputAncestral(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas,
                     std::map<std::string, std::string> &reference1Genome, std::map<std::string, std::string> &reference2Genome, std::string &outPutFile);


#endif //ANCHORWAVE_READMAFFORMSA_H
