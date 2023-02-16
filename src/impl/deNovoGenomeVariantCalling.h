//
// Created by Baoxing song on 20.10.18.
//

#pragma once

#include "getSubsequence.h"
#include "readFastaFile.h"
#include "../model/AlignmentMatch.h"
#include "../myImportandFunction/alignSlidingWindow.h"
#include "../util/myutil.h"

#include <atomic>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <unistd.h>

void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> &map_v_am,
                                      const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                                      const int32_t &windowWidth, const std::string &outPutMafFile,
                                      const std::string &outPutFragedFile, const int32_t &matchingScore,
                                      const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1,
                                      const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                      const int &maxThread);

void genomeAlignment(std::vector<std::vector<AlignmentMatch>> &v_v_am,
                     const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                     const int32_t &windowWidth,
                     const std::string &outPutMafFile, const std::string &outPutFragedFile,
                     const int32_t &matchingScore, const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                     const int32_t &extendGapPenalty1,
                     const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                     const int &maxThread);
