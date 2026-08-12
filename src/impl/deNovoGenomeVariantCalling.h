//
// Created by Baoxing song on 20.10.18.
//

#pragma once

#include "getSubsequence.h"
#include "AlignmentBlockBuffer.h"
#include "readFastaFile.h"
#include "../model/AlignmentMatch.h"
#include "../myImportandFunction/alignSlidingWindow.h"
#include "../myImportandFunction/AlignmentAlgorithmSelector.h"
#include "../myImportandFunction/AlignmentResourceScheduler.h"
#include "../service/AlignmentGapPlanner.h"
#include "../service/AnchorTaskExecutor.h"
#include "../util/myutil.h"

#include <atomic>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <unistd.h>

void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> &map_v_am,
                                      const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                                      const int32_t &windowWidth, const std::string &outPutMafFile,
                                      const std::string &outPutFragedFile, const std::string &outPutBedFile, const int32_t &matchingScore,
                                      const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1,
                                      const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                      const int &maxThread,
                                      const uint64_t &maxProcessMemoryBytes);

void genomeAlignment(std::vector<std::vector<AlignmentMatch>> &v_v_am,
                     const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                     const int32_t &windowWidth,
                     const std::string &outPutMafFile, const std::string &outPutFragedFile, const std::string &outPutBedFile,
                     const int32_t &matchingScore, const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                     const int32_t &extendGapPenalty1,
                     const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                     const int &maxThread,
                     const uint64_t &maxProcessMemoryBytes);
