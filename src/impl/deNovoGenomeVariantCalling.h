//
// Created by Baoxing song on 20.10.18.
//

#ifndef PROALI_DENOVOGENOMEVARIANTCALLING_H
#define PROALI_DENOVOGENOMEVARIANTCALLING_H
#include <ctime>
#include "../model/model.h"
#include "../util/util.h"
#include "./readFastaFile.h"
#include "./readGffFileWithEverything.h"
#include "../myImportandFunction/myImportantFunction.h"
#include "./SequenceCharToUInt8.h"
#include <iomanip>

void deNovoGenomeVariantCalling( std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                 const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                                 const size_t & widownWidth, const std::string & outPutMafFile, const std::string & outPutVcfFile, const std::string & outPutFragedFile);


void genomeAlignment( std::vector<std::vector<OrthologPair2>> & alignmentMatchsMap,
                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                      const size_t & widownWidth, const std::string & outPutFilePath, bool & outPutAlignmentForEachInterval,
                      const bool & localAlignment, const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
                      const int32_t & extendGapPenalty1, int32_t & seed_window_size, const  int32_t & mini_cns_score, const int32_t & step_size,
                      const int32_t & matrix_boundary_distance, const  int32_t & scoreThreshold, const  int32_t & w, const  int32_t & xDrop);



#endif //PROALI_DENOVOGENOMEVARIANTCALLING_H
