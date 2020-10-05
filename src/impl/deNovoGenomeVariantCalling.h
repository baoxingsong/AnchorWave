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


void deNovoGenomeVariantCalling( std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                 const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                                 const size_t & widownWidth, const std::string & outPutMafFile, const std::string & outPutVcfFile, const std::string & outPutFragedFile);
void genomeAlignment( std::vector<std::vector<OrthologPair2>> & alignmentMatchsMap,
                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                      const size_t & widownWidth, const std::string & outPutFilePath, bool & outPutAlignmentForEachInterval);

#endif //PROALI_DENOVOGENOMEVARIANTCALLING_H
