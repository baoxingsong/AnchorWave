//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_GETSEQUENCESFROMGFF_H
#define ANNOTATIONLIFTOVER_GETSEQUENCESFROMGFF_H

#include "../model/model.h"
#include "../util/util.h"
#include "readGffFile.h"
#include "readFastaFile.h"
#include "checkOrfState.h"
#include "TranscriptUpdateInformation.h"

void getSequences(const std::string& gffFile, const std::string& genome,
                  const std::string& outputCdsSequences,  std::map<std::string, std::string>& parameters, const int & minIntron);

#endif //ANNOTATIONLIFTOVER_GETSEQUENCESFROMGFF_H
