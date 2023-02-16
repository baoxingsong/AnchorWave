//
// Created by baoxing on 10/10/17.
//

#pragma once

#include "readGffFile.h"
#include "readFastaFile.h"
#include "TranscriptUpdateInformation.h"

void getSequences(const std::string &gffFile, const std::string &genome, const std::string &outputCdsSequences, const int &minIntron, const bool &exonModel);
