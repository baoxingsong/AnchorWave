//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_TRANSCRIPTUPDATEINFORMATION_H
#define ANNOTATIONLIFTOVER_TRANSCRIPTUPDATEINFORMATION_H

#include "../model/model.h"
#include "getSubsequence.h"


void TranscriptUpdateCdsInformation(Transcript & transcript, std::map<std::string, Fasta>& genome);

#endif //ANNOTATIONLIFTOVER_TRANSCRIPTUPDATEINFORMATION_H
