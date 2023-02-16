//
// Created by baoxing on 10/10/17.
//

#pragma once

#include "getSubsequence.h"
#include "../model/Transcript.h"

void TranscriptUpdateCdsInformation(Transcript &transcript, std::map<std::string, std::tuple<std::string, long, long, int> > &genome);
