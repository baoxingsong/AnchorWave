//
// Created by song on 8/19/18.
//

#ifndef ZSDP_CHECKANDUPDATETRANSCRIPTSENDS_H
#define ZSDP_CHECKANDUPDATETRANSCRIPTSENDS_H

#include "../model/model.h"
#include "TranscriptUpdateInformation.h"
#include "../util/util.h"
#include "checkOrfState.h"

void CheckAndUpdateTranscriptsEnds(std::map<std::string, Transcript> & Transcripts,
         std::map<std::string, Fasta>& sequences,  NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix, const int & minIntron);

void CheckAndUpdateTranscriptsEnds(std::map<std::string, std::vector<Transcript> > & transcriptHashSet, std::map<std::string, Fasta>& sequences,
                                   NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix, const int & minIntron);

#endif //ZSDP_CHECKANDUPDATETRANSCRIPTSENDS_H
