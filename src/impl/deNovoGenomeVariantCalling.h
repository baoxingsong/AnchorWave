//
// Created by Baoxing song on 20.10.18.
//

#ifndef PROALI_DENOVOGENOMEVARIANTCALLING_H
#define PROALI_DENOVOGENOMEVARIANTCALLING_H
#include <ctime>
#include "../model/model.h"
#include "../util/util.h"
#include "./readFastaFile.h"
#include "../myImportandFunction/myImportantFunction.h"
#include "./SequenceCharToUInt8.h"
#include <iomanip>

#include <atomic>
#include <mutex>
#include <unistd.h>
#include <thread>
#include <iostream>
#include <chrono>



void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                                      const int32_t & widownWidth, const int32_t & wfaSize, const int32_t & wfaSize2, const std::string & outPutMafFile, const std::string & outPutVcfFile,
                                      const std::string & outPutFragedFile, /*std::string & outPutLocalalignmentFile,*/ const int32_t & matchingScore,
                                      const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1,
                                      const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                                      const int32_t & min_wavefront_length, const int32_t & max_distance_threshold,
                                      int32_t & seed_window_size, const int32_t & mini_cns_score, const int32_t & step_size,
                                      const int32_t & matrix_boundary_distance, const int32_t & scoreThreshold, const  int32_t & w, const  int32_t & xDrop, const int & maxThread, std::map<std::string, std::string>& parameters );

void genomeAlignment( std::vector<std::vector<AlignmentMatch>> & alignmentMatchsMap,
                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                      const int32_t & widownWidth, const int32_t & wfaSize, const int32_t & wfaSize2,
                      const std::string & outPutMafFile, const std::string & outPutFragedFile, /*std::string & outPutLocalalignmentFile,*/
                      const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
                      const int32_t & extendGapPenalty1,
                      const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2, int32_t & seed_window_size, const int32_t & mini_cns_score, const int32_t & step_size,
                      const int32_t & matrix_boundary_distance, const  int32_t & scoreThreshold, const  int32_t & w, const  int32_t & xDrop,
                      const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const int & maxThread, std::map<std::string, std::string>& parameters );

void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::vector<Variant> & sdiRecordsThisOne, std::string chr, std::string & refSequence, int32_t refLetterNumber);
void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::map <std::string, std::string> & refSequences);

void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::map <std::string, std::string> & refSequences, int32_t refLetterNumber, const bool & gvcf);
void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::map <std::string, std::string> & refSequences, int32_t refLetterNumber);
void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::vector<Variant> & sdiRecordsThisOne, std::string chr, std::map <std::string, std::string> & refSequences, int32_t refLetterNumber);
void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::string & refSequence, int32_t refLetterNumber);
void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::string & refSequence, int32_t refLetterNumber, const bool & gvcf);

#endif //PROALI_DENOVOGENOMEVARIANTCALLING_H
