//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_CHECKORFSTATE_H
#define ANNOTATIONLIFTOVER_CHECKORFSTATE_H

#include "../model/model.h"
#include "../util/util.h"
#include "getSubsequence.h"
#include "TranscriptUpdateInformation.h"
#include <atomic>
#include <thread>
#include <mutex>


void checkOrfState( Transcript& targetTranscript, std::map<std::string, Fasta>& targetGenome,
                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix, const int& minIntron);

bool ifSpliceSitesOk(Transcript& targetTranscript, Transcript &referenceTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     std::map<std::string, Fasta>& referenceGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool checkSpliceSites( std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);


bool checkSpliceSites( std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);




std::string agIUPACcodesTranslation(std::string& ag);
std::string gtIUPACcodesTranslation(std::string& gt);
bool ifLengthDivisibleByThree(std::string& sequence);
bool ifNewStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool ifEndWithStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool ifStartWithStartCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);




bool ifSpliceSitesOk(Transcript& targetTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);

bool checkSpliceSites(std::string& s1, std::string& s2, std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool ifIntronEnoughLarge(Transcript& targetTranscript, const int& minIntron);
void checkOrfPversionBin(int const  index, int const & binSize, std::vector<Transcript>   & transcripts, std::map<std::string, Fasta> & genome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,  int const & minIntron, std::atomic_int & number_of_runing_threads);
void checkOrfPversion(Transcript& transcript, std::map<std::string, Fasta>& genome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,  int& minIntron, std::atomic_int & number_of_runing_threads);

#endif //ANNOTATIONLIFTOVER_CHECKORFSTATE_H
