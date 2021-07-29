//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_GETSUBSEQUENCE_H
#define ANNOTATIONLIFTOVER_GETSUBSEQUENCE_H

#include "../model/model.h"
#include "GetReverseComplementary.h"

std::string getSubsequence( std::map<std::string, std::string>& sequences, const std::string& seqName, const int& start, const int& end);
std::string getSubsequence( std::map<std::string, std::string>& sequences, const std::string& seqName, const int& start, const int& end, const STRAND& strand);
std::string getSubsequence( const std::string & sequence, const int& _start, const int& _end);
std::string getSubsequence( const std::string & sequence, const int& _start, const int& _end, const STRAND& strand);
#endif //ANNOTATIONLIFTOVER_GETSUBSEQUENCE_H
