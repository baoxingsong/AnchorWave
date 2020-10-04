//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_READFASTAFILE_H
#define ANNOTATIONLIFTOVER_READFASTAFILE_H
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include "../model/model.h"
#include <sstream>

void readFastaFile( const std::string& filePath, std::map<std::string, Fasta>& sequences);


#endif //ANNOTATIONLIFTOVER_READFASTAFILE_H
