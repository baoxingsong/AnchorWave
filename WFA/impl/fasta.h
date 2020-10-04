//
// Created by bs674 on 9/14/20.
//

#ifndef TESTWFA_FASTA_H
#define TESTWFA_FASTA_H


#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <vector>
#include <algorithm>
#include <regex>

void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences);

#endif //TESTWFA_FASTA_H
