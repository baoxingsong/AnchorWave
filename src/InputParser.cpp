/*
 * =====================================================================================
 *
 *       Filename:  InputParser.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 01:13:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

#include "InputParser.h"

InputParser::InputParser(int &argc, char **argv) {
    for (int i = 1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const std::string InputParser::getCmdOption(std::string &option) {
    std::vector<std::string>::iterator itr;
    itr = std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
        return *itr;
    }

    return "";
}

std::string InputParser::getCmdOption(const char *o) {
    std::string option = o;
    return getCmdOption(option);
}

bool InputParser::cmdOptionExists(std::string &option) {
    return std::find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
}

bool InputParser::cmdOptionExists(const char *o) {
    std::string option = o;
    return cmdOptionExists(option);
}

void usage() {
    std::cout << "Program " << PROGRAMNAME << " v" << VERSION << std::endl <<
              "Usage: " << PROGRAMNAME << " <command> [options]" << std::endl <<
              "Commands:" << std::endl <<
              "    gff2seq     get the longest full-length CDS for each gene" << std::endl <<
              "    genoAli     whole chromosome global alignment and variant calling" << std::endl <<
              "    proali      genome alignment with relocation variation, chromosome fusion or whole genome duplication" << std::endl <<
              "    ali         perform global alignment for a pair of sequences using the 2-piece affine gap cost strategy" << std::endl
            ;
}
