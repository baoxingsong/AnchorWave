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
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "InputParser.h"

InputParser::InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const std::string InputParser::getCmdOption( std::string &option) {
    std::vector<std::string>::iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    return "";
}

std::string InputParser::getCmdOption( const char* o) {
    std::string option = o;
    return getCmdOption(option);
}

bool InputParser::cmdOptionExists(std::string &option) {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}

bool InputParser::cmdOptionExists( const char* o){
    std::string option = o;
    return cmdOptionExists(option);
}

void usage( ){
    std::string progName = "proali";
    std::cout << "Program " << progName << std::endl <<
    "Usage: "<<progName<<" <command> [options]"<< std::endl <<
    "Commands:"<< std::endl <<
        "    gff2seq     get the protein/CDS/gene sequence from GFF/GTF file" << std::endl<<
        "    genoAli     global alignment and variant calling for whole chromosome" <<  std::endl <<
        "    proali      whole genome alignment with rearrangements" << std::endl;
}
