//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../impl/readFastaFile.h"
TEST(readFastaFile, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/alignb73againstmo17/Zea_mays.AGPv4.dna.toplevel.fa";
    std::map<std::string, std::string> sequences;
//    readFastaFile( fastaFilePath, sequences);
    for( std::map<std::string, std::string>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
        std::cout << it->first << std::endl;
    }
    ASSERT_EQ(0, 0);
}
