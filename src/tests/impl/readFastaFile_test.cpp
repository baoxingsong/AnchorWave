//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../impl/readFastaFile.h"
TEST(readFastaFile, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/netscratch/dep_tsiantis/grp_gan/song/geneStructureAnnotation/denovel/ler0/tair10.fa";
    std::map<std::string, Fasta> sequences;
    readFastaFile( fastaFilePath, sequences);
    for( std::map<std::string, Fasta>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
        std::cout << it->first << std::endl;
    }
    ASSERT_EQ(0, 0);
}
