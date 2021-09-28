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
    std::string fastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/alignb73againstmo17/Zea_mays.AGPv4.dna.toplevel.fa";
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    for( std::map<std::string, std::string>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
        std::cout << it->first << std::endl;
    }
    ASSERT_EQ(0, 0);
}


TEST(revc, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa";
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::ofstream ofile;
    ofile.open("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_rc.fa");
    for( std::map<std::string, std::string>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
        int lineWidth=60;
        writeFasta(ofile, it->first, getReverseComplementary(it->second),  lineWidth);
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}

TEST(revc, c2){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion.fa";
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::ofstream ofile;
    ofile.open("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_rc.fa");
    for( std::map<std::string, std::string>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
        int lineWidth=60;
        writeFasta(ofile, it->first, getReverseComplementary(it->second),  lineWidth);
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}

TEST(revc, c3){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2.fa";
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::ofstream ofile;
    ofile.open("/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc.fa");
    for( std::map<std::string, std::string>::iterator it=sequences.begin(); it!=sequences.end(); ++it ){
        int lineWidth=60;
        writeFasta(ofile, it->first, getReverseComplementary(it->second),  lineWidth);
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}
