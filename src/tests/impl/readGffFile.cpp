//
// Created by song on 8/8/18.
//

#include "include/gtest/gtest.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../service/TransferGffWithNucmerResult.h"
TEST(readGffFileWithEveryThing, c1){ // just to make sure that every line has been analysed
    std::string gffFilePath = "/home/bs674/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3";
    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene > geneHashMap;
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    int minExon = 20;
    std::string regex = "([\\s\\S]*)Parent=([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:_-]+)";
//    readGffFile (gffFilePath, transcriptHashSet, regex, minExon);

    for( std::map<std::string, std::vector<Transcript> >::iterator it=transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it ){
        std::cout << it->first << " " << std::endl;
        for( Transcript transcriptId : it->second ){
            std::cout << transcriptId.getName() << " " << transcriptId.getPStart() << " " << transcriptId.getPEnd()
                      << " " << transcriptId.getCdsVector().size() << " CDSs " << transcriptId.getExonVector().size() << " exons" << std::endl;
        }
    }

    ASSERT_EQ(0, 0);
}
