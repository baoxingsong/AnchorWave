//
// Created by song on 8/8/18.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../service/TransferGffWithNucmerResult.h"
TEST(readGffFileWithEveryThing, c1){ // just to make sure that every line has been analysed
    std::string gffFilePath = "/netscratch/dep_tsiantis/grp_gan/song/zbl/haiwang/Zea_mays.AGPv4.38.gtf";
    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene > geneHashMap;
    std::map<std::string, Transcript> transcriptHashMap;
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);
    ASSERT_EQ(0, 0);
}
TEST(readGffFileWithEveryThing, c2){// just to make sure that every line has been analysed
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene > geneHashMap;
    std::map<std::string, Transcript> transcriptHashMap;
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);

    for( std::map<std::string, std::vector<std::string> >::iterator it=geneNameMap.begin(); it != geneNameMap.end(); ++it ){
        for( std::vector<std::string>::iterator it2=it->second.begin(); it2 != it->second.end(); ++it2 ){
            std::cout << geneHashMap[*it2].getName() << " " << geneHashMap[*it2].getStart() << " " << geneHashMap[*it2].getEnd() << std::endl;
            for( std::string transcriptId : geneHashMap[*it2].getTranscriptVector() ){
                std::cout << transcriptId << " " << transcriptHashMap[transcriptId].getPStart() << " " << transcriptHashMap[transcriptId].getPEnd()
                          << " " << transcriptHashMap[transcriptId].getCdsVector().size() << " CDSs " << transcriptHashMap[transcriptId].getExonVector().size() << " exons" << std::endl;
            }
        }
    }
    ASSERT_EQ(0, 0);
}
TEST(readGffFileWithEveryThing, c3){// just to make sure that every line has been analysed
    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes_transposons.gff";
    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene > geneHashMap;
    std::map<std::string, Transcript> transcriptHashMap;
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);
    ASSERT_EQ(0, 0);
}
