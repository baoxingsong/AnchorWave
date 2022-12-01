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
    std::string gffFilePath = "/home/bs674/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3";
    std::map<std::string, std::vector<std::string> > geneNameMap;
    std::map<std::string, Gene > geneHashMap;
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    int minExon = 20;
    std::string regex = "([\\s\\S]*)Parent=([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:_-]+)";
    readGffFile (gffFilePath, transcriptHashSet, regex, minExon);

    for( std::map<std::string, std::vector<Transcript> >::iterator it=transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it ){
        std::cout << it->first << " " << std::endl;
        for( Transcript transcriptId : it->second ){
            std::cout << transcriptId.getName() << " " << transcriptId.getPStart() << " " << transcriptId.getPEnd()
                      << " " << transcriptId.getCdsVector().size() << " CDSs " << transcriptId.getExonVector().size() << " exons" << std::endl;
        }
    }

    ASSERT_EQ(0, 0);
}
//
//TEST(readGffFileWithEveryThing, c2){// just to make sure that every line has been analysed
//    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes.gff";
//    std::map<std::string, std::vector<std::string> > geneNameMap;
//    std::map<std::string, Gene > geneHashMap;
//    std::map<std::string, Transcript> transcriptHashMap;
//    int minExon = 20;
//    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap, minExon);
//
//    for( std::map<std::string, std::vector<std::string> >::iterator it=geneNameMap.begin(); it != geneNameMap.end(); ++it ){
//        for( std::vector<std::string>::iterator it2=it->second.begin(); it2 != it->second.end(); ++it2 ){
//            std::cout << geneHashMap[*it2].getName() << " " << geneHashMap[*it2].getStart() << " " << geneHashMap[*it2].getEnd() << std::endl;
//            for( std::string transcriptId : geneHashMap[*it2].getTranscriptVector() ){
//                std::cout << transcriptId << " " << transcriptHashMap[transcriptId].getPStart() << " " << transcriptHashMap[transcriptId].getPEnd()
//                          << " " << transcriptHashMap[transcriptId].getCdsVector().size() << " CDSs " << transcriptHashMap[transcriptId].getExonVector().size() << " exons" << std::endl;
//            }
//        }
//    }
//    ASSERT_EQ(0, 0);
//}
//TEST(readGffFileWithEveryThing, c3){// just to make sure that every line has been analysed
//    std::string gffFilePath = "/biodata/dep_tsiantis/grp_gan/song/gennepredication/TAIR10annotationclassification/unmodifyed/TAIR10_GFF3_genes_transposons.gff";
//    std::map<std::string, std::vector<std::string> > geneNameMap;
//    std::map<std::string, Gene > geneHashMap;
//    std::map<std::string, Transcript> transcriptHashMap;
//    int minExon = 20;
//    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap, minExon);
//    ASSERT_EQ(0, 0);
//}
