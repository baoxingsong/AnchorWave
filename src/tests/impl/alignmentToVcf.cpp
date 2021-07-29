//
// Created by song on 9/19/18.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include "../../service/service.h"
#include "../../util/myutil.h"
#include <string>
#include <iostream>
#include <map>
#include <regex>
#include <cstdlib>
#include "../../controlLayer.h"

TEST(mafToSdi, c1){ // just to make sure that every line has been analysed
//    std::string outputovcffile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm23/outputvcf";
//    std::string fastaFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/Zea_mays.AGPv4.dna.toplevel.fa";
//    std::string mafFile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm23/subtract.maf";

//    std::string outputovcffile = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/output_bur0.vcf";
//    std::string fastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/tair10.fa";
//    std::string mafFile = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/tsu_0.maf";

//    mafTovcf(fastaFilePath, mafFile, fastaFilePath,  outputovcffile );

    std::string outputovcffile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm20/outputvcf";
    std::string fastaFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string mafFile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm20/subtract1.maf";
    mafTovcf( mafFile, fastaFilePath,  outputovcffile );

    fastaFilePath = "/media/bs674/ppi8t/testWAF/GSAlign_dataset/hg38.fa";
    outputovcffile = "/media/bs674/ppi8t/testWAF/GSAlign_dataset/v2.vcf";
    mafFile = "/media/bs674/ppi8t/testWAF/GSAlign_dataset/v2hg38-5.maf";
//    mafTovcf(fastaFilePath, mafFile, fastaFilePath,  outputovcffile );
    outputovcffile = "/media/bs674/ppi8t/testWAF/GSAlign_dataset/v3.vcf";
    mafFile = "/media/bs674/ppi8t/testWAF/GSAlign_dataset/v3hg38-5.maf";
 //   mafTovcf(fastaFilePath, mafFile, fastaFilePath,  outputovcffile );
    ASSERT_EQ(0, 0);
}
