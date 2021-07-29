//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include "../../service/service.h"
#include <string>
#include <iostream>
#include <map>
#include <regex>
#include <cstdlib>




TEST(sdiToMaf, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/tair10.fa";
    std::string targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/can_0.fa";
    std::string sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/can_0.v7c.sdi";
    std::string outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_can.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/bur_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/bur_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_bur.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/ct_1.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/ct_1.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_ct.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/edi_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/edi_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_edi.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/hi_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/hi_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_hi.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/kn_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/kn_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_kn.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/ler_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/ler_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_ler.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/mt_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/mt_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_mt.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/no_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/no_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_no.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/oy_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/oy_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_oy.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/po_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/po_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_po.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/rsch_4.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/rsch_4.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_rsch.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/sf_2.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/sf_2.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_sf.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/tsu_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/tsu_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_tsu.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/wil_2.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/wil_2.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_wil.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/ws_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/ws_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_ws.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/wu_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/wu_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_wu.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );

    targetFastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/zu_0.fa";
    sdiFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/zu_0.v7c.sdi";
    outfile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/col_zu.maf";
    sdiToMaf(fastaFilePath, targetFastaFilePath, sdiFile, outfile );
    ASSERT_EQ(0, 0);
}




TEST(samToMaf, c1){ // just to make sure that every line has been analysed
    std::string referenceFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/tair10.fa";
    std::string queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/bur_0.fa";
    std::string samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_bur_0.sam";
    std::string outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_bur_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/can_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_can_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_can_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/ct_1.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_ct_1.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_ct_1.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/edi_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_edi_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_edi_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/hi_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_hi_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_hi_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/kn_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_kn_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_kn_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/ler_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_ler_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_ler_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mt_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_mt_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_mt_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/no_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_no_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_no_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/oy_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_oy_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_oy_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/po_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_po_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_po_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/rsch_4.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_rsch_4.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_rsch_4.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/sf_2.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_sf_2.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_sf_2.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/tsu_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_tsu_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_tsu_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/wil_2.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_wil_2.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_wil_2.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/ws_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_ws_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_ws_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/wu_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_wu_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_wu_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/zu_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_zu_0.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/minimap2_zu_0.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    ASSERT_EQ(0, 0);
}



TEST(samToMaf, c2){ // just to make sure that every line has been analysed
    std::string referenceFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/tair10.fa";
    std::string queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/bur_0.fa";
    std::string samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.bur_0.short.sam";
    std::string outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.bur_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.bur_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.bur_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);


    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/can_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.can_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.can_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.can_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.can_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/ct_1.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ct_1.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ct_1.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ct_1.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ct_1.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/edi_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.edi_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.edi_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.edi_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.edi_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/hi_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.hi_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.hi_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.hi_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.hi_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/kn_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.kn_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.kn_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.kn_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.kn_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/ler_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ler_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ler_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ler_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ler_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mt_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.mt_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.mt_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.mt_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.mt_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/no_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.no_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.no_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.no_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.no_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/oy_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.oy_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.oy_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.oy_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.oy_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/po_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.po_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.po_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.po_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.po_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/rsch_4.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.rsch_4.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.rsch_4.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.rsch_4.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.rsch_4.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/sf_2.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.sf_2.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.sf_2.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.sf_2.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.sf_2.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/tsu_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.tsu_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.tsu_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.tsu_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.tsu_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/wil_2.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wil_2.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wil_2.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wil_2.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wil_2.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/ws_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ws_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ws_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ws_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.ws_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/wu_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wu_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wu_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wu_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.wu_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    queryFastaFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/zu_0.fa";
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.zu_0.short.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.zu_0.short.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);
    samFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.zu_0.long.sam";
    outputFilePath = "/media/bs674/ppi8t/testWAF/ara/colagainstSimulation/mumer.zu_0.long.maf";
    samToMaf( samFilePath, referenceFastaFilePath, queryFastaFilePath, outputFilePath);

    ASSERT_EQ(0, 0);
}
