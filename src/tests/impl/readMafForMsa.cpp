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

TEST(mafToSdi, c1){
    //std::string mafFile = "/media/bs674/ppi8t/anchorwavepublishingresults/arabidopsis/wu_0.f.maf";
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/mexicana_tripsacum.maf";
    mafFile = "/media/bs674/ppi8t/ZeaMSA/Zea_luxurians_maize.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/Tripsacum_dactyloides.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/Tripsacum_luxurians.maf";
    //mafFile = "/media/bs674/ppi8t/ZeaMSA/mexicana_tripsacum.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/Zm-B73-REFERENCE-NAM-5.0.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/Zea_luxurians_chromosomes-allmaps_v2.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
//    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/maize_luxurians_Ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}
//31m22s372ms



TEST(mafToSdi, c2){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/sorghum_coix/coix_sorghum.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/sorghum_coix/Andropogon_sorghum.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/sorghum_coix/Andropogon_coix.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/sorghum_coix/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/sorghum_coix/coix_lacryma.FASTA", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/sorghum_coix/AncestralSeqs.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}





TEST(mafToSdi, c3){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/Tripsacum_Zea_ancestral.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/sorghumCoix_ancestral_Zea_ancestral.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/sorghumCoix_ancestral.Tripsacum.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/Zea_ancestral.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/Tripsacum_dactyloides-southern_hifiasm-bionano_scaffolds_v1.0.fasta", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/Zea.Tripsacum_ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}





TEST(nostructureVariation, c1){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair10_bur_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}




TEST(ancestralNoInversion, c1){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair_bur_0_chr1_inversion.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair10_burIv_tair10.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}


TEST(ancestralTwoInversion, c1){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_bur_0_chr1_inversion_bur_0_chr1_inversion.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair10_burIv_burIv.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}


TEST(ref2differentStrand, c1){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/burrc_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_burrc.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_rc.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair10_burrc_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}





TEST(ref2querydifferentStrand, c1){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/burrc_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/burrc_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/burrc_burrc.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_rc.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair10_burrc_burrc.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}




TEST(ref2querydifferentStrandInversion, c1){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_rc_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur_0_chr1_inversion_rc.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_rc.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair10_burIrc_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}


TEST(ref2querydifferentStrandInversion, c2){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_rc_bur.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair_bur.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair_bur_0_chr1_inversion_rc.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1_inversion_rc.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_burIrc_tair10.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}



TEST(ref2queryInversion, c2){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2_tair.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_tair.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur_iv2.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/TAIR10_chr1.fas", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/tair_bur_iv2_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}




TEST(ref2queryInversion, c3){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2_bur.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur_iv2.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur_iv2_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}



TEST(ref2queryInversion, c4){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc_bur.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_buriv2rc.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_buriv2rc_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}


TEST(ref2queryInversion, c5){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_bur_iv2.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2_bur_iv2.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2_bur.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_iv2_bur_iv2_bur.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}


TEST(ref2queryInversion, c6){
    std::string mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/bur_buriv2rc.maf";
    std::vector<BlocksForMsa> ref1Ref2BlocksForMsas;
    readMafForMsa( mafFile, ref1Ref2BlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc_buriv2rc.maf";
    std::vector<BlocksForMsa> ref1QueryBlocksForMsas;
    readMafForMsa( mafFile, ref1QueryBlocksForMsas );

    mafFile = "/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc_bur.maf";
    std::vector<BlocksForMsa> ref2QueryBlocksForMsas;
    readMafForMsa( mafFile, ref2QueryBlocksForMsas );

    std::map<std::string, std::string> reference1Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc.fa", reference1Genome);
    std::map<std::string, std::string> reference2Genome;
    readFastaFile("/media/bs674/ppi8t/ZeaMSA/debug/bur_0_chr1.fa", reference2Genome);
    ancestorInversion( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas);
    ancestorLink( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas); //for this dataset, there is relocation, so no need to run this runction
    generateMsa( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome);
    std::string outPutFile = "/media/bs674/ppi8t/ZeaMSA/debug/buriv2rc_bur_buriv2rc.ancestral.fa";
    outputAncestral( ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome, reference2Genome, outPutFile);

    ASSERT_EQ(0, 0);
}
