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

void gffToMaf(std::string & fastaFilePath,  std::string & gffFile, std::string & outfile, std::string newFastaFilePath ){
    std::map<std::string, std::string> referenceSequences;
    readFastaFileWorkWithIUPACcode( fastaFilePath, referenceSequences);

    std::map<std::string, std::string> teRemovedSequences;
    readFastaFileWorkWithIUPACcode( fastaFilePath, teRemovedSequences);

    std::ifstream infile(gffFile);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << gffFile << std::endl;
        exit (1);
    }

    std::regex reg0("^(\\S*)\t([\\s\\S]*)\tLTR_retrotransposon\t(\\S*)\t(\\S*)\t.*Parent=(\\w+);");
    std::regex reg1("^(\\S*)\t([\\s\\S]*)\ttarget_site_duplication\t(\\S*)\t(\\S*)\t.*Parent=(\\w+);");
    std::string line;
    std::string currentParent;
    while (std::getline(infile, line)){
        std::smatch match0;
        std::smatch match1;
        regex_search(line, match0, reg0);
        regex_search(line, match1, reg1);
        if( !match0.empty() && line[0]!='#' && line.size()>9){
            int start = stoi(match0[3]);
            int end = stoi(match0[4]);
            if(start>end){
                int temp=start;
                start = end;
                end = temp;
            }
            int length = end - start + 1;
            std::string newString = std::string(length, '-');
            std::string chr = match0[1];
            if( teRemovedSequences.find(chr) != teRemovedSequences.end() ){
                teRemovedSequences[chr].replace(start-1, length, newString);
            }
            currentParent = match0[5];
        } else if( !match1.empty() && line[0]!='#' && line.size()>9 && match1[5].compare(currentParent) != 0 ){
            int start = stoi(match1[3]);
            int end = stoi(match1[4]);
            if(start>end){
                int temp=start;
                start = end;
                end = temp;
            }
            int length = end - start + 1;
            std::string newString = std::string(length, '-');
            std::string chr = match1[1];
            if( teRemovedSequences.find(chr) != teRemovedSequences.end() ){
                teRemovedSequences[chr].replace(start-1, length, newString);
            }
        }
    }
    infile.close();

    std::ofstream omaffile;
    omaffile.open(outfile);
    omaffile << "##maf version=1" << std::endl;
    for( std::map<std::string, std::string>::iterator iterFasta=teRemovedSequences.begin(); iterFasta!=teRemovedSequences.end(); iterFasta++ ) {
        omaffile << "a\tscore=>" << 0 << std::endl;
        omaffile << "s\t" << std::left << std::setw(20) <<  "Zea_mays.AGPv4.dna.toplevel.fa." + iterFasta->first << "\t" << std::right
                 << std::setw(9) << 0 << "\t" << std::setw(9)
                 << referenceSequences[iterFasta->first].size() << "\t+\t" << referenceSequences[iterFasta->first].size() << "\t"
                 << referenceSequences[iterFasta->first] << std::endl;
        omaffile << "s\t" << std::left << std::setw(20) << "B73V4.pseudomolecule.subtract1.fa." + iterFasta->first << "\t" << std::right
                 << std::setw(9) << 0 << "\t" << std::setw(9)
                 << teRemovedSequences[iterFasta->first].size() << "\t+\t" << teRemovedSequences[iterFasta->first].size() << "\t"
                 << teRemovedSequences[iterFasta->first] << std::endl;
        omaffile << std::endl;
    }
    omaffile.close();

    std::vector<std::string> chrs;
    chrs.push_back("10");
    chrs.push_back("1");
    chrs.push_back("3");
    chrs.push_back("2");
    chrs.push_back("5");
    chrs.push_back("4");
    chrs.push_back("7");
    chrs.push_back("6");
    chrs.push_back("9");
    chrs.push_back("8");

    int linewidth = 70;
    omaffile.open(newFastaFilePath);
    for( std::string chr : chrs ) {
        omaffile << ">" << chr << std::endl;
        teRemovedSequences[chr].erase(std::remove(teRemovedSequences[chr].begin(), teRemovedSequences[chr].end(), '-'), teRemovedSequences[chr].end());
        std::size_t pos = 0;
        while (pos < teRemovedSequences[chr].length()) {
            omaffile << teRemovedSequences[chr].substr(pos, linewidth) << std::endl;
            pos += linewidth;
        }
    }
    omaffile.close();
}


TEST(gffToMaf, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/Zea_mays.AGPv4.dna.toplevel.fa";
    std::string newFastaFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/Zea_mays.AGPv4.dna.toplevel_withTe_removed.fa";
    std::string gffFile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/B73V4.pseudomolecule.ltrharvest.contignames.gff3.contigpositions.gff3";
    std::string outfile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/b73tebenckhmark.maf";
    gffToMaf(fastaFilePath, gffFile, outfile, newFastaFilePath );
    ASSERT_EQ(0, 0);
}



TEST(gff2Vcf, c1){ // just to make sure that every line has been analysed
    std::vector<std::string> chrs;
    chrs.push_back("10");
    chrs.push_back("1");
    chrs.push_back("3");
    chrs.push_back("2");
    chrs.push_back("5");
    chrs.push_back("4");
    chrs.push_back("7");
    chrs.push_back("6");
    chrs.push_back("9");
    chrs.push_back("8");
    std::string fastaFilePath = "/media/bs674/ppi8t/anchorwavepublishingresults/maizeTE/Zea_mays.AGPv4.dna.toplevel.fa";
    std::map<std::string, std::string> referenceGenome;
    readFastaFile( fastaFilePath, referenceGenome);
    std::cout << "genome reading done" << std::endl;

    std::ofstream ofile;
    ofile.open("/media/bs674/ppi8t/anchorwavepublishingresults/maizeTE/B73V4.pseudomolecule.ltrharvest.contignames.gff3.contigpositions.vcf");
    std::string gffFile = "/media/bs674/ppi8t/anchorwavepublishingresults/maizeTE/B73V4.pseudomolecule.ltrharvest.contignames.gff3.contigpositions.gff3";

    for( std::string chr : chrs ){
        std::vector<Variant> benchmarkVariants;
        gffToVariant(fastaFilePath,  gffFile, chr, benchmarkVariants );
        for(  Variant v : benchmarkVariants ){
            ofile << v.getChromosome() << "\t" << v.getPosition() << "\t" << v.getChanginglength() << "\t" << v.getReference() << "\t" << v.getAlternative() << std::endl;
        }
    }
    ofile.close();
    ASSERT_EQ(0, 0);
}




TEST(CompareWithVcf, c1){ // just to make sure that every line has been analysed
    std::vector<std::string> chrs;
    chrs.push_back("10");
    chrs.push_back("1");
    chrs.push_back("3");
    chrs.push_back("2");
    chrs.push_back("5");
    chrs.push_back("4");
    chrs.push_back("7");
    chrs.push_back("6");
    chrs.push_back("9");
    chrs.push_back("8");
    std::string fastaFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/Zea_mays.AGPv4.dna.toplevel.fa";
    std::map<std::string, std::string> referenceGenome;
    readFastaFile( fastaFilePath, referenceGenome);
    std::cout << "genome reading done" << std::endl;


    std::ofstream ofile2;
    ofile2.open("/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm23/gffTovariant.vcf");

    for( std::string chr : chrs ){
        std::ofstream ofile;

        ofile.open("/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm23/" + chr + "_subtract.vcf.errors");
        std::string vcfFilePath = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/lm23/subtract.vcf";


        std::map< std::string, std::vector<Variant>> variants;
        vcfToVariant(vcfFilePath, variants, chr );
        std::cout << "vcf to variant done" << std::endl;

        std::string gffFile = "/media/bs674/ppi8t/testWAF/maizeTeAlignment/B73V4.pseudomolecule.ltrharvest.contignames.gff3.contigpositions.gff3";
        std::vector<Variant> benchmarkVariants;
        gffToVariant(fastaFilePath,  gffFile, chr, benchmarkVariants );


        std::cout << "gff to variant done" << std::endl;
        for(  Variant v : benchmarkVariants ){
            ofile2 << v.getChromosome() << "\t" << v.getPosition() << "\t" << v.getChanginglength() << "\t" << v.getReference() << "\t" << v.getAlternative() << std::endl;
            if( variantWith(v, variants, referenceGenome) ){
                ofile << "good:" << v.getChromosome() << "\t" << v.getPosition() << "\t" << v.getChromosome() << "_" << v.getPosition() <<"\t"+ v.getReference() << "\t" << v.getAlternative() << std::endl;
            }else{
                ofile << "bad:" << v.getChromosome() << "\t" << v.getPosition() << "\t" << v.getChromosome() << "_" << v.getPosition() <<"\t"+ v.getReference() << "\t" << v.getAlternative() << std::endl;
            }
        }
        ofile.close();
    }
    ofile2.close();
    ASSERT_EQ(0, 0);
}

