//
// Created by bs674 on 12/19/20.
//

#include "evaluation.h"



void vcfToVariant( std::string & vcfFilePath,  std::map< std::string, std::vector<Variant>> & variants , std::string & chrTotest){
    std::ifstream infile(vcfFilePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << vcfFilePath << std::endl;
        exit (1);
    }
    variants[chrTotest] = std::vector<Variant>();
    std::string line;
    char delim = '\t';
    while (std::getline(infile, line)){
        std::vector<std::string> elems;
        split(line, delim, elems);
        if( elems.size()>5 && line[0]!='#' && line.size()>9){
            std::string chr = elems[0];
            std::string result = elems[4];
            if( chr == chrTotest) { // this is very specific for this purpose
                int32_t start = stoi(elems[1]);
                std::string ori = elems[3];
                Variant mapSingleRecord = Variant(chr, start, ori, result);
//                std::cout << chr << "\r" << std::to_string(start) << "\t" << ori << "\t" << result << std::endl;
                variants[chr].push_back(mapSingleRecord);
            }
        }
    }
    infile.close();
}

bool equalVariant( Variant & v1, Variant & v2,  std::map<std::string, std::string>  & referenceGenome ){

    if( v1.getChromosome() == v2.getChromosome() && referenceGenome.find(v1.getChromosome())!=referenceGenome.end() && v1.getChanginglength() == v2.getChanginglength() ){
        if( v1.getPosition() == v2.getPosition() && v1.getReference() == v2.getReference()  ){
            return true;
        }
        std::string referenceSequence = referenceGenome[v1.getChromosome()];
        int totalSize = referenceSequence.size();
        std::string reserveString;
        reserveString.reserve(totalSize * 1.5);

        std::stringstream sequencestream1(reserveString);
        int currentPosition1 = 1;
        sequencestream1 << referenceSequence.substr( currentPosition1-1, v1.getPosition() - 1);
        sequencestream1 << v1.getAlternative();
        currentPosition1 = v1.getPosition() + v1.getReference().size();
        sequencestream1 << referenceSequence.substr(currentPosition1 - 1, totalSize - currentPosition1 + 1);


        std::stringstream sequencestream2(reserveString);
        int currentPosition2 = 1;
        sequencestream2 << referenceSequence.substr( currentPosition2-1, v2.getPosition() - 1);
        sequencestream2 << v2.getAlternative();
        currentPosition2 = v2.getPosition() + v2.getReference().size();
        sequencestream2 << referenceSequence.substr(currentPosition2 - 1, totalSize - currentPosition2 + 1);



        std::string seq1 = sequencestream1.str();
        std::string seq2 = sequencestream2.str();
        if( seq1 == seq2 ){
            return true;
        }
        return false;
    }
    return false;
}

bool variantWith( Variant & v, std::map< std::string, std::vector<Variant>> variants, std::map<std::string, std::string>  & referenceGenome ){
//
    if( variants.find(v.getChromosome()) != variants.end() ){
        for( Variant v2 : variants[v.getChromosome()] ){
            if( equalVariant(v, v2, referenceGenome) ){
                return true;
            }
        }
        return false;
    }
    std::cout << v.getChromosome() << std::endl;
    return false;
}



void gffToVariant(std::string & fastaFilePath,  std::string & gffFile, std::string & chrTotest, std::vector<Variant> & sdiRecordsThisOne ){

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
            std::string chr = match0[1];
            if(chr == chrTotest ){
                int length = end - start + 1;
                std::string newString = std::string(length, '-');
                if( teRemovedSequences.find(chr) != teRemovedSequences.end() ){
                    teRemovedSequences[chr].replace(start-1, length, newString);
                }
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
            std::string chr = match1[1];
            if(chr == chrTotest ) {
                int length = end - start + 1;
                std::string newString = std::string(length, '-');
                if (teRemovedSequences.find(chr) != teRemovedSequences.end()) {
                    teRemovedSequences[chr].replace(start - 1, length, newString);
                }
            }
        }
    }
    infile.close();

    std::string refAlignSeq = referenceSequences[chrTotest];
    alignmentToVcf(teRemovedSequences[chrTotest], refAlignSeq, sdiRecordsThisOne, chrTotest, referenceSequences, 0);
}
