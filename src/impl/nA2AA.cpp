//
// Created by baoxing on 10/10/17.
//

#include "nA2AA.h"

std::string nA2AA( std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    std::string pattern="-";
    std::string pattern2="";
    seq=songStrReplaceAll(seq, pattern, pattern2);
    std::stringstream ssSb;
    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
    return ssSb.str();
}

void na2aaLocal(  std::string& seq, std::stringstream& ssSb, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    seq=songStrRemoveBlank(seq);
    transform(seq.begin(), seq.end(), seq.begin(),::toupper);
    std::string pattern="U";
    std::string pattern2="T";
    seq = songStrReplaceAll(seq, pattern, pattern2);

    for (size_t jj = 0; jj <= seq.length() - 3; jj += 3) {
        std::string codeSeq = seq.substr(jj, 3);
        BEGINMIDDLEEND position=MIDDLE;
        if( jj==0  ){
            position=BEGIN;
        }else if( jj==seq.length() - 3  ){
            position=END;
        }
        ssSb << nucleotideCodeSubstitutionMatrix.getGeneticCode( codeSeq, position );
    }
}
