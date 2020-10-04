/*
 * =====================================================================================
 *
 *       Filename:  nucleotideCodeSubstitutionMatrix.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/03/2017 23:47:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include "nucleotideCodeSubstitutionMatrix.h"
#include "parameters.h"
#include "myutil.h"

#include <fstream>
#include <regex>
NucleotideCodeSubstitutionMatrix::NucleotideCodeSubstitutionMatrix(  ){

}
NucleotideCodeSubstitutionMatrix::NucleotideCodeSubstitutionMatrix( std::map<std::string, std::string>& parameters ){
    /*
    int8_t m[16][16] = {
        //  A,     T,     C,     G,     U,    N,      R,     Y,     S,      W,     K,     M,     B,     D,      H,      V
        {   36,   -36,   -36,   -36,   -36,    0,     18,   -36,   -36,     18,   -36,    18,   -36,    12,    12,    12} ,   //A
        {  -36,    36,   -36,   -36,    36,    0,    -36,    18,   -36,     18,    18,   -36,    12,    12,    12,   -36} ,   //T
        {  -36,   -36,    36,   -36,   -36,    0,    -36,    18,    18,    -36,   -36,    18,    12,   -36,    12,    12} ,   //C
        {  -36,   -36,   -36,    36,   -36,    0,     18,   -36,    18,    -36,    18,   -36,    12,    12,   -36,    12} ,   //G
        {  -36,    36,   -36,   -36,    36,    0,    -36,    18,   -36,     18,    18,   -36,    12,    12,    12,   -36} ,   //U
        {    0,     0,     0,     0,     0,    0,      0,     0,     0,      0,     0,     0,     0,     0,     0,     0} ,   //N
        {   18,   -36,   -36,    18,   -36,    0,     36,   -36,     9,      9,     9,     9,     6,    12,     6,    12} ,   //R
        {  -36,    18,    18,   -36,    18,    0,    -36,    36,     9,      9,     9,     9,    12,     6,    12,     6} ,   //Y
        {  -36,   -36,    18,    18,   -36,    0,      9,     9,    36,    -36,     9,     9,    12,     6,     6,    12} ,   //S
        {   18,    18,   -36,   -36,    18,    0,      9,     9,   -36,     36,     9,     9,     6,    12,    12,     6} ,   //W
        {  -36,    18,   -36,    18,    18,    0,      9,     9,     9,      9,    36,   -36,    12,    12,     6,     6} ,   //K
        {   18,   -36,    18,   -36,   -36,    0,      9,     9,     9,      9,   -36,    36,     6,     6,    12,    12} ,   //M
        {  -36,    12,    12,    12,    12,    0,      6,    12,    12,      6,    12,     6,    36,     8,     8,     8} ,   //B
        {   12,    12,   -36,    12,    12,    0,     12,     6,     6,     12,    12,     6,     8,    36,     8,     8} ,   //D
        {   12,    12,    12,   -36,    12,    0,      6,    12,     6,     12,     6,    12,     8,     8,    36,     8} ,   //H
        {   12,   -36,    12,    12,   -36,    0,     12,     6,    12,      6,     6,    12,     8,     8,     8,    36}     //V
    };
    for( size_t i=0; i<16; ++i ){
        for( size_t j=i; j<16; ++j ){
           nucleotide_substitution_matrix[i][j]=m[i][j];
           nucleotide_substitution_matrix[j][i]=m[j][i];
        }
    }

    _exon_subsitition_matrix = new int16_t*[16];
    _intron_subsitition_matrix = new int16_t*[16];
    _start_stop_codon_subsitition_matrix = new int16_t*[16];
    _splice_sites_subsitition_matrix = new int16_t*[16];
    for( int i=0; i<16; ++i ){
        _exon_subsitition_matrix[i] = new int16_t[16];
        _intron_subsitition_matrix[i] = new int16_t[16];
        _start_stop_codon_subsitition_matrix[i] = new int16_t[16];
        _splice_sites_subsitition_matrix[i] = new int16_t[16];
    }

    for( int i=0; i<16; i++){
        for( int j=0; j<16; j++){
            if( this->nucleotide_substitution_matrix[i][j] > 0 ){
                this->_exon_subsitition_matrix[i][j] = this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentExonMatchP", parameters));
                this->_intron_subsitition_matrix[i][j] = this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentIntronMatchP", parameters));
                this->_splice_sites_subsitition_matrix[i][j] = this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
                this->_start_stop_codon_subsitition_matrix[i][j] = this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
            }else{
                this->_exon_subsitition_matrix[i][j] = -(this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentExonMismatchP", parameters)));
                this->_intron_subsitition_matrix[i][j] = -(this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentIntronMismatchP", parameters)));
                this->_splice_sites_subsitition_matrix[i][j] = -(this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters)));
                this->_start_stop_codon_subsitition_matrix[i][j] = -(this->nucleotide_substitution_matrix[i][j] * stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters)));
            }
        }
    }
    dna_acid_map['A'] = 0;
    dna_acid_map['T'] = 1;
    dna_acid_map['C'] = 2;
    dna_acid_map['G'] = 3;
    dna_acid_map['U'] = 4;
    dna_acid_map['N'] = 5;
    dna_acid_map['R'] = 6;
    dna_acid_map['Y'] = 7;
    dna_acid_map['S'] = 8;
    dna_acid_map['W'] = 9;
    dna_acid_map['K'] = 10;
    dna_acid_map['M'] = 11;
    dna_acid_map['B'] = 12;
    dna_acid_map['D'] = 13;
    dna_acid_map['H'] = 14;
    dna_acid_map['V'] = 15;
     */

    int8_t m[5][5] = {
            //  A,     T,     C,     G,     N
            {    2,    -4,    -4,    -4,    0} ,   //A
            {   -4,     2,    -4,    -4,    0} ,   //T
            {   -4,    -4,     2,    -4,    0} ,   //C
            {   -4,    -4,    -4,     2,    0},     //G
            {    0,     0,     0,     0,    0}     //N
    };
    int i, j;
    for( i=0; i<5;++i ){
        for( j=0; j<5;++j ){
            nucleotide_substitution_matrix[i][j] = m[i][j];
        }
    }

    _exon_subsitition_matrix = new int8_t*[5];
    _intron_subsitition_matrix = new int8_t*[5];
    _start_stop_codon_subsitition_matrix = new int8_t*[5];
    _splice_sites_subsitition_matrix = new int8_t*[5];
    for( i=0; i<5; ++i ){
        _exon_subsitition_matrix[i] = new int8_t[5];
        _intron_subsitition_matrix[i] = new int8_t[5];
        _start_stop_codon_subsitition_matrix[i] = new int8_t[5];
        _splice_sites_subsitition_matrix[i] = new int8_t[5];
    }

    for( i=0; i<5; i++){
        for( j=0; j<5; j++){
            if( i==4 || j==4 ){
                this->_exon_subsitition_matrix[i][j] = 0;
            }else if( i==j ){
                this->_exon_subsitition_matrix[i][j] = stoi(get_parameters("alignmentExonMatchP", parameters));
                this->_intron_subsitition_matrix[i][j] = stoi(get_parameters("alignmentIntronMatchP", parameters));
                this->_splice_sites_subsitition_matrix[i][j] = stoi(get_parameters("alignmentSpliceSitesMatchP", parameters));
                this->_start_stop_codon_subsitition_matrix[i][j] = stoi(get_parameters("alignmentStartStopCodonMatchP", parameters));
            }else{
                this->_exon_subsitition_matrix[i][j] = stoi(get_parameters("alignmentExonMismatchP", parameters));
                this->_intron_subsitition_matrix[i][j] = stoi(get_parameters("alignmentIntronMismatchP", parameters));
                this->_splice_sites_subsitition_matrix[i][j] = stoi(get_parameters("alignmentSpliceSitesMismatchP", parameters));
                this->_start_stop_codon_subsitition_matrix[i][j] = stoi(get_parameters("alignmentStartStopCodonMismatchP", parameters));
            }
        }
    }

    dna_acid_map['A']=0;
    dna_acid_map['a']=0;
    dna_acid_map['T']=1;
    dna_acid_map['t']=1;
    dna_acid_map['U']=1;
    dna_acid_map['u']=1;
    dna_acid_map['C']=2;
    dna_acid_map['c']=2;
    dna_acid_map['G']=3;
    dna_acid_map['g']=3;

//
//
//    allStartCodons["TTG"]=0.6;
//    allStartCodons["ATG"]=1;
//    allStartCodons["CTG"]=0.6;
//    allStartCodons["YTG"]=0.6;
//    allStartCodons["WTG"]=0.8;
//    allStartCodons["MTG"]=0.8;
//    allStartCodons["HTG"]=0.8;
//    allStartCodons["RTG"]=0.5;
//    allStartCodons["DTG"]=0.5;
//    allStartCodons["VTG"]=0.5;
//    allStartCodons["NTG"]=0.5;
//    allStartCodons["KTG"]=0.4;
//    allStartCodons["KTG"]=0.3;
    
//    allStopCodons["TAA"]=1;
//    allStopCodons["TAG"]=1;
//    allStopCodons["TGA"]=1;
//    allStopCodons["TRA"]=1;
//    allStopCodons["TAR"]=1;

    mustStartCodons.insert("TTG");
    mustStartCodons.insert("CTG");
    mustStartCodons.insert("ATG");
    mustStartCodons.insert("YTG");
    mustStartCodons.insert("WTG");
    mustStartCodons.insert("MTG");
    mustStartCodons.insert("HTG");

    mustStopCodons.insert("TAA");
    mustStopCodons.insert("TAG");
    mustStopCodons.insert("TGA");
    mustStopCodons.insert("TAR");
    mustStopCodons.insert("TRA");

    basicStartCodons.insert("ATG");
    basicStartCodons.insert("TTG");
    basicStartCodons.insert("CTG");

    basicStopCodons.insert("TAA");
    basicStopCodons.insert("TAG");
    basicStopCodons.insert("TGA");


    dnaIupacCode["A"]=std::set<std::string>();
    dnaIupacCode["A"].insert("A");

    dnaIupacCode["C"]=std::set<std::string>();
    dnaIupacCode["C"].insert("C");

    dnaIupacCode["G"]=std::set<std::string>();
    dnaIupacCode["G"].insert("G");

    dnaIupacCode["T"]=std::set<std::string>();
    dnaIupacCode["T"].insert("T");
    dnaIupacCode["T"].insert("U");

    dnaIupacCode["U"]=std::set<std::string>();
    dnaIupacCode["U"].insert("U");
    dnaIupacCode["U"].insert("T");
    dnaIupacCode["U"].insert("U");

    dnaIupacCode["R"]=std::set<std::string>();
    dnaIupacCode["R"].insert("A");
    dnaIupacCode["R"].insert("G");

    dnaIupacCode["Y"]=std::set<std::string>();
    dnaIupacCode["Y"].insert("C");
    dnaIupacCode["Y"].insert("T");
    dnaIupacCode["Y"].insert("U");

    dnaIupacCode["S"]=std::set<std::string>();
    dnaIupacCode["S"].insert("G");
    dnaIupacCode["S"].insert("C");

    dnaIupacCode["W"]=std::set<std::string>();
    dnaIupacCode["W"].insert("A");
    dnaIupacCode["W"].insert("T");
    dnaIupacCode["W"].insert("U");

    dnaIupacCode["K"]=std::set<std::string>();
    dnaIupacCode["K"].insert("G");
    dnaIupacCode["K"].insert("T");
    dnaIupacCode["K"].insert("U");

    dnaIupacCode["M"]=std::set<std::string>();
    dnaIupacCode["M"].insert("A");
    dnaIupacCode["M"].insert("C");

    dnaIupacCode["B"]=std::set<std::string>();
    dnaIupacCode["B"].insert("C");
    dnaIupacCode["B"].insert("T");
    dnaIupacCode["B"].insert("U");
    dnaIupacCode["B"].insert("G");

    dnaIupacCode["D"]=std::set<std::string>();
    dnaIupacCode["D"].insert("A");
    dnaIupacCode["D"].insert("T");
    dnaIupacCode["D"].insert("U");
    dnaIupacCode["D"].insert("G");

    dnaIupacCode["H"]=std::set<std::string>();
    dnaIupacCode["H"].insert("A");
    dnaIupacCode["H"].insert("T");
    dnaIupacCode["H"].insert("U");
    dnaIupacCode["H"].insert("C");

    dnaIupacCode["V"]=std::set<std::string>();
    dnaIupacCode["V"].insert("A");
    dnaIupacCode["V"].insert("G");
    dnaIupacCode["V"].insert("C");

    dnaIupacCode["N"]=std::set<std::string>();
    dnaIupacCode["N"].insert("A");
    dnaIupacCode["N"].insert("G");
    dnaIupacCode["N"].insert("C");
    dnaIupacCode["N"].insert("T");
    dnaIupacCode["N"].insert("U");


    for(std::map<std::string, std::set<std::string> >::iterator itIupac=dnaIupacCode.begin(); itIupac!=dnaIupacCode.end(); ++itIupac ){
        std::string iupac1 = itIupac->first;
        for( std::set<std::string>::iterator itIupac1 = itIupac->second.begin(); itIupac1!=itIupac->second.end(); ++itIupac1 ){
            std::string iupac2 = (*itIupac1);
            if( revDnaIupacCode.find(iupac2) == revDnaIupacCode.end() ){
                revDnaIupacCode[iupac2] = std::set<std::string>();
            }
            revDnaIupacCode[iupac2].insert(iupac1);
        }
    }


    for( std::set<std::string>::iterator it1=basicStartCodons.begin(); it1!=basicStartCodons.end(); ++it1){
        std::string a=(*it1).substr(0,1);
        std::string b=(*it1).substr(1,1);
        std::string c=(*it1).substr(2,1);
        //std::cout << "198 " << a << b << c << std::endl;
        for( std::map<std::string, std::set<std::string> >::iterator it2=dnaIupacCode.begin(); it2!=dnaIupacCode.end(); ++it2){
            if( it2->second.find(a) != it2->second.end() ){
                for( std::map<std::string, std::set<std::string> >::iterator it3=dnaIupacCode.begin(); it3!=dnaIupacCode.end(); ++it3){
                    if( it3->second.find(b) != it3->second.end() ){
                        for( std::map<std::string, std::set<std::string> >::iterator it4=dnaIupacCode.begin(); it4!=dnaIupacCode.end(); ++it4){
                            if( it4->second.find(c) != it4->second.end() ){
                                std::string threeNa = it2->first+it3->first+it4->first;
                                //std::cout << "206 " << threeNa << std::endl;
                                this->possibleStartCodons.insert(threeNa);
                            }
                        }
                    }
                }
            }
        }
    }

    for( std::set<std::string>::iterator it1=basicStopCodons.begin(); it1!=basicStopCodons.end(); ++it1){
        std::string a=(*it1).substr(0,1);
        std::string b=(*it1).substr(1,1);
        std::string c=(*it1).substr(2,1);
        //std::cout << "220 " << a << b << c << std::endl;
        for( std::map<std::string, std::set<std::string> >::iterator it2=dnaIupacCode.begin(); it2!=dnaIupacCode.end(); ++it2){
            if( it2->second.find(a) != it2->second.end() ){
                for( std::map<std::string, std::set<std::string> >::iterator it3=dnaIupacCode.begin(); it3!=dnaIupacCode.end(); ++it3){
                    if( it3->second.find(b) != it3->second.end() ){
                        for( std::map<std::string, std::set<std::string> >::iterator it4=dnaIupacCode.begin(); it4!=dnaIupacCode.end(); ++it4){
                            if( it4->second.find(c) != it4->second.end() ){
                                std::string threeNa = it2->first+it3->first+it4->first;
                                //std::cout << "227 " << threeNa << std::endl;
                                this->possibleStopCodons.insert(threeNa);
                            }
                        }
                    }
                }
            }
        }
    }


    ////////////////////*
    middleStandardGeneticCode["TAA"]='*';
    middleStandardGeneticCode["TAG"]='*';
    middleStandardGeneticCode["TGA"]='*';
    /////
    middleStandardGeneticCode["TRA"]='*';
    middleStandardGeneticCode["TAR"]='*';

    ////////////////////F
    middleStandardGeneticCode["TTT"]='F';
    middleStandardGeneticCode["TTC"]='F';
    /////
    middleStandardGeneticCode["TTY"]='F';

    ////////////////////L
    middleStandardGeneticCode["TTA"]='L';
    middleStandardGeneticCode["TTG"]='L';

    middleStandardGeneticCode["CTC"]='L';
    middleStandardGeneticCode["CTT"]='L';
    middleStandardGeneticCode["CTA"]='L';
    middleStandardGeneticCode["CTG"]='L';
    /////
    middleStandardGeneticCode["TTR"]='L';

    middleStandardGeneticCode["YTA"]='L';
    middleStandardGeneticCode["YTG"]='L';
    middleStandardGeneticCode["YTR"]='L';

    middleStandardGeneticCode["CTR"]='L';
    middleStandardGeneticCode["CTY"]='L';
    middleStandardGeneticCode["CTS"]='L';
    middleStandardGeneticCode["CTW"]='L';
    middleStandardGeneticCode["CTK"]='L';
    middleStandardGeneticCode["CTM"]='L';
    middleStandardGeneticCode["CTB"]='L';
    middleStandardGeneticCode["CTD"]='L';
    middleStandardGeneticCode["CTH"]='L';
    middleStandardGeneticCode["CTV"]='L';
    middleStandardGeneticCode["CTN"]='L';

    ////////////////////S
    middleStandardGeneticCode["TCC"]='S';
    middleStandardGeneticCode["TCT"]='S';
    middleStandardGeneticCode["TCA"]='S';
    middleStandardGeneticCode["TCG"]='S';

    middleStandardGeneticCode["TCR"]='S';
    middleStandardGeneticCode["TCY"]='S';
    middleStandardGeneticCode["TCS"]='S';
    middleStandardGeneticCode["TCW"]='S';
    middleStandardGeneticCode["TCK"]='S';
    middleStandardGeneticCode["TCM"]='S';
    middleStandardGeneticCode["TCB"]='S';
    middleStandardGeneticCode["TCD"]='S';
    middleStandardGeneticCode["TCH"]='S';
    middleStandardGeneticCode["TCV"]='S';
    middleStandardGeneticCode["TCN"]='S';

    middleStandardGeneticCode["AGT"]='S';
    middleStandardGeneticCode["AGC"]='S';
    middleStandardGeneticCode["AGY"]='S';

    ////////////////////Y
    middleStandardGeneticCode["TAT"]='Y';
    middleStandardGeneticCode["TAC"]='Y';
    middleStandardGeneticCode["TAY"]='Y';

    ////////////////////C
    middleStandardGeneticCode["TGT"]='C';
    middleStandardGeneticCode["TGC"]='C';
    middleStandardGeneticCode["TGY"]='C';

    ///////////////////W
    middleStandardGeneticCode["TGG"]='W';

    ///////////////////P
    middleStandardGeneticCode["CCT"]='P';
    middleStandardGeneticCode["CCC"]='P';
    middleStandardGeneticCode["CCA"]='P';
    middleStandardGeneticCode["CCG"]='P';

    middleStandardGeneticCode["CCR"]='P';
    middleStandardGeneticCode["CCY"]='P';
    middleStandardGeneticCode["CCS"]='P';
    middleStandardGeneticCode["CCW"]='P';
    middleStandardGeneticCode["CCK"]='P';
    middleStandardGeneticCode["CCM"]='P';
    middleStandardGeneticCode["CCB"]='P';
    middleStandardGeneticCode["CCD"]='P';
    middleStandardGeneticCode["CCH"]='P';
    middleStandardGeneticCode["CCV"]='P';
    middleStandardGeneticCode["CCN"]='P';

    ///////////////////H
    middleStandardGeneticCode["CAT"]='H';
    middleStandardGeneticCode["CAC"]='H';
    middleStandardGeneticCode["CAY"]='H';

    ///////////////////Q
    middleStandardGeneticCode["CAA"]='Q';
    middleStandardGeneticCode["CAG"]='Q';
    middleStandardGeneticCode["CAR"]='Q';

    ///////////////////R
    middleStandardGeneticCode["CGT"]='R';
    middleStandardGeneticCode["CGC"]='R';
    middleStandardGeneticCode["CGA"]='R';
    middleStandardGeneticCode["CGG"]='R';

    middleStandardGeneticCode["CGR"]='R';
    middleStandardGeneticCode["CGY"]='R';
    middleStandardGeneticCode["CGS"]='R';
    middleStandardGeneticCode["CGW"]='R';
    middleStandardGeneticCode["CGK"]='R';
    middleStandardGeneticCode["CGM"]='R';
    middleStandardGeneticCode["CGB"]='R';
    middleStandardGeneticCode["CGD"]='R';
    middleStandardGeneticCode["CGH"]='R';
    middleStandardGeneticCode["CGV"]='R';
    middleStandardGeneticCode["CGN"]='R';

    middleStandardGeneticCode["AGA"]='R';
    middleStandardGeneticCode["AGG"]='R';
    middleStandardGeneticCode["AGR"]='R';

    middleStandardGeneticCode["MGA"]='R';
    middleStandardGeneticCode["MGG"]='R';
    middleStandardGeneticCode["MGR"]='R';



    ///////////////////I
    middleStandardGeneticCode["ATT"]='I';
    middleStandardGeneticCode["ATC"]='I';
    middleStandardGeneticCode["ATA"]='I';

    middleStandardGeneticCode["ATY"]='I';
    middleStandardGeneticCode["ATW"]='I';
    middleStandardGeneticCode["ATM"]='I';

    ///////////////////M
    middleStandardGeneticCode["ATG"]='M';

    ///////////////////T
    middleStandardGeneticCode["ACT"]='T';
    middleStandardGeneticCode["ACC"]='T';
    middleStandardGeneticCode["ACA"]='T';
    middleStandardGeneticCode["ACG"]='T';

    middleStandardGeneticCode["ACR"]='T';
    middleStandardGeneticCode["ACY"]='T';
    middleStandardGeneticCode["ACS"]='T';
    middleStandardGeneticCode["ACW"]='T';
    middleStandardGeneticCode["ACK"]='T';
    middleStandardGeneticCode["ACM"]='T';
    middleStandardGeneticCode["ACB"]='T';
    middleStandardGeneticCode["ACD"]='T';
    middleStandardGeneticCode["ACH"]='T';
    middleStandardGeneticCode["ACV"]='T';
    middleStandardGeneticCode["ACN"]='T';

    ///////////////////N
    middleStandardGeneticCode["AAT"]='N';
    middleStandardGeneticCode["AAC"]='N';
    middleStandardGeneticCode["AAY"]='N';

    ///////////////////K
    middleStandardGeneticCode["AAA"]='K';
    middleStandardGeneticCode["AAG"]='K';

    middleStandardGeneticCode["AAR"]='K';

    ///////////////////V
    middleStandardGeneticCode["GTT"]='V';
    middleStandardGeneticCode["GTC"]='V';
    middleStandardGeneticCode["GTA"]='V';
    middleStandardGeneticCode["GTG"]='V';

    middleStandardGeneticCode["GTR"]='V';
    middleStandardGeneticCode["GTY"]='V';
    middleStandardGeneticCode["GTS"]='V';
    middleStandardGeneticCode["GTW"]='V';
    middleStandardGeneticCode["GTK"]='V';
    middleStandardGeneticCode["GTM"]='V';
    middleStandardGeneticCode["GTB"]='V';
    middleStandardGeneticCode["GTD"]='V';
    middleStandardGeneticCode["GTH"]='V';
    middleStandardGeneticCode["GTV"]='V';
    middleStandardGeneticCode["GTN"]='V';

    ///////////////////A
    middleStandardGeneticCode["GCT"]='A';
    middleStandardGeneticCode["GCC"]='A';
    middleStandardGeneticCode["GCA"]='A';
    middleStandardGeneticCode["GCG"]='A';

    middleStandardGeneticCode["GCR"]='A';
    middleStandardGeneticCode["GCY"]='A';
    middleStandardGeneticCode["GCS"]='A';
    middleStandardGeneticCode["GCW"]='A';
    middleStandardGeneticCode["GCK"]='A';
    middleStandardGeneticCode["GCM"]='A';
    middleStandardGeneticCode["GCB"]='A';
    middleStandardGeneticCode["GCD"]='A';
    middleStandardGeneticCode["GCH"]='A';
    middleStandardGeneticCode["GCV"]='A';
    middleStandardGeneticCode["GCN"]='A';

    ///////////////////D
    middleStandardGeneticCode["GAT"]='D';
    middleStandardGeneticCode["GAC"]='D';
    middleStandardGeneticCode["GAY"]='D';

    ///////////////////E
    middleStandardGeneticCode["GAA"]='E';
    middleStandardGeneticCode["GAG"]='E';
    middleStandardGeneticCode["GAR"]='E';

    ///////////////////G
    middleStandardGeneticCode["GGT"]='G';
    middleStandardGeneticCode["GGC"]='G';
    middleStandardGeneticCode["GGA"]='G';
    middleStandardGeneticCode["GGG"]='G';

    middleStandardGeneticCode["GGR"]='G';
    middleStandardGeneticCode["GGY"]='G';
    middleStandardGeneticCode["GGS"]='G';
    middleStandardGeneticCode["GGW"]='G';
    middleStandardGeneticCode["GGK"]='G';
    middleStandardGeneticCode["GGM"]='G';
    middleStandardGeneticCode["GGB"]='G';
    middleStandardGeneticCode["GGD"]='G';
    middleStandardGeneticCode["GGH"]='G';
    middleStandardGeneticCode["GGV"]='G';
    middleStandardGeneticCode["GGN"]='G';

    middleStandardGeneticCode["---"]='-';



    legalNasString.insert("A");
    legalNasString.insert("a");
    legalNasString.insert("T");
    legalNasString.insert("t");
    legalNasString.insert("C");
    legalNasString.insert("c");
    legalNasString.insert("G");
    legalNasString.insert("g");
    legalNasString.insert("U");
    legalNasString.insert("u");

    legalNasChar.insert('A');
    legalNasChar.insert('a');
    legalNasChar.insert('T');
    legalNasChar.insert('t');
    legalNasChar.insert('C');
    legalNasChar.insert('c');
    legalNasChar.insert('G');
    legalNasChar.insert('g');
    legalNasChar.insert('U');
    legalNasChar.insert('u');

    std::string exepath = getexepath(parameters);
    std::string spliceSitesFile = exepath+"/SpliceSites";
    std::ifstream infile(spliceSitesFile);
    if( ! infile.good()){
        std::cerr << "error in opening splice sites file " << spliceSitesFile << std::endl;
        exit (1);
    }

    std::regex reg("^(\\S+)\\s+(\\S+)");
    std::string line="";
    while (std::getline(infile, line)) {
        if (line.compare(0, 1, "#") == 0) {
            continue;
        }
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty() ) {
        } else {
            std::string dornor = match[1];
            std::string acceptor = match[2];
            dornors.emplace_back(dornor);
            acceptors.emplace_back(acceptor);
            std::set<std::string> allDornors;
            std::set<std::string> allAcceptors;

            getAllPossibleWithIupac(dornor, allDornors);
            getAllPossibleWithIupac(acceptor, allAcceptors);

            allPossiableDornors[dornor] = allDornors;
            allPossiableAcceptors[acceptor] = allAcceptors;
        }
    }
}

NucleotideCodeSubstitutionMatrix::~NucleotideCodeSubstitutionMatrix(){
/*    int row = sizeof(nucleotide_substitution_matrix);
    int col = sizeof(nucleotide_substitution_matrix[0]);
    for (int i = 0; i < row; i++)
    {
        for( int j=0; j< col; j++ ){
            delete[] nucleotide_substitution_matrix[i][j];
        }
        delete[] this->nucleotide_substitution_matrix[i];
    }
    delete[] this->nucleotide_substitution_matrix;
*/
    for( int i=0; i<5; ++i ){
        delete[] _exon_subsitition_matrix[i];
        delete[] _intron_subsitition_matrix[i];
        delete[] _start_stop_codon_subsitition_matrix[i];
        delete[] _splice_sites_subsitition_matrix[i];
    }
    delete[] _exon_subsitition_matrix;
    delete[] _intron_subsitition_matrix;
    delete[] _start_stop_codon_subsitition_matrix;
    delete[] _splice_sites_subsitition_matrix;
}
int8_t  NucleotideCodeSubstitutionMatrix::get_dna_acid_map( const char & c ){
    if( dna_acid_map.find(c) == dna_acid_map.end() ){
        return 4;
    }else{
        return  dna_acid_map[c];
    }
}

//std::map<char, int>& NucleotideCodeSubstitutionMatrix::get_dna_acid_map(){
//    return dna_acid_map;
//}
//std::map<std::string, double>& NucleotideCodeSubstitutionMatrix::get_allStartCodons(){
//    return allStartCodons;
//}
//std::map<std::string, double>& NucleotideCodeSubstitutionMatrix::get_allStopCodons(){
//    return allStopCodons;
//}
std::set<std::string>& NucleotideCodeSubstitutionMatrix::getMustStartCodons() {
    return mustStartCodons;
}
std::set<std::string>& NucleotideCodeSubstitutionMatrix::getMustStopCodons() {
    return mustStopCodons;
}

std::set<std::string>& NucleotideCodeSubstitutionMatrix::getPossibleStartCodons() {
    return possibleStartCodons;
}
std::set<std::string>& NucleotideCodeSubstitutionMatrix::getPossibleStopCodons() {
    return possibleStopCodons;
}
std::set<std::string>& NucleotideCodeSubstitutionMatrix::getBasicStartCodons() {
    return basicStartCodons;
}
std::set<std::string>& NucleotideCodeSubstitutionMatrix::getBasicStopCodons() {
    return basicStopCodons;
}
std::map<std::string, std::set<std::string> >& NucleotideCodeSubstitutionMatrix::getDnaIupacCode(){
    return this->dnaIupacCode;
}

//
//void NucleotideCodeSubstitutionMatrix::getAllPossibleWithIupac(const std::string& seq, std::set<std::string>& currentPossibleCombinations) {
//
//    for( size_t i=0; i<seq.size(); ++i ){
//        if( i > 0 ) {
//            std::set<std::string> tempPossibleCombinations;
//            for (std::set<std::string>::iterator it = revDnaIupacCode[seq.substr(i,1)].begin(); it != revDnaIupacCode[seq.substr(i,1)].end(); ++it) {
//                for( std::set<std::string>::iterator it1=currentPossibleCombinations.begin(); it1!=currentPossibleCombinations.end(); ++it1  ){
//                    tempPossibleCombinations.insert( (*it1) + (*it));
//                }
//            }
//            currentPossibleCombinations = tempPossibleCombinations;
//        }else{
//            for (std::set<std::string>::iterator it = revDnaIupacCode[seq.substr(i,1)].begin();
//                 it != revDnaIupacCode[seq.substr(i,1)].end(); ++it) {
//                currentPossibleCombinations.insert((*it));
//            }
//        }
//    }
//    return;
//}


void NucleotideCodeSubstitutionMatrix::getAllPossibleWithIupac(const std::string& seq, std::set<std::string>& currentPossibleCombinations) {
    std::vector< std::set<std::string> > allPossibleForEachChar;
    for( size_t i=0; i < seq.size(); ++i){
        std::set<std::string> allPossibleForThisChar;
        std::string thisChar = seq.substr(i, 1);
        for( std::map<std::string, std::set<std::string> >::iterator it1=dnaIupacCode.begin(); it1!=dnaIupacCode.end(); ++it1){
            if( it1->second.find(thisChar) != it1->second.end() ){
                std::string thisIupac = it1->first;
                allPossibleForThisChar.insert(thisIupac);
            }
        }
        allPossibleForEachChar.push_back(allPossibleForThisChar);
    }


    for( size_t i=0; i<allPossibleForEachChar.size(); ++i ){
        if( currentPossibleCombinations.size() > 0 ) {
            std::set<std::string> tempPossibleCombinations;
            for (std::set<std::string>::iterator it = allPossibleForEachChar[i].begin(); it != allPossibleForEachChar[i].end(); ++it) {
                for( std::set<std::string>::iterator it1=currentPossibleCombinations.begin(); it1!=currentPossibleCombinations.end(); ++it1  ){
                    tempPossibleCombinations.insert( (*it1) + (*it));
                }
            }
            currentPossibleCombinations = tempPossibleCombinations;
        }else{
            for (std::set<std::string>::iterator it = allPossibleForEachChar[i].begin();
                it != allPossibleForEachChar[i].end(); ++it) {
                currentPossibleCombinations.insert((*it));
            }
        }
    }
    return;
}


char NucleotideCodeSubstitutionMatrix::getGeneticCode(const std::string & coding) {
    return getGeneticCode(coding, MIDDLE);
}
char NucleotideCodeSubstitutionMatrix::getGeneticCode(const std::string & coding, const BEGINMIDDLEEND & beginmiddleend ){
    switch (beginmiddleend){
        case MIDDLE:
            if(middleStandardGeneticCode.find(coding)!=middleStandardGeneticCode.end()){
                return middleStandardGeneticCode[coding];
            }else{
                return 'X';
            }
            break;
        case BEGIN:
            if( possibleStartCodons.find(coding)!=possibleStartCodons.end() ){
                return 'M';
            }else if(middleStandardGeneticCode.find(coding)!=middleStandardGeneticCode.end()){
                return middleStandardGeneticCode[coding];
            }else{
                return 'X';
            }
            break;
        case END:
            if( possibleStopCodons.find(coding)!=possibleStopCodons.end() ){
                return '*';
            }else if(middleStandardGeneticCode.find(coding)!=middleStandardGeneticCode.end()){
                return middleStandardGeneticCode[coding];
            }else{
                return 'X';
            }
            break;
        default:
            return 'X';
            break;
    }
}

int8_t** NucleotideCodeSubstitutionMatrix::get_exon_subsitition_matrix() const{
    return _exon_subsitition_matrix;
}
int8_t** NucleotideCodeSubstitutionMatrix::get_intron_subsitition_matrix() const{
    return _intron_subsitition_matrix;
}
int8_t** NucleotideCodeSubstitutionMatrix::get_start_stop_codon_subsitition_matrix() const{
    return _start_stop_codon_subsitition_matrix;
}
int8_t** NucleotideCodeSubstitutionMatrix::get_splice_sites_subsitition_matrix() const{
    return _splice_sites_subsitition_matrix;
}
std::set<std::string>& NucleotideCodeSubstitutionMatrix::getLegalNasString() {
    return legalNasString;
}
std::set<char>& NucleotideCodeSubstitutionMatrix::getLegalNasChar() {
    return legalNasChar;
}
std::vector<std::string>& NucleotideCodeSubstitutionMatrix::getDornors() {
    return dornors;
}
std::vector<std::string>& NucleotideCodeSubstitutionMatrix::getAcceptors() {
    return acceptors;
}

std::set<std::string> NucleotideCodeSubstitutionMatrix::getAllPossiableDornors(std::string& dornor){
    return allPossiableDornors[dornor];
}
std::set<std::string> NucleotideCodeSubstitutionMatrix::getAllPossiableAcceptors(std::string& acceptor){
    return allPossiableAcceptors[acceptor];
}

/*  * IUPAC codes
 *  DNA:
 *
 *  Nucleotide Code:  Base:
 *  ----------------  -----
 *  A.................Adenine
 *  C.................Cytosine
 *  G.................Guanine
 *  T (or U)..........Thymine (or Uracil)
 *  R.................A or G
 *  Y.................C or T
 *  S.................G or C
 *  W.................A or T
 *  K.................G or T
 *  M.................A or C
 *  B.................C or G or T
 *  D.................A or G or T
 *  H.................A or C or T
 *  V.................A or C or G
 *  N.................any base
 *  . or -............gap */
