/*
 * =====================================================================================
 *
 *       Filename:  nucleotideCodeSubstitutionMatrix.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/03/2017 23:33:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

#ifndef _NUCLEOTIDECODESUBSTITUTIONMATRIX_H
#define _NUCLEOTIDECODESUBSTITUTIONMATRIX_H

#include <map>
#include <set>
#include <string>
#include <vector>

enum VARIANTCATEGORY {
    SNP, INSERTION, DELETION, SNPORINSERTION, SNPORDELETION, INSERTIONORDELETION, SNPORINSERTIONORDELETION
};


enum BEGINMIDDLEEND //  FOR DNA TO PROTEIN TRANSLATION
{
    BEGIN, MIDDLE, END
};

class NucleotideCodeSubstitutionMatrix {
private:
    std::map<char, int> dna_acid_map;
    std::set<std::string> mustStartCodons;
    std::set<std::string> mustStopCodons;
    std::set<std::string> possibleStartCodons;
    std::set<std::string> possibleStopCodons;
    std::set<std::string> basicStartCodons;
    std::set<std::string> basicStopCodons;

    std::vector<std::string> dornors;
    std::vector<std::string> acceptors;


    std::map<std::string, std::set<std::string> > dnaIupacCode;
    std::map<std::string, std::set<std::string> > revDnaIupacCode;
    std::map<std::string, char> middleStandardGeneticCode;


    std::set<std::string> legalNasString;
    std::set<char> legalNasChar;
public:
    NucleotideCodeSubstitutionMatrix();

    void initial();

    NucleotideCodeSubstitutionMatrix(const int32_t &matchingScore, const int32_t &mismatchingPenalty);

    int32_t get_dna_acid_map(const char &c);

    //std::map<char, int>& get_dna_acid_map();
    //std::map<std::string, double>& get_allStartCodons();
    //std::map<std::string, double>& get_allStopCodons();
    std::set<std::string> &getMustStartCodons();

    std::set<std::string> &getMustStopCodons();

    std::set<std::string> &getPossibleStartCodons();

    std::set<std::string> &getPossibleStopCodons();

    std::set<std::string> &getBasicStartCodons();

    std::set<std::string> &getBasicStopCodons();

    //give a sequence, return all the possiable iupac combinations
    void getAllPossibleWithIupac(const std::string &seq, std::set<std::string> &currentPossibleCombinations);

    std::map<std::string, std::set<std::string> > &getDnaIupacCode();

    int32_t nucleotide_substitution_matrix[5][5];

    char getGeneticCode(const std::string &coding, const BEGINMIDDLEEND &beginmiddleend);

    char getGeneticCode(const std::string &coding);

    std::set<std::string> &getLegalNasString();

    std::set<char> &getLegalNasChar();

    std::vector<std::string> &getDornors();

    std::vector<std::string> &getAcceptors();
};

#endif
//std::unordered_set<std::string> stopCodons = {"ATG", "TTG", "CTG"};// = {"TAA", "TAG", "TGA"}; // TO DO, use iupac code


