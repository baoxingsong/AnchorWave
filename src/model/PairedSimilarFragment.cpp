//
// Created by Baoxing song on 2019-01-02.
//

#include "PairedSimilarFragment.h"
PairedSimilarFragment::PairedSimilarFragment(){
    start1=0;
    start2=0;
    end1=0;
    end2=0;
    score=0;
}

PairedSimilarFragment::PairedSimilarFragment(uint32_t start1, uint32_t end1, uint32_t start2, uint32_t end2,
        uint32_t score, const std::vector<uint32_t> &cigar, double & pValue, double & eValue): start1(start1),
    end1(end1), start2(start2), end2(end2), score(score), cigar(cigar), pValue(pValue), eValue(eValue){

}

uint32_t PairedSimilarFragment::getStart1() const {
    return start1;
}
void PairedSimilarFragment::setStart1(uint32_t start1) {
    PairedSimilarFragment::start1 = start1;
}
uint32_t PairedSimilarFragment::getEnd1() const {
    return end1;
}
void PairedSimilarFragment::setEnd1(uint32_t end1) {
    PairedSimilarFragment::end1 = end1;
}
uint32_t PairedSimilarFragment::getStart2() const {
    return start2;
}
void PairedSimilarFragment::setStart2(uint32_t start2) {
    PairedSimilarFragment::start2 = start2;
}
uint32_t PairedSimilarFragment::getEnd2() const {
    return end2;
}
void PairedSimilarFragment::setEnd2(uint32_t end2) {
    PairedSimilarFragment::end2 = end2;
}
uint32_t PairedSimilarFragment::getScore() const {
    return score;
}
uint32_t PairedSimilarFragment::getLength() const{
    return end2 - start2 + 1;
}
void PairedSimilarFragment::setScore(uint32_t score) {
    PairedSimilarFragment::score = score;
}
const std::vector<uint32_t> &PairedSimilarFragment::getCigar() const {
    return cigar;
}
void PairedSimilarFragment::setCigar(const std::vector<uint32_t> &cigar) {
    PairedSimilarFragment::cigar = cigar;
}
double PairedSimilarFragment::getPValue() const {
    return pValue;
}

void PairedSimilarFragment::setPValue(double pValue) {
    PairedSimilarFragment::pValue = pValue;
}

double PairedSimilarFragment::getEValue() const {
    return eValue;
}

void PairedSimilarFragment::setEValue(double eValue) {
    PairedSimilarFragment::eValue = eValue;
}

const std::string &PairedSimilarFragment::getAlignment1() const {
    return alignment1;
}

void PairedSimilarFragment::setAlignment1(const std::string &alignment1) {
    PairedSimilarFragment::alignment1 = alignment1;
}

const std::string &PairedSimilarFragment::getAlignment2() const {
    return alignment2;
}

void PairedSimilarFragment::setAlignment2(const std::string &alignment2) {
    PairedSimilarFragment::alignment2 = alignment2;
}

PairedSimilarFragment2::PairedSimilarFragment2(const std::string &species, const std::string &queryChr, uint32_t start1,
                                               uint32_t end1, uint32_t start2, uint32_t end2,
                                               const std::string &alignment1, const std::string &alignment2,
                                               int8_t strand) : species(species), queryChr(queryChr), start1(start1),
                                                                end1(end1), start2(start2), end2(end2),
                                                                alignment1(alignment1), alignment2(alignment2),
                                                                strand(strand) {}

const std::string &PairedSimilarFragment2::getSpecies() const {
    return species;
}

void PairedSimilarFragment2::setSpecies(const std::string &species) {
    PairedSimilarFragment2::species = species;
}

const std::string &PairedSimilarFragment2::getQueryChr() const {
    return queryChr;
}

void PairedSimilarFragment2::setQueryChr(const std::string &queryChr) {
    PairedSimilarFragment2::queryChr = queryChr;
}

uint32_t PairedSimilarFragment2::getStart1() const {
    return start1;
}

void PairedSimilarFragment2::setStart1(uint32_t start1) {
    PairedSimilarFragment2::start1 = start1;
}

uint32_t PairedSimilarFragment2::getEnd1() const {
    return end1;
}

void PairedSimilarFragment2::setEnd1(uint32_t end1) {
    PairedSimilarFragment2::end1 = end1;
}

uint32_t PairedSimilarFragment2::getStart2() const {
    return start2;
}

void PairedSimilarFragment2::setStart2(uint32_t start2) {
    PairedSimilarFragment2::start2 = start2;
}

uint32_t PairedSimilarFragment2::getEnd2() const {
    return end2;
}

void PairedSimilarFragment2::setEnd2(uint32_t end2) {
    PairedSimilarFragment2::end2 = end2;
}

const std::string &PairedSimilarFragment2::getAlignment1() const {
    return alignment1;
}

void PairedSimilarFragment2::setAlignment1(const std::string &alignment1) {
    PairedSimilarFragment2::alignment1 = alignment1;
}

const std::string &PairedSimilarFragment2::getAlignment2() const {
    return alignment2;
}

void PairedSimilarFragment2::setAlignment2(const std::string &alignment2) {
    PairedSimilarFragment2::alignment2 = alignment2;
}

int8_t PairedSimilarFragment2::getStrand() const {
    return strand;
}

void PairedSimilarFragment2::setStrand(int8_t strand) {
    PairedSimilarFragment2::strand = strand;
}
