// =====================================================================================
//
//       Filename:  Transcript.cpp
//
//    Description:
//
//        Version:  1.0
//        Created:  04/12/2017 02:13:56 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//
// =====================================================================================

#include "Transcript.h"

Transcript::Transcript(const std::string &name, const std::string &chromeSomeName, const STRAND &strand) {
    this->_name = name;
    this->_chromeSomeName = chromeSomeName;
    this->_strand = strand;
    this->_start = std::numeric_limits<int>::max();
    this->_end = 0;
}

Transcript::Transcript() {
    this->_strand = POSITIVE;
    this->_start = std::numeric_limits<int>::max();
    this->_end = 0;
}

const std::string &Transcript::getName() const {
    return _name;
}

void Transcript::setName(const std::string &name) {
    this->_name = name;
}

std::set<GenomeBasicFeature> &Transcript::getCdsHashSet() {
    return _cdsHashSet;
}

std::vector<GenomeBasicFeature> &Transcript::getCdsVector() {
    return _cdsVector;
}

std::set<GenomeBasicFeature> &Transcript::getExonHashSet() {
    return this->_exonHashSet;
}

std::vector<GenomeBasicFeature> &Transcript::getExonVector() {
    return this->_exonVector;
}

std::vector<GenomeBasicFeature> &Transcript::getThreePrimerUtr() {
    return this->_threePrimerUtrVector;
}

void Transcript::addThreePrimerUtr(const GenomeBasicFeature &threePrimerUtr) {
    this->_threePrimerUtrVector.push_back(threePrimerUtr);
    this->_ifThreePrimerUtr = true;
}

std::vector<GenomeBasicFeature> &Transcript::getFivePrimerUtr() {
    return this->_fivePrimerUtrVector;
}

void Transcript::addFivePrimerUtr(const GenomeBasicFeature &fivePrimerUtr) {
    this->_fivePrimerUtrVector.push_back(fivePrimerUtr);
    this->_ifFivePrimerUtr = true;
}

const STRAND &Transcript::getStrand() const {
    return _strand;
}

const std::string &Transcript::getChromeSomeName() const {
    return _chromeSomeName;
}

void Transcript::addCds(const GenomeBasicFeature &cds) {
    this->_cdsHashSet.insert(cds);
    this->_cdsVector.push_back(cds);
}

void Transcript::addExon(const GenomeBasicFeature &exon) {
    this->_exonHashSet.insert(exon);
    this->_exonVector.push_back(exon);
}

void Transcript::setStart(const int &start) {
    this->_start = start;
}

void Transcript::setEnd(const int &end) {
    this->_end = end;
}

const int &Transcript::getStart() const {
    return this->_start;
}

const int &Transcript::getEnd() const {
    return this->_end;
}

const int &Transcript::getPStart() const {
    return _pStart;
}

const int &Transcript::getPEnd() const {
    return _pEnd;
}

void Transcript::setPStart(const int &pStart) {
    this->_pStart = pStart;
}

void Transcript::setPEnd(const int &pEnd) {
    this->_pEnd = pEnd;
}

void Transcript::updateInforCds() {
    std::sort(this->_cdsVector.begin(), this->_cdsVector.end(), [](GenomeBasicFeature a, GenomeBasicFeature b) {
        return a < b;
    });
    this->_pStart = this->_cdsVector[0].getStart();
    this->_pEnd = this->_cdsVector[this->_cdsVector.size() - 1].getEnd();
}

const std::string &Transcript::getGeneomeSequence() const {
    return this->_geneomeSequence;
}

void Transcript::setGeneomeSequence(const std::string &genomeSequence) {
    this->_geneomeSequence = genomeSequence;
}

const std::string &Transcript::getCdsSequence() const {
    return this->_cdsSequence;
}

void Transcript::setCdsSequence(const std::string &cdsSequence) {
    this->_cdsSequence = cdsSequence;
}

const std::string &Transcript::getExonSequence() const {
    return this->_exonSequence;
}

void Transcript::setExonSequence(const std::string &exonSequence) {
    this->_exonSequence = exonSequence;
}

void Transcript::setMetaInformation(const std::string &metaInformation) {
    this->_metaInformation = metaInformation;
}

const std::string &Transcript::getMetaInformation() const {
    return this->_metaInformation;
}

const bool &Transcript::getIfOrfShift() const {
    return _ifOrfShift;
}

const std::string &Transcript::getScore() const {
    return this->_score;
}

void Transcript::setScore(const std::string &score) {
    this->_score = score;
}

const std::string &Transcript::getLastColumnInformation() const {
    return this->_lastColumnInformation;
}

void Transcript::setLastColumnInformation(const std::string &lastColumnInformation) {
    this->_lastColumnInformation = lastColumnInformation;
}

const std::string &Transcript::getType() const {
    return this->_type;
}

void Transcript::setType(const std::string &type) {
    this->_type = type;
}

void Transcript::cleanCdsVector() {
    this->_cdsVector.clear();
}

bool Transcript::ifOverLapIgnorStrand(const Transcript &transcript) const {
    if (this->_pStart <= transcript._pStart &&
        transcript._pStart <= this->_pEnd) {
        return true;
    }
    if (this->_pStart <= transcript._pEnd &&
        transcript._pEnd <= this->_pEnd) {
        return true;
    }
    if (transcript._pStart <= this->_pStart &&
        this->_pStart <= transcript._pEnd) {
        return true;
    }
    if (transcript._pStart <= this->_pEnd &&
        this->_pEnd <= transcript._pEnd) {
        return true;
    }
    return false;
}

bool Transcript::ifOverLap(const Transcript &transcript) const {
    if (this->_strand == transcript._strand && this->_chromeSomeName.compare(transcript.getChromeSomeName()) == 0) {
        if (this->_pStart <= transcript._pStart &&
            transcript._pStart <= this->_pEnd) {
            return true;
        }
        if (this->_pStart <= transcript._pEnd &&
            transcript._pEnd <= this->_pEnd) {
            return true;
        }
        if (transcript._pStart <= this->_pStart &&
            this->_pStart <= transcript._pEnd) {
            return true;
        }
        if (transcript._pStart <= this->_pEnd &&
            this->_pEnd <= transcript._pEnd) {
            return true;
        }
    }
    return false;
}

void Transcript::setIfOrfShift(const bool &ifOrfShift) {
    this->_ifOrfShift = ifOrfShift;
}

const std::string &Transcript::getSource() const {
    return this->_source;
}

void Transcript::setSource(const std::string source) {
    this->_source = source;
}

bool Transcript::is_ifFivePrimerUtr() const {
    return _ifFivePrimerUtr;
}

void Transcript::set_ifFivePrimerUtr(bool _ifFivePrimerUtr) {
    Transcript::_ifFivePrimerUtr = _ifFivePrimerUtr;
}

bool Transcript::is_ifThreePrimerUtr() const {
    return _ifThreePrimerUtr;
}

void Transcript::set_ifThreePrimerUtr(bool _ifThreePrimerUtr) {
    Transcript::_ifThreePrimerUtr = _ifThreePrimerUtr;
}
