/*
 * =====================================================================================
 *
 *       Filename:  GenomeBasicFeature.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2017 15:47:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************

no matter which strand it is on, it always start with the smaller coordinate

 ************************************************************************/

#include "GenomeBasicFeature.h"

GenomeBasicFeature::GenomeBasicFeature(const int &start, const int &end) {
    if (start < end) {
        this->_start = start;
        this->_end = end;
    } else {
        this->_end = start;
        this->_start = end;
    }
}

GenomeBasicFeature::GenomeBasicFeature(const int &start, const int &end, const std::string &score, const std::string &codonFrame, const std::string &lastColumnInformation) {
    if (start < end) {
        this->_start = start;
        this->_end = end;
    } else {
        this->_end = start;
        this->_start = end;
    }
    this->_score = score;
    this->_codonFrame = codonFrame;
    this->_lastColumnInformation = lastColumnInformation;
}

GenomeBasicFeature::GenomeBasicFeature() {

}

const int &GenomeBasicFeature::getStart() const {
    return _start;
}

void GenomeBasicFeature::setStart(const int &start) {
    this->_start = start;
}

const int &GenomeBasicFeature::getEnd() const {
    return _end;
}

void GenomeBasicFeature::setEnd(const int &end) {
    this->_end = end;
}

const std::string &GenomeBasicFeature::getScore() const {
    return _score;
}

void GenomeBasicFeature::setScore(const std::string &score) {
    this->_score = score;
}

const std::string &GenomeBasicFeature::getCodonFrame() const {
    return _codonFrame;
}

void GenomeBasicFeature::setCodonFrame(const std::string &codonFrame) {
    this->_codonFrame = codonFrame;
}

const std::string &GenomeBasicFeature::getLastColumnInformation() const {
    return _lastColumnInformation;
}

void GenomeBasicFeature::setLastColumnInformation(const std::string &lastColumnInformation) {
    this->_lastColumnInformation = lastColumnInformation;
}

const std::string &GenomeBasicFeature::getType() const {
    return type;
}

void GenomeBasicFeature::setType(const std::string &type) {
    GenomeBasicFeature::type = type;
}
