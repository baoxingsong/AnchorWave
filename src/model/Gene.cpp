/*
 * =====================================================================================
 *
 *       Filename:  Gene.cpp
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




**************************************************************************/


#include "Gene.h"
//gene begin
Gene::Gene(const std::string& name, const STRAND& strand){
    this->_name=name;
    this->_strand=strand;
}
Gene::Gene(const std::string& name, const std::string & chromeSomeName, const STRAND& strand, const int & start, const int & end){
    this->_start=start;
    this->_end=end;
    this->_name=name;
    this->_strand=strand;
    this->_chromeSomeName = chromeSomeName;
}
Gene::Gene(const std::string& name, const std::string & chromeSomeName, const STRAND& strand){
    this->_name=name;
    this->_strand=strand;
    this->_chromeSomeName = chromeSomeName;
}
Gene::Gene() {
}
const std::string & Gene::getName(){
    return this->_name;
}
const STRAND & Gene::getStrand(){
    return this->_strand;
}
const int & Gene::getStart() const{
    return this->_start;
}
void Gene::setStart( const int & start){
    this->_start = start;
}
const int & Gene::getEnd() const{
    return this->_end;
}
void Gene::setEnd( const int & end){
    this->_end=end;
}
const std::string & Gene::getChromeSomeName() const{
    return this->_chromeSomeName;
}
std::vector<Transcript>& Gene::getTranscripts(){
    return this->_transcripts;
}
std::vector<std::string>& Gene::getTranscriptVector(){
    return this->_transcriptVector;
}

void Gene::addTranscript(const Transcript& transcript){
    this->_transcriptVector.push_back(transcript.getName());
    this->updateStartEnd(transcript);
}
void Gene::updateStartEnd( const Transcript& transcript ){
    if( _start > transcript.getStart()  ){
        _start = transcript.getStart();
    }
    if(_end < transcript.getEnd()){
        _end=transcript.getEnd();
    }
}
std::string & Gene::getSource(){
    return this->_source;
}
void Gene::setSource(const std::string& source){
    this->_source=source;
}

const std::string &Gene::getScore() const {
    return _score;
}
void Gene::setScore(const std::string &score) {
    this->_score = score;
}

const std::string &Gene::getCodonFrame() const {
    return _codonFrame;
}
void Gene::setCodonFrame(const std::string &codonFrame) {
    this->_codonFrame = codonFrame;
}

const std::string &Gene::getLastColumnInformation() const {
    return _lastColumnInformation;
}

void Gene::setLastColumnInformation(const std::string &lastColumnInformation) {
    this->_lastColumnInformation = lastColumnInformation;
}
//gene end
