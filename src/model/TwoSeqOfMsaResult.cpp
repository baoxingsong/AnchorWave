// =====================================================================================
//
//       Filename:  TwoSeqOfMsaResult.cpp
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
/*************************************************************************





*************************************************************************/




#include "TwoSeqOfMsaResult.h"

TwoSeqOfMsaResult::TwoSeqOfMsaResult(const int & _refStart, const int&  _refEnd, const std::string & _refSeq, const int & _resultStart, const int & _resultEnd, const std::string & _resultSeq){
    this->refStart = _refStart;
    this->refEnd =_refEnd;
    this->refSeq = _refSeq;
    this->resultStart = _resultStart;
    this->resultEnd = _resultEnd;
    this->resultSeq = _resultSeq;
}
const int & TwoSeqOfMsaResult::getRefStart() const{
    return this->refStart;
}
void TwoSeqOfMsaResult::setRefStart( const int & _refStart){
    this->refStart=_refStart;
}
const int & TwoSeqOfMsaResult::getRefEnd() const{
    return this->refEnd;
}
void TwoSeqOfMsaResult::setRefEnd( const int & _refEnd){
    this->refEnd = _refEnd;
}
const std::string & TwoSeqOfMsaResult::getRefSeq() const{
    return this->refSeq;
}
void TwoSeqOfMsaResult::setRefSeq(const std::string & _refSequence){
    this->refSeq = _refSequence;
}
const int & TwoSeqOfMsaResult::getResultEnd() const{
    return this->resultEnd;
}
void TwoSeqOfMsaResult::setResultEnd( const int & _resultEnd){
    this->resultEnd = _resultEnd;
}
const int & TwoSeqOfMsaResult::getResultStart() const{
    return this->resultStart;
}
void TwoSeqOfMsaResult::setResultStart(const int & _resultStart){
    this->resultStart=_resultStart;
}
const std::string & TwoSeqOfMsaResult::getResultSeq() const{
    return this->resultSeq;
}
void TwoSeqOfMsaResult::setResultSeq(const std::string & _resultSequence){
    this->resultSeq=_resultSequence;
}
// msa link data sructure

