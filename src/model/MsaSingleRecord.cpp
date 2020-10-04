// =====================================================================================
//
//       Filename:  MsaSingleRecord.cpp
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



#include "MsaSingleRecord.h"

MsaSingleRecord::MsaSingleRecord(const int & _start, const int & _end, const std::string & _accesionName, const std::string & _sequence){
    this->start = _start;
    this->end = _end;
    this->accesionName = _accesionName;
    this->sequence = _sequence;
}

MsaSingleRecord::MsaSingleRecord(){

}
const int & MsaSingleRecord::getStart() const {
    return this->start;
}
const int & MsaSingleRecord::getEnd() const{
    return this->end;
}
void MsaSingleRecord::setStart( const int & _start ){
    this->start = _start;
}
void MsaSingleRecord::setEnd( const int & _end ){
    this->end = _end;
}
const std::string & MsaSingleRecord::getSequence(){
    return this->sequence;
}