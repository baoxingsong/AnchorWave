// =====================================================================================
//
//       Filename:  MsaFileRecord.h
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

The MSA based variants recalling pipeline, cut the genome sequence into windows firstly,
    and them perform MSA, and recall variants

The class is designed for the MSA result. It has a windows start and windows end.

 windows start - overlap size = reference start
 windows end + overlap size = reference end

The class has window start and window end, and a group of single sequence records.
 This reference start, reference end and non-reference start and non-reference end are variantble in single sequnce records


 ************************************************************************/


#ifndef ANNOTATIONLIFTOVER_MSAFILERECORD_H
#define ANNOTATIONLIFTOVER_MSAFILERECORD_H

#include "MsaSingleRecord.h"
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>


class MsaFileRecord{
    private:
        int start;
        int end;
        std::map<std::string, MsaSingleRecord> msaSingleRecordRecords;
    public:
    //        MsaFileRecord();
        MsaFileRecord(const int & _start, const int & _end);
        const int & getStart();
        const int & getEnd();
        void setStart( const int & _start );
        void setEnd( const int & _end );
        void addMsaSingleRecord(std::string & lineName, MsaSingleRecord & msaSingleRecord);
        std::map<std::string, MsaSingleRecord>& getMsaSingleRecordRecords();
        bool operator<( const MsaFileRecord& msaFileRecord ) const {
            if( this->start < msaFileRecord.start ){
                return true;
            }else if( (this->start == msaFileRecord.start) && (this->end < msaFileRecord.end) ){
                return true;
            }else{
                return false;
            }
        }
        bool operator>( const MsaFileRecord& msaFileRecord ) const {
            if( this->start > msaFileRecord.start ){
                return true;
            }else if( (this->start == msaFileRecord.start) && (this->end > msaFileRecord.end) ){
                return true;
            }else{
                return false;
            }
        }
        bool operator == ( const MsaFileRecord& msaFileRecord ) const {
            return ( (this->start == msaFileRecord.start) && (this->end == msaFileRecord.end) );
        }
        bool operator != ( const MsaFileRecord& msaFileRecord ) const {
            if( (this->start == msaFileRecord.start) && (this->end == msaFileRecord.end) ){
                return false;
            }else{
                return true;
            }
        }
};
void msaFileRead(MsaFileRecord & msaFileRecord, std::string & fileLocation, std::map<std::string, std::string>& sdiFiles);
void msaFileRead(MsaFileRecord & msaFileRecord, std::string & fileLocation, std::set<std::string>& accessionNames);

#endif //ANNOTATIONLIFTOVER_MSAFILERECORD_H
