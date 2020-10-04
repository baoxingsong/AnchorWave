// =====================================================================================
//
//       Filename:  MsaSingleRecord.h
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

The class is designed for the MSA result.
For each MSA result, is has a groups sequence from different accessions.

This class if designed for the sequence record of each accession


*************************************************************************/



#ifndef ANNOTATIONLIFTOVER_MSASINGLERECORD_H
#define ANNOTATIONLIFTOVER_MSASINGLERECORD_H
#include <string>


class MsaSingleRecord {
    private:
        int start;
        int end;
        std::string accesionName;
        std::string sequence;
    public:
        MsaSingleRecord(const int & _start, const int & _end, const std::string & _accesionName, const std::string & _sequence);
    //        MsaSingleRecord(std::string _accesionName, std::string _sequence);
        MsaSingleRecord();
        const int & getStart() const;
        const int & getEnd() const;
        void setStart( const int & _start);
        void setEnd( const int & _end);
        const std::string & getSequence();
};


#endif //ANNOTATIONLIFTOVER_MSASINGLERECORD_H
