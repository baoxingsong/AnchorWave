// =====================================================================================
//
//       Filename:  TwoSeqOfMsaResult.h
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

When recall variants, we compare each non-reference sequence with reference pair-wisely.

 This class is designed to store pair-of sequence




*************************************************************************/



#ifndef ANNOTATIONLIFTOVER_TWOSEQOFMSARESULT_H
#define ANNOTATIONLIFTOVER_TWOSEQOFMSARESULT_H

#include <string>


class TwoSeqOfMsaResult{
    private:
        int refStart;
        int refEnd;
        std::string refSeq;
        int resultStart;
        int resultEnd;
        std::string resultSeq;
    public:
        TwoSeqOfMsaResult( const int & _refStart, const int & _refEnd, const std::string & _refSeq, const int & _resultStart, const int & _resultEnd, const std::string & _resultSeq);
        bool operator<( const TwoSeqOfMsaResult& twoSeqOfMsaResult ) const {
            return this->refStart < twoSeqOfMsaResult.refStart;
        }
        bool operator>( const TwoSeqOfMsaResult& twoSeqOfMsaResult ) const {
            return this->refStart > twoSeqOfMsaResult.refStart;
        }
        const int & getRefStart() const;
        void setRefStart( const int & _refStart);
        const int & getRefEnd() const;
        void setRefEnd( const int & _refEnd);
        const std::string & getRefSeq() const;
        void setRefSeq(const std::string & _refSequence);
        const int & getResultEnd() const;
        void setResultEnd( const int & _resultEnd);
        const int & getResultStart() const;
        void setResultStart(const int & _resultStart);
        const std::string & getResultSeq() const;
        void setResultSeq(const std::string & _resultSequence);
}; // for MSA parse end

#endif //ANNOTATIONLIFTOVER_TWOSEQOFMSARESULT_H
