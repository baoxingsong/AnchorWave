// =====================================================================================
//
//       Filename:  STRAND.h
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

 Stand information of CDS transcript or gene

*************************************************************************/


#ifndef ANNOTATIONLIFTOVER_STRAND_H
#define ANNOTATIONLIFTOVER_STRAND_H

#include "map"
#include <ostream>

enum STRAND {
    POSITIVE = 0, NEGATIVE
};

std::ostream &operator<<(std::ostream &out, const STRAND value);


#endif //ANNOTATIONLIFTOVER_STRAND_H
