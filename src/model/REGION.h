// =====================================================================================
//
//       Filename:  REGION.h
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

 If we perform gene structure based MSA.
 It is necessary to know a gene region window is from CDS or INTRON
*************************************************************************/


#ifndef ANNOTATIONLIFTOVER_REGION_H
#define ANNOTATIONLIFTOVER_REGION_H

enum REGION {
    CDS, CDSINTRON
};

#endif //ANNOTATIONLIFTOVER_REGION_H
