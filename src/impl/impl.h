// =====================================================================================
//
//       Filename:  impl.h
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

the methods under this folder only rely on code under model and util
 ************************************************************************/


#ifndef ANNOTATIONLIFTOVER_IMPL_H
#define ANNOTATIONLIFTOVER_IMPL_H

#include "../model/model.h"
#include "GetReverseComplementary.h"
#include "getSequencesFromGff.h"
#include "getSubsequence.h"
#include "readFastaFile.h"
#include "readGffFile.h"
#include "TranscriptUpdateInformation.h"
#include "CheckAndUpdateTranscriptsEnds.h"
#include "deNovoGenomeVariantCalling.h"
#include "geneSyntenic.h"
#include "readMafForMsa.h"
#include "WriteFasta.h"
#include "sequencesToKmerDIstanceMatrix.h"
#include "dynmatrix.h"
//#include "findSimilarFragmentsForPairedSequence.h"
//#include "sequenceAlignment.h"

#endif //ANNOTATIONLIFTOVER_IMPL_H
