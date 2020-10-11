//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
#define ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include <map>
#include <regex>

void setupAnchorsWithSpliceAlignmentResult( const std::string & gffFilePath, const std::string & samFile, std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,  double & inversion_PENALTY, double & MIN_ALIGNMENT_SCORE, bool & considerInversion);

void setupAnchorsWithSpliceAlignmentResultQuota( const std::string & gffFilePath, const std::string & samFile, std::vector<std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                                 double & INDEL_SCORE, double & GAP_OPEN_PENALTY, double & MIN_ALIGNMENT_SCORE,
                                                 int & MAX_DIST_BETWEEN_MATCHES, int & refMaximumTimes, int & queryMaximumTimes,
                                                 double & calculateIndelDistance);

#endif //ANNOTATIONLIFTOVER_TRANSFERGFFWITHNUCMERRESULT_H
