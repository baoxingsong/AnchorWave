//
// Created by Baoxing Song on 2019-03-13.
//

#include "AlignmentMatch.h"

AlignmentMatch::AlignmentMatch(const uint32_t &refStartPos, const uint32_t &refEndPos,
                               const uint32_t &queryStartPos, const uint32_t &queryEndPos, const double &score,
                               const STRAND &strand) : refStartPos(refStartPos), refEndPos(refEndPos),
                                                       queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand) {
}

AlignmentMatch::AlignmentMatch() {

}

uint32_t AlignmentMatch::getRefStartPos() const {
    return refStartPos;
}

void AlignmentMatch::setRefStartPos(uint32_t refStartPos) {
    AlignmentMatch::refStartPos = refStartPos;
}

uint32_t AlignmentMatch::getRefEndPos() const {
    return refEndPos;
}

void AlignmentMatch::setRefEndPos(uint32_t refEndPos) {
    AlignmentMatch::refEndPos = refEndPos;
}

uint32_t AlignmentMatch::getQueryStartPos() const {
    return queryStartPos;
}

void AlignmentMatch::setQueryStartPos(uint32_t queryStartPos) {
    AlignmentMatch::queryStartPos = queryStartPos;
}

uint32_t AlignmentMatch::getQueryEndPos() const {
    return queryEndPos;
}

void AlignmentMatch::setQueryEndPos(uint32_t queryEndPos) {
    AlignmentMatch::queryEndPos = queryEndPos;
}

double AlignmentMatch::getScore() const {
    return score;
}

void AlignmentMatch::setScore(double score) {
    AlignmentMatch::score = score;
}

STRAND AlignmentMatch::getStrand() const {
    return strand;
}

void AlignmentMatch::setStrand(STRAND &strand) {
    AlignmentMatch::strand = strand;
}

AlignmentMatch::AlignmentMatch(const std::string &refChr, const std::string &queryChr,
                               const uint32_t &refStartPos, const uint32_t &refEndPos,
                               const uint32_t &queryStartPos, const uint32_t &queryEndPos, const double &score,
                               const STRAND &strand, const int &refId, const int &queryId, const std::string &referenceGeneName,
                               const std::string &queryGeneName) : refChr(refChr), queryChr(queryChr),
                                                                   refStartPos(refStartPos), refEndPos(refEndPos),
                                                                   queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand), refId(refId),
                                                                   queryId(queryId), referenceGeneName(referenceGeneName), queryGeneName(queryGeneName) {
}

AlignmentMatch::AlignmentMatch(const std::string &refChr, const std::string &queryChr,
                               const uint32_t &refStartPos, const uint32_t &refEndPos,
                               const uint32_t &queryStartPos, const uint32_t &queryEndPos, const double &score,
                               const STRAND &strand, const std::string &referenceGeneName,
                               const std::string &queryGeneName) : refChr(refChr), queryChr(queryChr),
                                                                   refStartPos(refStartPos), refEndPos(refEndPos),
                                                                   queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand),
                                                                   referenceGeneName(referenceGeneName), queryGeneName(queryGeneName) {
}

AlignmentMatch::AlignmentMatch(const AlignmentMatch &orthologPair) {
    refChr = orthologPair.getRefChr();
    queryChr = orthologPair.getQueryChr();
    refStartPos = orthologPair.getRefStartPos();
    refEndPos = orthologPair.getRefEndPos();
    queryStartPos = orthologPair.getQueryStartPos();
    queryEndPos = orthologPair.getQueryEndPos();
    score = orthologPair.getScore();
    strand = orthologPair.getStrand();
    refId = orthologPair.getRefId();
    queryId = orthologPair.getQueryId();
    referenceGeneName = orthologPair.getReferenceGeneName();
    queryGeneName = orthologPair.getQueryGeneName();
}

const std::string &AlignmentMatch::getRefChr() const {
    return refChr;
}

void AlignmentMatch::setRefChr(const std::string &refChr) {
    AlignmentMatch::refChr = refChr;
}

const std::string &AlignmentMatch::getQueryChr() const {
    return queryChr;
}

void AlignmentMatch::setQueryChr(const std::string &queryChr) {
    AlignmentMatch::queryChr = queryChr;
}

int AlignmentMatch::getRefId() const {
    return refId;
}

void AlignmentMatch::setRefId(int refId) {
    AlignmentMatch::refId = refId;
}

int AlignmentMatch::getQueryId() const {
    return queryId;
}

void AlignmentMatch::setQueryId(int queryId) {
    AlignmentMatch::queryId = queryId;
}

const std::string &AlignmentMatch::getReferenceGeneName() const {
    return referenceGeneName;
}

void AlignmentMatch::setReferenceGeneName(const std::string &referenceGeneName) {
    AlignmentMatch::referenceGeneName = referenceGeneName;
}

const std::string &AlignmentMatch::getQueryGeneName() const {
    return queryGeneName;
}

void AlignmentMatch::setQueryGeneName(const std::string &queryGeneName) {
    AlignmentMatch::queryGeneName = queryGeneName;
}
