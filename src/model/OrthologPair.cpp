//
// Created by Baoxing Song on 2019-03-13.
//

#include "OrthologPair.h"

OrthologPair2::OrthologPair2(const uint32_t & refStartPos, const uint32_t & refEndPos,
        const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
        const STRAND & strand) : refStartPos(refStartPos), refEndPos(refEndPos),
        queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand) {
//    refMiddlePos = (refStartPos+refEndPos)/2;
//    queryMiddlePos = (queryStartPos+queryEndPos)/2;
}
OrthologPair2::OrthologPair2(){

}

uint32_t OrthologPair2::getRefStartPos() const {
    return refStartPos;
}

void OrthologPair2::setRefStartPos(uint32_t refStartPos) {
    OrthologPair2::refStartPos = refStartPos;
}

uint32_t OrthologPair2::getRefEndPos() const {
    return refEndPos;
}

void OrthologPair2::setRefEndPos(uint32_t refEndPos) {
    OrthologPair2::refEndPos = refEndPos;
}
//
//uint32_t OrthologPair2::getRefMiddlePos() const {
//    return refMiddlePos;
//}
//
//void OrthologPair2::setRefMiddlePos(uint32_t refMiddlePos) {
//    OrthologPair2::refMiddlePos = refMiddlePos;
//}

uint32_t OrthologPair2::getQueryStartPos() const {
    return queryStartPos;
}

void OrthologPair2::setQueryStartPos(uint32_t queryStartPos) {
    OrthologPair2::queryStartPos = queryStartPos;
}

uint32_t OrthologPair2::getQueryEndPos() const {
    return queryEndPos;
}

void OrthologPair2::setQueryEndPos(uint32_t queryEndPos) {
    OrthologPair2::queryEndPos = queryEndPos;
}
//
//uint32_t OrthologPair2::getQueryMiddlePos() const {
//    return queryMiddlePos;
//}
//
//void OrthologPair2::setQueryMiddlePos(uint32_t queryMiddlePos) {
//    OrthologPair2::queryMiddlePos = queryMiddlePos;
//}

double OrthologPair2::getScore() const {
    return score;
}

void OrthologPair2::setScore(double score) {
    OrthologPair2::score = score;
}

STRAND OrthologPair2::getStrand() const {
    return strand;
}

void OrthologPair2::setStrand(STRAND & strand) {
    OrthologPair2::strand = strand;
}


OrthologPair2::OrthologPair2(const std::string &refChr, const std::string &queryChr,
            const uint32_t & refStartPos, const uint32_t & refEndPos,
            const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
            const STRAND & strand, const int & refId, const int & queryId, const std::string & referenceGeneName,
            const std::string & queryGeneName) : refChr(refChr), queryChr(queryChr),
            refStartPos(refStartPos), refEndPos(refEndPos),
            queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand), refId(refId),
            queryId(queryId), referenceGeneName(referenceGeneName), queryGeneName(queryGeneName) {
//    refMiddlePos = (refStartPos+refEndPos)/2;
//    queryMiddlePos = (queryStartPos+queryEndPos)/2;
}

OrthologPair2::OrthologPair2(const std::string &refChr, const std::string &queryChr,
                             const uint32_t & refStartPos, const uint32_t & refEndPos,
                             const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
                             const STRAND & strand, const std::string & referenceGeneName,
                             const std::string & queryGeneName) : refChr(refChr), queryChr(queryChr),
                                                                  refStartPos(refStartPos), refEndPos(refEndPos),
                                                                  queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand),
                                                                  referenceGeneName(referenceGeneName), queryGeneName(queryGeneName) {
//    refMiddlePos = (refStartPos+refEndPos)/2;
//    queryMiddlePos = (queryStartPos+queryEndPos)/2;
}

OrthologPair2::OrthologPair2( const OrthologPair2 & orthologPair )  {
    refChr=orthologPair.getRefChr();
    queryChr=orthologPair.getQueryChr();
    refStartPos=orthologPair.getRefStartPos();
    refEndPos=orthologPair.getRefEndPos();
    queryStartPos=orthologPair.getQueryStartPos();
    queryEndPos=orthologPair.getQueryEndPos();
    score=orthologPair.getScore();
    strand=orthologPair.getStrand();
    refId=orthologPair.getRefId();
    queryId=orthologPair.getQueryId();
    referenceGeneName=orthologPair.getReferenceGeneName();
    queryGeneName=orthologPair.getQueryGeneName();
//    refMiddlePos=orthologPair.getRefMiddlePos();
//    queryMiddlePos=orthologPair.getQueryMiddlePos();
}

const std::string &OrthologPair2::getRefChr() const {
    return refChr;
}

void OrthologPair2::setRefChr(const std::string &refChr) {
    OrthologPair2::refChr = refChr;
}

const std::string &OrthologPair2::getQueryChr() const {
    return queryChr;
}

void OrthologPair2::setQueryChr(const std::string &queryChr) {
    OrthologPair2::queryChr = queryChr;
}

int OrthologPair2::getRefId() const {
    return refId;
}

void OrthologPair2::setRefId(int refId) {
    OrthologPair2::refId = refId;
}

int OrthologPair2::getQueryId() const {
    return queryId;
}

void OrthologPair2::setQueryId(int queryId) {
    OrthologPair2::queryId = queryId;
}

const std::string &OrthologPair2::getReferenceGeneName() const {
    return referenceGeneName;
}
void OrthologPair2::setReferenceGeneName(const std::string &referenceGeneName) {
    OrthologPair2::referenceGeneName = referenceGeneName;
}
const std::string &OrthologPair2::getQueryGeneName() const {
    return queryGeneName;
}
void OrthologPair2::setQueryGeneName(const std::string &queryGeneName) {
    OrthologPair2::queryGeneName = queryGeneName;
}










OrthologPair::OrthologPair(const int & queryIndex, const uint32_t & refStartPos, const uint32_t & refEndPos,
                           const uint32_t & queryStartPos, const uint32_t & queryEndPos, const double & score,
                           const STRAND & strand) : queryIndex(queryIndex), refStartPos(refStartPos), refEndPos(refEndPos),
                                                    queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand) {
    refMiddlePos = (refStartPos+refEndPos)/2;
    queryMiddlePos = (queryStartPos+queryEndPos)/2;
}

uint32_t OrthologPair::getRefStartPos() const {
    return refStartPos;
}

void OrthologPair::setRefStartPos(uint32_t refStartPos) {
    OrthologPair::refStartPos = refStartPos;
}

uint32_t OrthologPair::getRefEndPos() const {
    return refEndPos;
}

void OrthologPair::setRefEndPos(uint32_t refEndPos) {
    OrthologPair::refEndPos = refEndPos;
}

uint32_t OrthologPair::getRefMiddlePos() const {
    return refMiddlePos;
}

void OrthologPair::setRefMiddlePos(uint32_t refMiddlePos) {
    OrthologPair::refMiddlePos = refMiddlePos;
}

uint32_t OrthologPair::getQueryStartPos() const {
    return queryStartPos;
}

void OrthologPair::setQueryStartPos(uint32_t queryStartPos) {
    OrthologPair::queryStartPos = queryStartPos;
}

uint32_t OrthologPair::getQueryEndPos() const {
    return queryEndPos;
}

void OrthologPair::setQueryEndPos(uint32_t queryEndPos) {
    OrthologPair::queryEndPos = queryEndPos;
}

uint32_t OrthologPair::getQueryMiddlePos() const {
    return queryMiddlePos;
}

void OrthologPair::setQueryMiddlePos(uint32_t queryMiddlePos) {
    OrthologPair::queryMiddlePos = queryMiddlePos;
}

double OrthologPair::getScore() const {
    return score;
}

void OrthologPair::setScore(double score) {
    OrthologPair::score = score;
}

STRAND OrthologPair::getStrand() const {
    return strand;
}

void OrthologPair::setStrand(STRAND & strand) {
    OrthologPair::strand = strand;
}

int OrthologPair::getQueryIndex() const {
    return queryIndex;
}

void OrthologPair::setQueryIndex(int queryIndex) {
    OrthologPair::queryIndex = queryIndex;
}