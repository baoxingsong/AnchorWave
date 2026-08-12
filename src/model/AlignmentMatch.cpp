//
// Created by Baoxing Song on 2019-03-13.
//

#include "AlignmentMatch.h"

#include <mutex>
#include <unordered_set>

namespace {

class AlignmentStringPool {
public:
    const std::string *intern(const std::string &value) {
        std::lock_guard<std::mutex> lock(mutex_);
        return &*values_.insert(value).first;
    }

private:
    std::mutex mutex_;
    std::unordered_set<std::string> values_;
};

const std::string *internAlignmentString(const std::string &value) {
    static AlignmentStringPool pool;
    return pool.intern(value);
}

const std::string *emptyAlignmentString() {
    static const std::string *empty = internAlignmentString(std::string());
    return empty;
}

}  // namespace

AlignmentMatch::AlignmentMatch(const uint32_t &refStartPos, const uint32_t &refEndPos,
                               const uint32_t &queryStartPos, const uint32_t &queryEndPos, const double &score,
                               const STRAND &strand)
        : refChr(emptyAlignmentString()), queryChr(emptyAlignmentString()),
          refStartPos(refStartPos), refEndPos(refEndPos),
          queryStartPos(queryStartPos), queryEndPos(queryEndPos),
          score(score), strand(strand), refId(0), queryId(0),
          referenceGeneName(emptyAlignmentString()),
          queryGeneName(emptyAlignmentString()) {
}

AlignmentMatch::AlignmentMatch()
        : refChr(emptyAlignmentString()), queryChr(emptyAlignmentString()),
          refStartPos(0), refEndPos(0), queryStartPos(0), queryEndPos(0),
          score(0.0), strand(POSITIVE), refId(0), queryId(0),
          referenceGeneName(emptyAlignmentString()),
          queryGeneName(emptyAlignmentString()) {}

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
                               const std::string &queryGeneName) : refChr(internAlignmentString(refChr)), queryChr(internAlignmentString(queryChr)),
                                                                   refStartPos(refStartPos), refEndPos(refEndPos),
                                                                   queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand), refId(refId),
                                                                   queryId(queryId), referenceGeneName(internAlignmentString(referenceGeneName)), queryGeneName(internAlignmentString(queryGeneName)) {
}

AlignmentMatch::AlignmentMatch(const std::string &refChr, const std::string &queryChr,
                               const uint32_t &refStartPos, const uint32_t &refEndPos,
                               const uint32_t &queryStartPos, const uint32_t &queryEndPos, const double &score,
                               const STRAND &strand, const std::string &referenceGeneName,
                               const std::string &queryGeneName) : refChr(internAlignmentString(refChr)), queryChr(internAlignmentString(queryChr)),
                                                                   refStartPos(refStartPos), refEndPos(refEndPos),
                                                                   queryStartPos(queryStartPos), queryEndPos(queryEndPos), score(score), strand(strand),
                                                                   refId(0), queryId(0), referenceGeneName(internAlignmentString(referenceGeneName)), queryGeneName(internAlignmentString(queryGeneName)) {
}

AlignmentMatch::AlignmentMatch(const AlignmentMatch &orthologPair) = default;

const std::string &AlignmentMatch::getRefChr() const {
    return *refChr;
}

void AlignmentMatch::setRefChr(const std::string &refChr) {
    AlignmentMatch::refChr = internAlignmentString(refChr);
}

const std::string &AlignmentMatch::getQueryChr() const {
    return *queryChr;
}

void AlignmentMatch::setQueryChr(const std::string &queryChr) {
    AlignmentMatch::queryChr = internAlignmentString(queryChr);
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
    return *referenceGeneName;
}

void AlignmentMatch::setReferenceGeneName(const std::string &referenceGeneName) {
    AlignmentMatch::referenceGeneName = internAlignmentString(referenceGeneName);
}

const std::string &AlignmentMatch::getQueryGeneName() const {
    return *queryGeneName;
}

void AlignmentMatch::setQueryGeneName(const std::string &queryGeneName) {
    AlignmentMatch::queryGeneName = internAlignmentString(queryGeneName);
}
