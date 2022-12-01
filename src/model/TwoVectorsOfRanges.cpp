//
// Created by bs674 on 7/5/21.
//

#include "TwoVectorsOfRanges.h"

const std::set<Range> &TwoVectorsOfRanges::getRanges1() const {
    return ranges1;
}

void TwoVectorsOfRanges::setRanges1(const std::set<Range> &ranges1) {
    TwoVectorsOfRanges::ranges1 = ranges1;
}

const std::set<Range> &TwoVectorsOfRanges::getRanges2() const {
    return ranges2;
}

void TwoVectorsOfRanges::setRanges2(const std::set<Range> &ranges2) {
    TwoVectorsOfRanges::ranges2 = ranges2;
}

int32_t TwoVectorsOfRanges::getStart() const {
    return start;
}

void TwoVectorsOfRanges::setStart(int32_t start) {
    TwoVectorsOfRanges::start = start;
}

int32_t TwoVectorsOfRanges::getAnEnd() const {
    return end;
}

void TwoVectorsOfRanges::setAnEnd(int32_t anEnd) {
    end = anEnd;
}

void TwoVectorsOfRanges::addRange1(Range &ran) {
    ranges1.insert(ran);
    if (1 == ranges1.size() && 0 == ranges2.size()) {
        this->start = ran.getStart();
        this->end = ran.getEnd();
        return;
    }
    if (ran.getStart() < this->start) {
        this->start = ran.getStart();
    }
    if (ran.getEnd() > this->end) {
        this->end = ran.getEnd();
    }
}

void TwoVectorsOfRanges::addRange2(Range &ran) {
    ranges2.insert(ran);
    if (0 == ranges1.size() && 1 == ranges2.size()) {
        this->start = ran.getStart();
        this->end = ran.getEnd();
        return;
    }
    if (ran.getStart() < this->start) {
        this->start = ran.getStart();
    }
    if (ran.getEnd() > this->end) {
        this->end = ran.getEnd();
    }
}
