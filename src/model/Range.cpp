//
// Created by bs674 on 7/5/21.
//

#include "Range.h"

Range::Range(int32_t start, int32_t end) : start(start), end(end) {}

int32_t Range::getStart() const {
    return start;
}

void Range::setStart(int32_t start) {
    Range::start = start;
}

int32_t Range::getEnd() const {
    return end;
}

void Range::setEnd(int32_t end) {
    Range::end = end;
}
