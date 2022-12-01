//
// Created by bs674 on 7/5/21.
//

#ifndef ANCHORWAVE_RANGE_H
#define ANCHORWAVE_RANGE_H

#include <iostream>

class Range {
private:
    int32_t start;
    int32_t end;
public:
    Range(int32_t start, int32_t end);

    int32_t getStart() const;

    void setStart(int32_t start);

    int32_t getEnd() const;

    void setEnd(int32_t end);

    bool operator==(const Range &rag) const {
        if (this->getStart() == rag.getStart() && this->getEnd() == rag.getEnd()) {
            return true;
        }
        return false;
    }

    bool operator!=(const Range &rag) const {
        if (this->getStart() == rag.getStart() && this->getEnd() == rag.getEnd()) {
            return false;
        }
        return true;
    }

    bool operator>(const Range &rag) const {
        return this->getStart() > rag.getStart();
    }

    bool operator<(const Range &rag) const {
        return this->getStart() < rag.getStart();
    }
};

#endif //ANCHORWAVE_RANGE_H
