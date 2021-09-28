//
// Created by bs674 on 7/5/21.
//
#include "./Range.h"
#include <vector>
#include <set>
#ifndef ANCHORWAVE_TWOVECTORSOFRANGES_H
#define ANCHORWAVE_TWOVECTORSOFRANGES_H


class TwoVectorsOfRanges {
private:
    std::set<Range> ranges1;
    std::set<Range> ranges2;
    int32_t start;
    int32_t end;
public:
    const std::set<Range> &getRanges1() const;

    void setRanges1(const std::set<Range> &ranges1);

    const std::set<Range> &getRanges2() const;

    void setRanges2(const std::set<Range> &ranges2);

    int32_t getStart() const;

    void setStart(int32_t start);

    int32_t getAnEnd() const;

    void setAnEnd(int32_t anEnd);

    void addRange1( Range & ran);

    void addRange2( Range & ran);
};
#endif //ANCHORWAVE_TWOVECTORSOFRANGES_H
