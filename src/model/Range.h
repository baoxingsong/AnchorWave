//
// Created by song on 8/4/18.
//

#ifndef ANNOTATIONLIFTOVER_RANGE_H
#define ANNOTATIONLIFTOVER_RANGE_H

#include <string>
#include "STRAND.h"
class Range {
private:
    std::string chr;
    std::string name;
    size_t start;
    size_t end;
    STRAND strand;
    size_t index;
public:
    Range(const std::string &chr, const size_t & start, const size_t & end, const STRAND & strand);
    Range();
    Range(const Range & range);
    const std::string &getName() const;
    void setName(const std::string &name);
    const std::string &getChr() const;
    void setChr(const std::string &chr);
    const size_t & getStart() const;
    void setStart(const size_t & start);
    const size_t & getEnd() const;
    void setEnd(const size_t & end);
    const STRAND & getStrand() const;
    void setStrand(const STRAND & strand);

    const size_t & getIndex() const;

    void setIndex(size_t index);
};


#endif //ANNOTATIONLIFTOVER_RANGE_H
