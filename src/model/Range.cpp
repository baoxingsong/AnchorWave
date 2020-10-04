//
// Created by song on 8/4/18.
//

#include "Range.h"

Range::Range(const std::string &chr, const size_t & start, const size_t & end, const STRAND & strand) : chr(chr), start(start), end(end),
                                                                                strand(strand) {}
Range::Range(){
    chr="";
    name ="";
    start=0;
    end=0;
    strand=POSITIVE;
    index=0;
}

Range::Range(const Range & range){
    chr=range.getChr();
    name =range.getName();
    start=range.getStart();
    end=range.getEnd();
    strand=range.getStrand();
    index=range.getIndex();

}
const std::string &Range::getChr() const {
    return chr;
}

void Range::setChr(const std::string &chr) {
    Range::chr = chr;
}

const size_t & Range::getStart() const {
    return start;
}

void Range::setStart(const size_t & start) {
    Range::start = start;
}

const size_t & Range::getEnd() const {
    return end;
}

void Range::setEnd(const size_t & end) {
    Range::end = end;
}

const STRAND & Range::getStrand() const {
    return strand;
}

void Range::setStrand(const STRAND & strand) {
    Range::strand = strand;
}

const size_t & Range::getIndex() const {
    return index;
}

void Range::setIndex(size_t index) {
    Range::index = index;
}

const std::string &Range::getName() const {
    return name;
}

void Range::setName(const std::string &name) {
    Range::name = name;
}
