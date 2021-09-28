//
// Created by bs674 on 6/8/21.
//

#include "AlignmentBlock.h"

AlignmentBlock::AlignmentBlock(const std::string &species, const std::string &chr, int32_t start, int32_t end,
                               const std::string &alignSeq, STRAND strand) : species(species), chr(chr), start(start),
                                                                             end(end), alignSeq(alignSeq),
                                                                             strand(strand) {}

const std::string &AlignmentBlock::getSpecies() const {
    return species;
}

void AlignmentBlock::setSpecies(const std::string &species) {
    AlignmentBlock::species = species;
}

const std::string &AlignmentBlock::getChr() const {
    return chr;
}

void AlignmentBlock::setChr(const std::string &chr) {
    AlignmentBlock::chr = chr;
}

int32_t AlignmentBlock::getStart() const {
    return start;
}

void AlignmentBlock::setStart(int32_t start) {
    AlignmentBlock::start = start;
}

int32_t AlignmentBlock::getEnd() const {
    return end;
}

void AlignmentBlock::setEnd(int32_t end) {
    AlignmentBlock::end = end;
}

const std::string &AlignmentBlock::getAlignSeq() const {
    return alignSeq;
}

void AlignmentBlock::setAlignSeq(const std::string &alignSeq) {
    AlignmentBlock::alignSeq = alignSeq;
}

STRAND AlignmentBlock::getStrand() const {
    return strand;
}

void AlignmentBlock::setStrand(STRAND strand) {
    AlignmentBlock::strand = strand;
}
