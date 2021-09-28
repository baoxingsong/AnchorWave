//
// Created by bs674 on 6/8/21.
//

#ifndef ANCHORWAVE_ALIGNMENTBLOCK_H
#define ANCHORWAVE_ALIGNMENTBLOCK_H
#include <string>
#include "STRAND.h"
class AlignmentBlock {
private:
    std::string species;
    std::string chr;
    int32_t start;
    int32_t end;
    std::string alignSeq;
    STRAND strand;
public:
    AlignmentBlock(const std::string &species, const std::string &chr, int32_t start, int32_t end,
                   const std::string &alignSeq, STRAND strand);

    const std::string &getSpecies() const;

    void setSpecies(const std::string &species);

    const std::string &getChr() const;

    void setChr(const std::string &chr);

    int32_t getStart() const;

    void setStart(int32_t start);

    int32_t getEnd() const;

    void setEnd(int32_t rnd);

    const std::string &getAlignSeq() const;

    void setAlignSeq(const std::string &alignSeq);

    STRAND getStrand() const;

    void setStrand(STRAND strand);
};


#endif //ANCHORWAVE_ALIGNMENTBLOCK_H
