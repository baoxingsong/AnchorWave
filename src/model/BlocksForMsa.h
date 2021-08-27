//
// Created by bs674 on 6/8/21.
//

#ifndef ANCHORWAVE_BLOCKSFORMSA_H
#define ANCHORWAVE_BLOCKSFORMSA_H
#include <string>
#include <vector>
#include "AlignmentBlock.h"
#include <iostream>
class BlocksForMsa {
private:
    std::vector<AlignmentBlock> alignmentBlocks;
    BlocksForMsa *next = NULL;
    BlocksForMsa *prev = NULL;
    std::string whichRef = "NA";
    std::string ancestralSequence = "";
public:
    std::vector<AlignmentBlock> &getAlignmentBlocks() ;

    const std::string &getWhichRef() const;

    const std::string &getAncestralSequence() const;

    void setAncestralSequence(const std::string &ancestralSequence);

    void setWhichRef(const std::string &whichRef);

    const std::string &getWhileOutPut() const;

    void setWhileOutPut(const std::string &whileOutPut);

    void setAlignmentBlocks(const std::vector<AlignmentBlock> &alignmentBlocks);
    void addAlignmentBlock(const AlignmentBlock & alignmentBlock);
    BlocksForMsa *getNext() const;
    void setNext(BlocksForMsa *next);
    BlocksForMsa *getPrev() const;
    void setPrev(BlocksForMsa *prev);
};


#endif //ANCHORWAVE_BLOCKSFORMSA_H
