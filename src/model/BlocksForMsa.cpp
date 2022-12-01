//
// Created by bs674 on 6/8/21.
//

#include "BlocksForMsa.h"


std::vector<AlignmentBlock> &BlocksForMsa::getAlignmentBlocks() {
    return alignmentBlocks;
}

void BlocksForMsa::setAlignmentBlocks(const std::vector<AlignmentBlock> &alignmentBlocks) {
    BlocksForMsa::alignmentBlocks = alignmentBlocks;
}

void BlocksForMsa::addAlignmentBlock(const AlignmentBlock &alignmentBlock) {
    this->alignmentBlocks.push_back(alignmentBlock);
}

BlocksForMsa *BlocksForMsa::getNext() const {
    return next;
}

void BlocksForMsa::setNext(BlocksForMsa *next) {
    BlocksForMsa::next = next;
}

BlocksForMsa *BlocksForMsa::getPrev() const {
    return prev;
}

void BlocksForMsa::setPrev(BlocksForMsa *prev) {
    BlocksForMsa::prev = prev;
}

const std::string &BlocksForMsa::getWhichRef() const {
    return whichRef;
}

void BlocksForMsa::setWhichRef(const std::string &whichRef) {
    BlocksForMsa::whichRef = whichRef;
}

const std::string &BlocksForMsa::getAncestralSequence() const {
    return ancestralSequence;
}

void BlocksForMsa::setAncestralSequence(const std::string &ancestralSequence) {
    BlocksForMsa::ancestralSequence = ancestralSequence;
}
