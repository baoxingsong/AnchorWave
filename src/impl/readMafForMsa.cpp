//
// Created by bs674 on 6/8/21.
//

#include "readMafForMsa.h"


// it assumes the input is a pairwise genome alignment
void readMafForMsa(const std::string &filePath, std::vector<BlocksForMsa> &blocksForMsas) {
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening MAF file " << filePath << std::endl;
        exit(1);
    }
    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;
    while (std::getline(infile, line)) {
        if (line[0] == 'a' && (line[1] == '\t' || line[1] == ' ')) {
            if (std::getline(infile, line)) {
                if (line[0] == 's' && (line[1] == '\t' || line[1] == ' ')) {
                    elems.clear();
                    split(line, delim, elems);
                    STRAND strand = POSITIVE;
                    if (elems[4][0] == '-') {
                        strand = NEGATIVE;
                        std::cerr << "we are have problem to read the MAF file, since the reference sequence is on the NEGATIVE strand" << std::endl;
                    }
                    std::string chr = elems[1];
                    trim(chr);
                    std::string alignSeq = elems[6];
                    trim(alignSeq);
                    BlocksForMsa blocksForMsa;
                    AlignmentBlock alignmentBlock("ref1", chr, stoi(elems[2]) + 1, stoi(elems[2]) + stoi(elems[3]), alignSeq, strand);
                    blocksForMsa.addAlignmentBlock(alignmentBlock);
                    blocksForMsas.push_back(blocksForMsa);

                } else {
                    std::cout << "line 39" << std::endl;
                }
            }
        } else if (line.size() > 1 && line[0] == 's' && (line[1] == '\t' || line[1] == ' ')) {
            elems.clear();
            split(line, delim, elems);
            STRAND strand = POSITIVE;
            if (elems[4][0] == '-') {
                strand = NEGATIVE;
            }
            std::string chr = elems[1];
            trim(chr);
            std::string alignSeq = elems[6];
            trim(alignSeq);
            AlignmentBlock alignmentBlock("ref2", chr, stoi(elems[2]) + 1, stoi(elems[2]) + stoi(elems[3]), alignSeq, strand);
            blocksForMsas[blocksForMsas.size() - 1].addAlignmentBlock(alignmentBlock);
        }
    }
    infile.close();
}


void ancestorInversion(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas) {
    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) { // merge neighborhood nodes together
        for (int32_t j = 0; j < ref1Ref2BlocksForMsas.size(); ++j) {
            if (i != j && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getChr()) == 0 &&
                ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getChr()) == 0 &&
                ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getEnd() < ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getStart()) {
                bool ifLink = true;
                for (int32_t k = 0; k < ref1Ref2BlocksForMsas.size(); ++k) {
                    if (k != i && k != j) {
                        if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getChr()) == 0 && (
                                (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() &&
                                 ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd()) ||
                                (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd()) ||
                                (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() &&
                                 ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd()) ||
                                (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd())

                        )) {
//                            std::cout << "line 73" << std::endl;
                            ifLink = false;
                            break;
                        }
                        if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getChr()) == 0 && (
                                (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() &&
                                 ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd()) ||
                                (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd()) ||
                                (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() &&
                                 ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd()) ||
                                (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd()) ||

                                (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() &&
                                 ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() &&
                                 ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd()) ||
                                (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd()))) {
                            //                          std::cout << "line 87" << std::endl;
                            ifLink = false;
                            break;
                        }
                    }
                }
                if (ifLink) {
                    ref1Ref2BlocksForMsas[i].setNext(&ref1Ref2BlocksForMsas[j]);
                    ref1Ref2BlocksForMsas[j].setPrev(&ref1Ref2BlocksForMsas[i]);
                    break;
                }
            }
        }
    }


    int index = 0;
    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        if (ref1Ref2BlocksForMsas[i].getPrev() == NULL) {
            BlocksForMsa *blocksForMsa = &(ref1Ref2BlocksForMsas[i]);

            while (blocksForMsa->getNext() != NULL) {
                blocksForMsa = blocksForMsa->getNext();
                if (blocksForMsa->getAlignmentBlocks()[1].getStrand() != blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) {
                    int ancestor = 0; // unknown

                    for (int32_t j = 0; j < ref1QueryBlocksForMsas.size(); ++j) {
                        if (ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr()) == 0 &&
                            ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= (*blocksForMsa->getPrev()).getAlignmentBlocks()[0].getEnd() && (*blocksForMsa->getPrev()).getAlignmentBlocks()[0].getEnd() <= ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd() &&
                            ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= blocksForMsa->getAlignmentBlocks()[0].getStart() && blocksForMsa->getAlignmentBlocks()[0].getStart() <= ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd()) {
//                            std::cout << "reference1 inversion is same with the outgroup" << std::endl;
                            ancestor = 1;
                        }
                    }

                    for (int32_t j = 0; j < ref2QueryBlocksForMsas.size(); ++j) {
                        if (ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr()) == 0 &&
                            ((ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= (*blocksForMsa->getPrev()).getAlignmentBlocks()[1].getEnd() &&
                              (*blocksForMsa->getPrev()).getAlignmentBlocks()[1].getEnd() <= ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd()) ||
                             (ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= (*blocksForMsa->getPrev()).getAlignmentBlocks()[1].getStart() &&
                              (*blocksForMsa->getPrev()).getAlignmentBlocks()[1].getStart() <= ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd())) &&
                            ((ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= blocksForMsa->getAlignmentBlocks()[1].getEnd() && blocksForMsa->getAlignmentBlocks()[1].getEnd() <= ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd()) ||
                             (ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= blocksForMsa->getAlignmentBlocks()[1].getStart() && blocksForMsa->getAlignmentBlocks()[1].getStart() <= ref2QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd()))) {
//                            std::cout << "reference2 inversion is same with the outgroup" << std::endl;
                            ancestor += 2;
                        }
                    }
                    if (2 == ancestor) {
//                        std::cout << "line 134: reference2 inversion is same with the outgroup" << std::endl;
                        blocksForMsa->setWhichRef("ref2");
                    } else if (1 == ancestor) {
//                        std::cout << "line 137: reference1 inversion is same with the outgroup" << std::endl;
                        blocksForMsa->setWhichRef("ref1");
                    }
                }
            }
        }
    }
}

void ancestorLink(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas) {
    int index = 0;
    // link using ref1
    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        if (ref1Ref2BlocksForMsas[i].getPrev() == NULL) {
            for (int32_t j = 0; j < ref1Ref2BlocksForMsas.size(); ++j) {
                if (i != j && ref1Ref2BlocksForMsas[j].getNext() == NULL) {
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getChr()) == 0) {
                        if (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() < ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart()) {
                            bool ifHasIntermedia = false;
                            for (int32_t k = 0; k < ref1Ref2BlocksForMsas.size(); ++k) {
                                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getChr()) == 0 && i != k && j != k) {
                                    if ((ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() &&
                                         ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart()) ||
                                        (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd() &&
                                         ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart()) ||
                                        (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() &&
                                         ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd()) ||
                                        (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() &&
                                         ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[0].getEnd())
                                            ) {
                                        ifHasIntermedia = true;
                                        break;
                                    }
                                }
                            }
                            // do not care about the ref2 here. If the ref1 and outgroup support they could be linked together, there maybe a rearrangement in the ref2
                            if (!ifHasIntermedia) {
                                for (int32_t l = 0; l < ref1QueryBlocksForMsas.size(); ++l) { // we do not have strand problem for ref1
                                    if (ref1QueryBlocksForMsas[l].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr()) == 0 &&
                                        ref1QueryBlocksForMsas[l].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[0].getEnd() &&
                                        ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1QueryBlocksForMsas[l].getAlignmentBlocks()[0].getEnd()) {
                                        ref1Ref2BlocksForMsas[i].setPrev(&(ref1Ref2BlocksForMsas[j]));
                                        ref1Ref2BlocksForMsas[i].setWhichRef("ref1Links");
                                        ref1Ref2BlocksForMsas[j].setNext(&(ref1Ref2BlocksForMsas[i]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // link using ref2
    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        if (ref1Ref2BlocksForMsas[i].getPrev() == NULL) {
            for (int32_t j = 0; j < ref1Ref2BlocksForMsas.size(); ++j) {
                if (i != j && ref1Ref2BlocksForMsas[j].getNext() == NULL && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStrand() == ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStrand()) { // it is very hard to deal with strand problem here
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getChr()) == 0) {
//                        if ( ( (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStrand() == POSITIVE && ref1Ref2BlocksForMsas[j].getPrev()== NULL) ||
//                            (ref1Ref2BlocksForMsas[j].getPrev()!= NULL && ref1Ref2BlocksForMsas[j].getPrev()->getAlignmentBlocks()[1].getEnd() < ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() ) ) &&
//                            ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() < ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart()) {
                        if (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStrand() == POSITIVE && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() < ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart()) {
                            bool ifHasIntermedia = false;
                            for (int32_t k = 0; k < ref1Ref2BlocksForMsas.size(); ++k) {
                                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getChr()) == 0 && i != k && j != k) {
                                    if ((ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() &&
                                         ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart()) ||
                                        (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() &&
                                         ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart()) ||
                                        (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() &&
                                         ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd()) ||
                                        (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() &&
                                         ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() <= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd())
                                            ) {
                                        ifHasIntermedia = true;
                                    }
                                }
                            }
                            if (!ifHasIntermedia) {
                                for (int32_t l = 0; l < ref2QueryBlocksForMsas.size(); ++l) {
                                    if (ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr()) == 0 &&
                                        ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() &&
                                        ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() <= ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getEnd()) {
                                        ref1Ref2BlocksForMsas[i].setPrev(&(ref1Ref2BlocksForMsas[j]));
                                        ref1Ref2BlocksForMsas[j].setNext(&(ref1Ref2BlocksForMsas[i]));
                                        ref1Ref2BlocksForMsas[i].setWhichRef("ref2Links");
                                    }
                                }
                            }
                        } else if (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStrand() == NEGATIVE && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) {
                            bool ifHasIntermedia = false;
                            for (int32_t k = 0; k < ref1Ref2BlocksForMsas.size(); ++k) {
                                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getChr()) == 0 && i != k && j != k) {
                                    if ((ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() >= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() &&
                                         ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() >= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                        (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() >= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() &&
                                         ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() >= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                        (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() >= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() &&
                                         ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() >= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart()) ||
                                        (ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() >= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() &&
                                         ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() >= ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart())
                                            ) {
                                        ifHasIntermedia = true;
                                    }
                                }
                            }
                            if (!ifHasIntermedia) {
                                for (int32_t l = 0; l < ref2QueryBlocksForMsas.size(); ++l) {
                                    if (ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr()) == 0 &&
                                        ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() &&
                                        ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() <= ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getEnd()) {
                                        ref1Ref2BlocksForMsas[i].setPrev(&(ref1Ref2BlocksForMsas[j]));
                                        ref1Ref2BlocksForMsas[j].setNext(&(ref1Ref2BlocksForMsas[i]));
                                        ref1Ref2BlocksForMsas[i].setWhichRef("ref2Links");
                                    }
                                }
                            }
                        }/*else if ( ( (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStrand() == NEGATIVE && ref1Ref2BlocksForMsas[j].getPrev()== NULL) || (ref1Ref2BlocksForMsas[j].getPrev()!= NULL && ref1Ref2BlocksForMsas[j].getPrev()->getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() )  )
                                   && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) {
                            bool ifHasIntermedia = false;
                            for (int32_t k = 0; k < ref1Ref2BlocksForMsas.size(); ++k) {
                                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getChr()) == 0 && i!=k && j!=k) {
                                    if( (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                        ( ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                        ( ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() > ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() >  ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() ) ||
                                        ( ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() >  ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() )
                                            ) {
                                        ifHasIntermedia = true;
                                    }
                                }
                            }
                            if (!ifHasIntermedia) {
                                for (int32_t l = 0; l < ref2QueryBlocksForMsas.size(); ++l) {
                                    if (ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr()) == 0 &&
                                        ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() &&
                                        ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() <= ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getEnd()) {
                                        ref1Ref2BlocksForMsas[i].setPrev(&(ref1Ref2BlocksForMsas[j]));
                                        ref1Ref2BlocksForMsas[j].setNext(&(ref1Ref2BlocksForMsas[i]));
                                        ref1Ref2BlocksForMsas[i].setWhichRef("ref2Links");
                                    }
                                }
                            }
                        }else if ( ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStrand() == NEGATIVE && ref1Ref2BlocksForMsas[i].getNext() == NULL  && ref1Ref2BlocksForMsas[j].getPrev()== NULL && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() < ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() ){
                            bool ifHasIntermedia = false;
                            for (int32_t k = 0; k < ref1Ref2BlocksForMsas.size(); ++k) {
                                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr().compare(ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getChr()) == 0 && i!=k && j!=k) {
                                    if( (ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                        ( ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() > ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd()) ||
                                        ( ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() > ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() && ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getStart() >  ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() ) ||
                                        ( ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getEnd() > ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getEnd() >  ref1Ref2BlocksForMsas[k].getAlignmentBlocks()[1].getStart() )
                                            ) {
                                        ifHasIntermedia = true;
                                    }
                                }
                            }
                            if (!ifHasIntermedia) {
                                for (int32_t l = 0; l < ref2QueryBlocksForMsas.size(); ++l) {
                                    if (ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getChr()) == 0 &&
                                        ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getStart() <=ref1Ref2BlocksForMsas[j].getAlignmentBlocks()[1].getEnd() && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getStart() <=ref2QueryBlocksForMsas[l].getAlignmentBlocks()[0].getEnd()) {
                                        ref1Ref2BlocksForMsas[i].setPrev(&(ref1Ref2BlocksForMsas[j]));
                                        ref1Ref2BlocksForMsas[j].setNext(&(ref1Ref2BlocksForMsas[i]));
                                        ref1Ref2BlocksForMsas[i].setWhichRef("ref2Links");
                                    }
                                }
                            }
                        }*/
                    }
                }
            }
        }
    }
}

void generateMsa(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas,
                 std::map<std::string, std::string> &reference1Genome, std::map<std::string, std::string> &reference2Genome) {
    int32_t matchingScore = 0;
    int32_t mismatchingPenalty = -4;
    int32_t open_gap_penalty1 = -4;
    int32_t extend_gap_penalty1 = -2;
    int32_t open_gap_penalty2 = -80;
    int32_t extend_gap_penalty2 = -1;
    int64_t slidingWindowSize = 500;
    int32_t wfaSize = 15000;
    int32_t min_wavefront_length = 10000;
    int32_t max_distance_threshold = 10000;

    int miniInsertionSize = 15;
    int maxDistance = 25;

    generateMsa(ref1Ref2BlocksForMsas, ref1QueryBlocksForMsas, ref2QueryBlocksForMsas, reference1Genome,
                reference2Genome, matchingScore, mismatchingPenalty, open_gap_penalty1, extend_gap_penalty1,
                open_gap_penalty2, extend_gap_penalty2, slidingWindowSize, wfaSize, min_wavefront_length, max_distance_threshold,
                miniInsertionSize, maxDistance);


}


void generateMsa(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas,
                 std::map<std::string, std::string> &reference1Genome, std::map<std::string, std::string> &reference2Genome,
                 int32_t &matchingScore, int32_t &mismatchingPenalty, int32_t &open_gap_penalty1, int32_t &extend_gap_penalty1,
                 int32_t &open_gap_penalty2, int32_t &extend_gap_penalty2, int64_t &slidingWindowSize, int32_t &wfaSize,
                 int32_t &min_wavefront_length, int32_t &max_distance_threshold, int &miniInsertionSize, int &maxDistance) {

    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        STRAND strand = POSITIVE;
        AlignmentBlock alignmentBlock("query", "chr", 0, 1, std::string(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().size(), '*'), strand); // here we use * as unaligned
        ref1Ref2BlocksForMsas[i].addAlignmentBlock(alignmentBlock);
    }

    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        for (int32_t j = 0; j < ref1QueryBlocksForMsas.size(); ++j) {
            if (ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getChr().compare(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getChr()) == 0 && (
                    (ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getEnd()) ||
                    (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart() <= ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() && ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart() <= ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getEnd())
            )) {
                int32_t ref1Ref2Position = ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getStart();
                int32_t ref1QueryPosition = ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getStart();
                int32_t ref1Ref2Index = 0;
                int32_t ref1QueryIndex = 0;
                std::string ref1;
                std::string ref2;
                std::string query;
                ref1.reserve(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().size());
                ref2.reserve(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().size());
                query.reserve(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().size());
                while (ref1Ref2Position < ref1QueryPosition) {
                    ref1 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index];
                    ref2 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[ref1Ref2Index];
                    query += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[ref1Ref2Index];
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index] != '-') {
                        ++ref1Ref2Position;
                    }
                    ++ref1Ref2Index;
                }
//                    std::cout << "line 236\t" << ref1Ref2Position << "\t" << ref1QueryPosition << std::endl;
                while (ref1QueryPosition < ref1Ref2Position) {
                    if (ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getAlignSeq()[ref1QueryIndex] != '-') {
                        ++ref1QueryPosition;
                    }
                    ++ref1QueryIndex;
                }
//                    std::cout << "line 243\t" << ref1Ref2Position << "\t" << ref1QueryPosition << std::endl;
                while (ref1Ref2Index < ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().size() &&
                       ref1QueryIndex < ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getAlignSeq().length()) {
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index] ==
                        ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getAlignSeq()[ref1QueryIndex]) {
                        ref1 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index];
                        ref2 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[ref1Ref2Index];
                        query += ref1QueryBlocksForMsas[j].getAlignmentBlocks()[1].getAlignSeq()[ref1QueryIndex];
                        ++ref1Ref2Index;
                        ++ref1QueryIndex;
                    } else if (
                            ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index] == '-' &&
                            ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getAlignSeq()[ref1QueryIndex] != '-') {
                        ref1 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index];
                        ref2 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[ref1Ref2Index];
                        query += '-';
                        ++ref1Ref2Index;
                    } else if (
                            ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[ref1Ref2Index] != '-' &&
                            ref1QueryBlocksForMsas[j].getAlignmentBlocks()[0].getAlignSeq()[ref1QueryIndex] == '-') {
                        ref1 += '-';
                        ref2 += '-';
                        query += ref1QueryBlocksForMsas[j].getAlignmentBlocks()[1].getAlignSeq()[ref1QueryIndex];
                        ++ref1QueryIndex;
                    } else {
                        std::cerr << "todo line 419" << std::endl;
                    }
                }
                ref1.append(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().substr(ref1Ref2Index));
                ref2.append(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq().substr(ref1Ref2Index));
                query.append(ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().substr(ref1Ref2Index));

                ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].setAlignSeq(ref1);
                ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].setAlignSeq(ref2);
                ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].setAlignSeq(query);
            }
        }
    }

//    std::cout << "line 374\t" << ref1Ref2BlocksForMsas.size() << std::endl;

    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {

        int32_t ref2Position = 0;
        int32_t ref2Length = 0;
        int32_t queryPosition = 0;
        int32_t queryLength = 0;

        std::vector<Range> rangesRef2;
        std::vector<Range> rangesQuery;

        for (int32_t pi = 0; pi < ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().length(); ++pi) {
            if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[pi] == '-') {
                if (ref2Length > 0 && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] != '-' && ref2Position + ref2Length == pi) {
                    ++ref2Length;
                } else {
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] == '-') {
                        if (ref2Length > miniInsertionSize) {
                            Range range(ref2Position + 1, ref2Position + ref2Length);
                            rangesRef2.push_back(range);
                        }
                    }
                    ref2Length = 0;
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] != '-') {
                        ref2Position = pi;
                        ref2Length = 1;
                    }
                }
            } else {
                if (ref2Length > miniInsertionSize) {
                    Range range(ref2Position + 1, ref2Position + ref2Length);
                    rangesRef2.push_back(range);
                }
                ref2Position = pi;
                ref2Length = 0;
            }

            if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[pi] == '-') {
                if (queryLength > 0 && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] != '-' && queryPosition + queryLength == pi) {
                    ++queryLength;
                } else {
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] == '-') {
                        if (queryLength > miniInsertionSize) {
                            Range range(queryPosition + 1, queryPosition + queryLength);
                            rangesQuery.push_back(range);
                        }
                    }
                    queryLength = 0;
                    if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] != '-') {
                        queryPosition = pi;
                        queryLength = 1;
                    }
                }
            } else {
                if (queryLength > miniInsertionSize) {
                    Range range(queryPosition + 1, queryPosition + queryLength);
                    rangesQuery.push_back(range);
                }
                queryPosition = pi;
                queryLength = 0;
            }
        }
//
//        std::cout << "line 400\t" << rangesRef2.size() <<  std::endl;
//        for( int32_t pi = 0; pi < rangesRef2.size(); ++pi ){
//            std::cout << rangesRef2[pi].getStart() << "\t" << rangesRef2[pi].getEnd() << std::endl;
//        }
//        std::cout << "line 404\t" << rangesQuery.size() << std::endl;
//        for( int32_t pi = 0; pi < rangesQuery.size(); ++pi ){
//            std::cout << rangesQuery[pi].getStart() << "\t" << rangesQuery[pi].getEnd() << std::endl;
//        }
        std::vector<TwoVectorsOfRanges> twoVectorsOfRangesVector;
        for (int32_t pi = 0; pi < rangesRef2.size(); ++pi) {
            TwoVectorsOfRanges twoVectorsOfRanges;
            twoVectorsOfRanges.addRange1(rangesRef2[pi]);
            for (int32_t pj = 0; pj < rangesQuery.size(); ++pj) {
                if (rangesQuery[pj].getEnd() < rangesRef2[pi].getStart()) {
                    if (rangesRef2[pi].getStart() - rangesQuery[pj].getEnd() < maxDistance) {
                        twoVectorsOfRanges.addRange2(rangesQuery[pj]);
                    }
                } else if (rangesQuery[pj].getStart() > rangesRef2[pi].getEnd()) {
                    if (rangesQuery[pj].getStart() - rangesRef2[pi].getEnd() < maxDistance) {
                        twoVectorsOfRanges.addRange2(rangesQuery[pj]);
                    }
                } else {
                    twoVectorsOfRanges.addRange2(rangesQuery[pj]);
                }
            }
            if (twoVectorsOfRanges.getRanges2().size() > 0) {
                twoVectorsOfRangesVector.push_back(twoVectorsOfRanges);
            }
        }
//        std::cout << "line 422\t" << twoVectorsOfRangesVector.size()  << std::endl;
//        for( int32_t pi = 0; pi < twoVectorsOfRangesVector.size(); ++pi ){
//            std::cout << "line 424:" << twoVectorsOfRangesVector[pi].getStart() << "\t" << twoVectorsOfRangesVector[pi].getAnEnd() << std::endl;
//        }
        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<int32_t> toremoves;
            for (int32_t pi = 0; pi < twoVectorsOfRangesVector.size(); ++pi) {
                for (int32_t pj = pi + 1; pj < twoVectorsOfRangesVector.size(); ++pj) {
                    bool overlap = false;
                    for (Range r: twoVectorsOfRangesVector[pj].getRanges2()) {
                        if (twoVectorsOfRangesVector[pi].getRanges2().find(r) != twoVectorsOfRangesVector[pi].getRanges2().end()) {
                            overlap = true;
                            break;
                        }
                    }
                    if (overlap) {
                        for (Range r: twoVectorsOfRangesVector[pj].getRanges2()) {
                            twoVectorsOfRangesVector[pi].addRange2(r);
                        }
                        for (Range r: twoVectorsOfRangesVector[pj].getRanges1()) {
                            twoVectorsOfRangesVector[pi].addRange1(r);
                        }
                        toremoves.push_back(pj);
                        changed = true;
                        pi = twoVectorsOfRangesVector.size();
                        break;
                    }
                }
            }
            std::reverse(toremoves.begin(), toremoves.end());
            for (int32_t tr: toremoves) {
                twoVectorsOfRangesVector.erase(twoVectorsOfRangesVector.begin() + tr);
            }
        }
        //std::cout << "line 454\t" << twoVectorsOfRangesVector.size()  << std::endl;
        std::sort(twoVectorsOfRangesVector.begin(), twoVectorsOfRangesVector.end(), [](TwoVectorsOfRanges a, TwoVectorsOfRanges b) {
            return a.getStart() < b.getStart();
        });

        std::stack<char> SQ;
        std::stack<char> SD;


        Scorei m(matchingScore, mismatchingPenalty);
        std::map<std::string, std::string> parameters;
        std::string alin0 = "";
        std::string alin1 = "";
        std::string alin2 = "";
        int lastEnd = 0;
        for (int32_t pi = 0; pi < twoVectorsOfRangesVector.size(); ++pi) {

            std::string _alignment_query;
            std::string _alignment_ref2;
            std::string _alignment_ref1;

            alin0 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().substr(lastEnd, twoVectorsOfRangesVector[pi].getStart() - lastEnd - 1);
            alin1 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq().substr(lastEnd, twoVectorsOfRangesVector[pi].getStart() - lastEnd - 1);
            alin2 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().substr(lastEnd, twoVectorsOfRangesVector[pi].getStart() - lastEnd - 1);
            std::string _dna_ref1 = ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().substr(twoVectorsOfRangesVector[pi].getStart() - 1, twoVectorsOfRangesVector[pi].getAnEnd() - twoVectorsOfRangesVector[pi].getStart() + 1);
            std::string _dna_ref2 = ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq().substr(twoVectorsOfRangesVector[pi].getStart() - 1, twoVectorsOfRangesVector[pi].getAnEnd() - twoVectorsOfRangesVector[pi].getStart() + 1);
            std::string _dna_query = ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().substr(twoVectorsOfRangesVector[pi].getStart() - 1, twoVectorsOfRangesVector[pi].getAnEnd() - twoVectorsOfRangesVector[pi].getStart() + 1);

//            std::cout << "line 514" << std::endl << _dna_ref1 << std::endl;
//            std::cout << _dna_ref2 << std::endl;
//            std::cout << _dna_query << std::endl;
            _dna_ref1.erase(remove(_dna_ref1.begin(), _dna_ref1.end(), '-'), _dna_ref1.end());
            _dna_ref2.erase(remove(_dna_ref2.begin(), _dna_ref2.end(), '-'), _dna_ref2.end());
            _dna_query.erase(remove(_dna_query.begin(), _dna_query.end(), '-'), _dna_query.end());
            _dna_query.erase(remove(_dna_query.begin(), _dna_query.end(), '*'), _dna_query.end());
            if (_dna_query.length() > 0) {
                int64_t thiScore = alignSlidingWindow(_dna_query, _dna_ref2, _alignment_query, _alignment_ref2,
                                                      slidingWindowSize, wfaSize, matchingScore, mismatchingPenalty,
                                                      open_gap_penalty1, extend_gap_penalty1, open_gap_penalty2,
                                                      extend_gap_penalty2, min_wavefront_length, max_distance_threshold, m, parameters);

                thiScore = alignSlidingWindow(_alignment_ref2, _alignment_query, _dna_ref1, _alignment_ref1, slidingWindowSize, matchingScore, mismatchingPenalty,
                                              open_gap_penalty1, extend_gap_penalty1, open_gap_penalty2,
                                              extend_gap_penalty2);
//                std::cout<< "line 530" << std::endl << _alignment_ref1 << std::endl;
//                std::cout << _alignment_ref2 << std::endl;
//                std::cout << _alignment_query << std::endl << std::endl;
                alin0 += _alignment_ref1; // to do
                alin1 += _alignment_ref2;
                alin2 += _alignment_query;
            } else {
                alin0 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().substr(twoVectorsOfRangesVector[pi].getStart() - 1, twoVectorsOfRangesVector[pi].getAnEnd() - twoVectorsOfRangesVector[pi].getStart() + 1);
                alin1 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq().substr(twoVectorsOfRangesVector[pi].getStart() - 1, twoVectorsOfRangesVector[pi].getAnEnd() - twoVectorsOfRangesVector[pi].getStart() + 1);
                alin2 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().substr(twoVectorsOfRangesVector[pi].getStart() - 1, twoVectorsOfRangesVector[pi].getAnEnd() - twoVectorsOfRangesVector[pi].getStart() + 1);
            }
//            std::cout <<  alin0 << "\t" << ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().substr(lastEnd, twoVectorsOfRangesVector[pi].getAnEnd()-lastEnd+1) << std::endl;
            lastEnd = twoVectorsOfRangesVector[pi].getAnEnd();
//            std::cout << "line 543\t" << alin0.size() << "\t" << lastEnd << std::endl;
            //pi = twoVectorsOfRangesVector.size();
        }
//        std::cout << "line 545\t" << ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().size() << "\t" << alin0.size() << std::endl;
        alin0 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().substr(lastEnd, ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().size() - lastEnd);
        alin1 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq().substr(lastEnd, ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().size() - lastEnd);
        alin2 += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq().substr(lastEnd, ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().size() - lastEnd);
        std::cout << "line 549\t" << ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().size() << "\t" << alin0.size() << std::endl;

        ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].setAlignSeq(alin0);
        ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].setAlignSeq(alin1);
        ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].setAlignSeq(alin2);
    }

    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        std::string seq;
        for (int32_t pi = 0; pi < ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq().length(); ++pi) {
            if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] == ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] ||
                (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] == '-' && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] == '*')) {
                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] != '-') {
                    seq += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi];
                }
            } else if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[pi] != '-') {
                seq += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[pi];
            } else if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[0].getAlignSeq()[pi] == '-') {
                if (ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi] != '-' && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] != '-' && ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[2].getAlignSeq()[pi] != '*') {
                    seq += ref1Ref2BlocksForMsas[i].getAlignmentBlocks()[1].getAlignSeq()[pi];
                }// else that is a - and we do not need to do anything here.
            }
        }
        ref1Ref2BlocksForMsas[i].setAncestralSequence(seq);
    }
}


void outputAncestral(std::vector<BlocksForMsa> &ref1Ref2BlocksForMsas, std::vector<BlocksForMsa> &ref1QueryBlocksForMsas, std::vector<BlocksForMsa> &ref2QueryBlocksForMsas,
                     std::map<std::string, std::string> &reference1Genome, std::map<std::string, std::string> &reference2Genome, std::string &outPutFile) {
    int index = 0;

    std::ofstream oCfile;
    oCfile.open(outPutFile);
    std::map<std::string, std::string> seqs;
    index = 1;
    for (int32_t i = 0; i < ref1Ref2BlocksForMsas.size(); ++i) {
        if (ref1Ref2BlocksForMsas[i].getPrev() == NULL) {
            std::string seq = "";
            BlocksForMsa *blocksForMsa = &(ref1Ref2BlocksForMsas[i]);
            bool usingRef1Strand = true;
            if (blocksForMsa->getNext() != NULL && blocksForMsa->getAlignmentBlocks()[1].getStrand() == POSITIVE && /*blocksForMsa->getNext()->getAlignmentBlocks()[1].getStrand() == NEGATIVE &&*/
                blocksForMsa->getAlignmentBlocks()[1].getStart() > blocksForMsa->getNext()->getAlignmentBlocks()[1].getEnd() && blocksForMsa->getNext()->getWhichRef().compare("ref2") == 0) {
                seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                usingRef1Strand = false;
                std::cout << "line 654" << std::endl;
            } else if (blocksForMsa->getNext() != NULL && blocksForMsa->getAlignmentBlocks()[1].getStrand() == NEGATIVE && /* blocksForMsa->getNext()->getAlignmentBlocks()[1].getStrand() == POSITIVE &&*/
                       blocksForMsa->getAlignmentBlocks()[1].getEnd() < blocksForMsa->getNext()->getAlignmentBlocks()[1].getStart() && blocksForMsa->getNext()->getWhichRef().compare("ref2") == 0) {//tested
                seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                usingRef1Strand = false;
                std::cout << "line 658" << std::endl;
            } else if (blocksForMsa->getNext() != NULL && blocksForMsa->getAlignmentBlocks()[1].getStrand() == NEGATIVE &&
                       blocksForMsa->getAlignmentBlocks()[1].getStart() > blocksForMsa->getNext()->getAlignmentBlocks()[1].getEnd() && blocksForMsa->getNext()->getWhichRef().compare("ref2Links") == 0) {
                seq += blocksForMsa->getAncestralSequence();
                blocksForMsa = blocksForMsa->getNext();
                seq += getReverseComplementary(getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() - 1, blocksForMsa->getAlignmentBlocks()[1].getEnd() + 1));
                seq += blocksForMsa->getAncestralSequence();
                usingRef1Strand = true;
//                std::cout << "line 730, this is not expected to happen. if you see this line, please report bug" << std::endl;
            } else { // tested
                std::cout << "line 668" << std::endl;
                seq += blocksForMsa->getAncestralSequence();
            }

            /*
            std::cout << blocksForMsa->getAlignmentBlocks()[0].getChr() << "\t"
                      << blocksForMsa->getAlignmentBlocks()[0].getStart() << "\t"
                      << blocksForMsa->getAlignmentBlocks()[0].getEnd() << "\t"
                      << blocksForMsa->getAlignmentBlocks()[1].getChr() << "\t"
                      << blocksForMsa->getAlignmentBlocks()[1].getStart() << "\t"
                      << blocksForMsa->getAlignmentBlocks()[1].getEnd() << "\t"
                      <<  blocksForMsa->getAlignmentBlocks()[1].getStrand() << "\t"
                      << blocksForMsa->getWhichRef() << std::endl;
*/
            while (blocksForMsa->getNext() != NULL) {
                blocksForMsa = blocksForMsa->getNext();
//                std::cout << "line 650\t" << usingRef1Strand << std::endl;
                if ((blocksForMsa->getWhichRef().compare("ref1") == 0 || blocksForMsa->getWhichRef().compare("NA") == 0) && usingRef1Strand && blocksForMsa->getAlignmentBlocks()[1].getStrand() == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) {
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {// tested
                        std::cout << "line 712" << std::endl;
                        seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                    } else {
                        std::cout << "line 462 unknown" << std::endl;
                    }
                    seq += blocksForMsa->getAncestralSequence();
                } else if ((blocksForMsa->getWhichRef().compare("ref1") == 0 || blocksForMsa->getWhichRef().compare("NA") == 0) && !usingRef1Strand && NEGATIVE == blocksForMsa->getAlignmentBlocks()[1].getStrand() &&
                           blocksForMsa->getAlignmentBlocks()[1].getStrand() == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) {
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[1].getStart()) { // tested
                        //                      std::cout << "line 722" << std::endl;
                        seq += getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1);
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[1].getStart()) {

                    } else {
                        std::cout << "line 471 unknown\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() << "\t" << blocksForMsa->getAlignmentBlocks()[1].getStart() << std::endl;
                    }
                    seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                } else if (blocksForMsa->getWhichRef().compare("ref1Links") == 0 && usingRef1Strand) {
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) { //tested
                        std::cout << "line 732" << std::endl;
                        seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                    } else {
                        std::cout << "line 480 unknown" << std::endl;
                    }
                    seq += blocksForMsa->getAncestralSequence();
                } else if (blocksForMsa->getWhichRef().compare("ref1Links") == 0 && !usingRef1Strand) {
                    std::cout << "line 741" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {
                        seq += getReverseComplementary(getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1));
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                    } else {
                        std::cout << "line 489 unknown\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() << "\t" << blocksForMsa->getAlignmentBlocks()[0].getStart() << std::endl;
                    }
                    seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                } else if (blocksForMsa->getWhichRef().compare("ref2Links") == 0 && blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() == NEGATIVE && usingRef1Strand) {
                    std::cout << "line 751" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() > blocksForMsa->getAlignmentBlocks()[1].getEnd() + 1) {
                        seq += getReverseComplementary(getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1));
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() == blocksForMsa->getAlignmentBlocks()[1].getEnd() + 1) {

                    } else {
                        std::cout << "line 757 unknown\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() << "\t" << blocksForMsa->getAlignmentBlocks()[0].getStart() << std::endl;
                    }
                    seq += (blocksForMsa->getAncestralSequence());
                } else if (blocksForMsa->getWhichRef().compare("ref2Links") == 0 && blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() == NEGATIVE && !usingRef1Strand) {
                    std::cout << "line 761" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() > blocksForMsa->getAlignmentBlocks()[1].getEnd() + 1) {
                        seq += (getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1));
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[1].getStart()) {

                    } else {
                        std::cout << "line 765 unknown\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() << "\t" << blocksForMsa->getAlignmentBlocks()[0].getStart() << std::endl;
                    }
                    seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                } else if (blocksForMsa->getWhichRef().compare("ref2Links") == 0 && blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() == POSITIVE) {
                    //std::cout << "line 771" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                        seq += getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1);
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[1].getStart()) {

                    } else {
                        std::cout << "line 774 unknown\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() << "\t" << blocksForMsa->getAlignmentBlocks()[0].getStart() << std::endl;
                    }
                    seq += blocksForMsa->getAncestralSequence();
                } else if (blocksForMsa->getWhichRef().compare("ref2") == 0 && blocksForMsa->getAlignmentBlocks()[1].getStrand() != blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() && POSITIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() &&
                           usingRef1Strand) {
                    //std::cout << "line 781" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[1].getStart()) { // tested
                        seq += getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1);
                        seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                        usingRef1Strand = false;
                        //std::cout << "line 694" << std::endl;
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                        seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                        usingRef1Strand = false;
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() > blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                        //std::cout << "line 736" << std::endl;
                        seq += getReverseComplementary(getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1));
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                    } else {
                        std::cout << "line 467 unknown" << std::endl;
                    }
                } else if (blocksForMsa->getWhichRef().compare("ref2") == 0 && blocksForMsa->getAlignmentBlocks()[1].getStrand() != blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() && NEGATIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() &&
                           !usingRef1Strand) {
                    std::cout << "line 799" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() > blocksForMsa->getAlignmentBlocks()[1].getEnd()) {
                        seq += getReverseComplementary(getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1));
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                        std::cout << "line 710" << std::endl;
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[1].getStart()) { // tested
                        seq += getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1);
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true; // tested
                        std::cout << "line 715" << std::endl;
                    }/*else if( blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd()+1 == blocksForMsa->getAlignmentBlocks()[1].getStart() ){
                        //seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                    }*/else {
                        std::cout << "line 487 unknown" << std::endl;
                    }
                } else if (blocksForMsa->getWhichRef().compare("ref2") == 0 && blocksForMsa->getAlignmentBlocks()[1].getStrand() == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() && usingRef1Strand) {
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {
                        std::cout << "line 761" << std::endl;
                        seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                    } else {
                        std::cout << "line 487 unknown" << std::endl;
                    }
                } else if (blocksForMsa->getWhichRef().compare("ref2") == 0 && NEGATIVE == blocksForMsa->getAlignmentBlocks()[1].getStrand() && blocksForMsa->getAlignmentBlocks()[1].getStrand() != blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() && !usingRef1Strand &&
                           blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() > blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                    std::cout << "line 769" << std::endl;
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() > blocksForMsa->getAlignmentBlocks()[1].getStart() + 1) {
                        seq += getReverseComplementary(getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getAlignmentBlocks()[1].getStart() + 1, blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() - 1));
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                    } else {
                        std::cout << "line 487 unknown" << std::endl;
                    }
                } else if ((blocksForMsa->getWhichRef().compare("ref1") == 0 || blocksForMsa->getWhichRef().compare("NA") == 0) && !usingRef1Strand && NEGATIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand() &&
                           NEGATIVE == blocksForMsa->getAlignmentBlocks()[1].getStrand() && blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() < blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                        std::cout << "line 737" << std::endl;
                        seq += getSubsequence(reference2Genome, blocksForMsa->getAlignmentBlocks()[1].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[1].getStart() - 1);
                        seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[1].getStart()) {
                        seq += getReverseComplementary(blocksForMsa->getAncestralSequence());
                    } else {
                        std::cout << "line 487 unknown" << std::endl;
                    }
                } else if ((blocksForMsa->getWhichRef().compare("ref1") == 0 || blocksForMsa->getWhichRef().compare("NA") == 0) && blocksForMsa->getAlignmentBlocks()[1].getStrand() != blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) {
                    if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() < blocksForMsa->getAlignmentBlocks()[1].getStart() && POSITIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) {
                        if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {// tested
                            std::cout << "line 746\t" << blocksForMsa->getAlignmentBlocks()[0].getChr() << "\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getStart() << "\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd()
                                      << "\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 << "\t" << blocksForMsa->getAlignmentBlocks()[0].getStart() - 1 << std::endl;
                            //seq = getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getStart(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd());
                            std::cout << seq.size() << std::endl;
                            seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                        } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                        } else {
                            std::cout << "line 548 unknown" << std::endl;
                        }
                        seq += blocksForMsa->getAncestralSequence();
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() < blocksForMsa->getAlignmentBlocks()[1].getStart() && NEGATIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) {
                        if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {// tested
                            std::cout << "line 756" << std::endl;
                            seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                        } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                        } else {
                            std::cout << "line 557 unknown" << std::endl;
                        }
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() > blocksForMsa->getAlignmentBlocks()[1].getEnd() && NEGATIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) { //tested
                        std::cout << "line 816" << std::endl;
                        if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {
                            seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                        } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                        } else {
                            std::cout << "line 567 unknown" << std::endl;
                        }
                        seq += blocksForMsa->getAncestralSequence();
                    } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStart() > blocksForMsa->getAlignmentBlocks()[1].getEnd() && POSITIVE == blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()) { //tested
                        std::cout << "line 826" << std::endl;
                        if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 < blocksForMsa->getAlignmentBlocks()[0].getStart()) {
                            seq += getSubsequence(reference1Genome, blocksForMsa->getAlignmentBlocks()[0].getChr(), blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1, blocksForMsa->getAlignmentBlocks()[0].getStart() - 1);
                        } else if (blocksForMsa->getPrev()->getAlignmentBlocks()[0].getEnd() + 1 == blocksForMsa->getAlignmentBlocks()[0].getStart()) {

                        } else {
                            std::cout << "line 576 unknown" << std::endl;
                        }
                        seq += blocksForMsa->getAncestralSequence();
                        usingRef1Strand = true;
                    } else {
                        std::cout << "line 581 unknown\t" << blocksForMsa->getWhichRef() << std::endl;
                    }

                } else {
                    std::cout << "line 585 unknown\t" << blocksForMsa->getWhichRef() << "\t" << blocksForMsa->getPrev()->getAlignmentBlocks()[1].getEnd() << "\t" << blocksForMsa->getAlignmentBlocks()[1].getStart() << "\t"
                              << blocksForMsa->getPrev()->getAlignmentBlocks()[1].getStrand()
                              << "\t" << blocksForMsa->getAlignmentBlocks()[1].getStrand() << "\t" << usingRef1Strand << std::endl;
                }
            }
            oCfile << ">Chr" << index << std::endl;
            oCfile << seq << std::endl;
            ++index;
        }
    }
    oCfile.close();
}
