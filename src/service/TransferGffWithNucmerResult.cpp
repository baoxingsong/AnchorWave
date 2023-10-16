//
// Created by song on 8/4/18.
//

#include "TransferGffWithNucmerResult.h"


void readSam(std::vector<AlignmentMatch> &alignmentMatchsMapT, std::ifstream &infile,
             std::map<std::string, Transcript> &transcriptHashMap, int &expectCopy, const double &minimumSimilarity,
             double &secondarySimilarity, std::set<std::string> &blackGeneList, const std::string &anchorSequenceFile,
             int32_t &matchingScore, int32_t &mismatchingPenalty, int32_t &openGapPenalty1,
             int32_t &extendGapPenalty1, int &k, bool &H, int &w,
             std::map<std::string, std::tuple<std::string, long, long, int> > &queryGenome) {

    std::map<std::string, std::tuple<std::string, long, long, int> > anchorSequences2;
    readFastaFile(anchorSequenceFile, anchorSequences2);

    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;

    std::string databaseChr;
    int32_t databaseStart;
    int32_t databaseEnd;
    std::string queryChr;
    int32_t queryStart;
    int32_t queryEnd;

    std::map<std::string, std::string> lastChr;
    std::map<std::string, int32_t> lastPosition;

    std::map<std::string, std::map<std::string, std::vector<double>> > geneScores; //fist key is gene name, second key is chr value is a vector of similarity

    while (std::getline(infile, line)) { // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if (line.substr(0, 3) == "@PG") {
            std::vector<std::string> elements;
            std::vector<std::string> elements2;
            char seperator = ' ';
            char seperator2 = ',';
            split(line, seperator, elements);
            if (elements[0].find("ID:minimap2") != std::string::npos) {
                std::cout << "using parameters detected from the input SAM file for novel anchors identification" << std::endl;
                for (size_t i = 0; i < elements.size(); ++i) {
                    std::string element = elements[i];
                    if (element == "-k" and (i + 1) < elements.size()) {
                        k = std::stoi(elements[i + 1]);
                        w = 0.666 * k;
                    }
                    if (element == "-A" and (i + 1) < elements.size()) {
                        matchingScore = std::stoi(elements[i + 1]);
                    }
                    if (element == "-B" and (i + 1) < elements.size()) {
                        mismatchingPenalty = std::stoi(elements[i + 1]);
                    }
                    if (element == "-O" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        openGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-E" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        extendGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-H") {
                        H = true;
                    }
                }
            }
        }

        if (line.size()>1 && line[0] != '@') { //ignore the header
//            elems.clear();
            std::vector<std::string> elems;
            split(line, delim, elems);
            queryStart = stoi(elems[3]);
            queryChr = elems[2];

            if (queryChr.compare("*") != 0 && transcriptHashMap.find(elems[0]) != transcriptHashMap.end() && queryGenome.find(queryChr) != queryGenome.end() ) { // ignore those none mapping records
                databaseChr = transcriptHashMap[elems[0]].getChromeSomeName();
                databaseStart = transcriptHashMap[elems[0]].getPStart();
                databaseEnd = transcriptHashMap[elems[0]].getPEnd();
                queryEnd = queryStart; // this 1 based position

                int samFlag = stoi(elems[1]);

                double score = 0;
                int32_t currentCDSPosition = 0; // 0 based position
                int32_t currentqueryPosition = queryStart - 1; // 0 based position

                std::string cdsSequence = getSubsequence2(anchorSequences2, elems[0]);

                if (16 == samFlag % 32) {
                    cdsSequence = getReverseComplementary(cdsSequence);
                }

                std::vector<std::string> cigarElems;
                splitCIGAR(elems[5], cigarElems);
                int headClipping = 0;
                int tailClipping = 0;
                int numberofMatch = 0;
                std::string qseq;

                for (size_t i = 0; i < cigarElems.size(); ++i) {
                    std::string cVal = cigarElems[i];
                    char cLetter = cVal[cVal.length() - 1];
                    int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                    if (i == cigarElems.size() - 1 && (cLetter == 'H' || cLetter == 'S')) { // ignore the last soft/hard clipping
                        tailClipping += cLen;
                        continue;
                    }

                    switch (cLetter) {
                        case 'H':
                            headClipping += cLen;
                            currentCDSPosition += cLen;
                            break;
                        case 'S':
                            headClipping += cLen;
                            currentCDSPosition += cLen;
                            break;
                        case 'M':
                            queryEnd += cLen;
                            numberofMatch += cLen;
                            qseq = getSubsequence2(queryGenome, queryChr, currentqueryPosition, currentqueryPosition + cLen - 1);

                            for (int32_t p = 0; p < cLen; p++) {
                                if (cdsSequence[currentCDSPosition] == qseq[p] ) {
                                    score += matchingScore;
                                } else {
                                    score += mismatchingPenalty;
                                }
                                ++currentCDSPosition;
                                ++currentqueryPosition;
                            }
                            break;
                        case '=':
                            queryEnd += cLen;
                            numberofMatch += cLen;
                            score += cLen * matchingScore;
                            currentCDSPosition += cLen;
                            currentqueryPosition += cLen;
                            break;
                        case 'X':
                            numberofMatch += cLen;
                            queryEnd += cLen;
                            score += cLen * mismatchingPenalty;
                            currentCDSPosition += cLen;
                            currentqueryPosition += cLen;
                            break;
                        case 'I':
                            score += openGapPenalty1 + cLen * extendGapPenalty1;
                            currentCDSPosition += cLen;
                            break;
                        case 'D':
                            queryEnd += cLen;
                            score += openGapPenalty1 + cLen * extendGapPenalty1;
                            currentqueryPosition += cLen;
                            break;
                        case 'N': // do not use intron for score calculation
                            queryEnd += cLen;
                            currentqueryPosition += cLen;
                            break;
                        case 'P':
                            break;
                        default:
                            std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                            break;
                    }
                }

                --queryEnd;

                std::map<int32_t, int32_t> positionsMap;
                size_t cdsSequenceLength = 0;
                if (transcriptHashMap[elems[0]].getStrand() == POSITIVE) {
                    int32_t cdsPosition = 0;
                    int32_t chromosomePosition = transcriptHashMap[elems[0]].getPStart() - 1;
                    for (size_t cdsIndex = 0; cdsIndex < transcriptHashMap[elems[0]].getCdsVector().size(); ++cdsIndex) {
                        for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart(); i <= transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd(); ++i) {
                            cdsPosition++;
                            chromosomePosition++;
                            cdsSequenceLength++;
                            positionsMap[cdsPosition] = chromosomePosition;
                        }
                        if (cdsIndex < transcriptHashMap[elems[0]].getCdsVector().size() - 1) { // for intron
                            for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd() + 1; i < transcriptHashMap[elems[0]].getCdsVector()[cdsIndex + 1].getStart(); ++i) {
                                chromosomePosition++;
                            }
                        }
                    }
                    assert(chromosomePosition == transcriptHashMap[elems[0]].getPEnd());
                } else {
                    int32_t cdsPosition = 0;
                    int32_t chromosomePosition = transcriptHashMap[elems[0]].getPEnd() + 1;
                    for (int32_t cdsIndex = transcriptHashMap[elems[0]].getCdsVector().size() - 1; cdsIndex >= 0; --cdsIndex) {
                        for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd(); i >= transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart(); --i) {
                            cdsPosition++;
                            chromosomePosition--;
                            cdsSequenceLength++;
                            positionsMap[cdsPosition] = chromosomePosition;
                        }
                        if (cdsIndex > 0) { // for intron
                            for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart() - 1; i > transcriptHashMap[elems[0]].getCdsVector()[cdsIndex - 1].getEnd(); --i) {
                                chromosomePosition--;
                            }
                        }
                    }
                    assert(transcriptHashMap[elems[0]].getPStart() == chromosomePosition);
                }

                if (0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == POSITIVE) {
                    databaseStart = positionsMap[headClipping + 1];
                    databaseEnd = positionsMap[cdsSequenceLength - tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if (0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == NEGATIVE) {
                    databaseEnd = positionsMap[headClipping + 1];
                    databaseStart = positionsMap[cdsSequenceLength - tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == POSITIVE) {
                    databaseEnd = positionsMap[cdsSequenceLength - headClipping];
                    databaseStart = positionsMap[1 + tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == NEGATIVE) {
                    databaseStart = positionsMap[cdsSequenceLength - headClipping];
                    databaseEnd = positionsMap[1 + tailClipping];
                    assert(databaseStart < databaseEnd);
                }

                if (lastChr.find(elems[0]) != lastChr.end() && lastChr[elems[0]] == queryChr &&
                    std::min((std::abs(lastPosition[elems[0]] - queryEnd)), std::abs(lastPosition[elems[0]] - queryStart)) < std::abs(transcriptHashMap[elems[0]].getPStart() - transcriptHashMap[elems[0]].getPEnd())) {
                    blackGeneList.insert(elems[0]);
                } // remove those genes generated weired alignment

                lastChr[elems[0]] = queryChr;
                lastPosition[elems[0]] = queryEnd;

                //double thisScore = 1.0 - (tailClipping+headClipping)/(double)cdsSequenceLength;
                double thisScore = (double) numberofMatch / (double) cdsSequenceLength;
                //double thisScore = score;
                if (thisScore > minimumSimilarity) {
                    if ((0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == POSITIVE)
                        || (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == NEGATIVE)) {

                        AlignmentMatch orthologPair(databaseChr, queryChr, databaseStart, databaseEnd, queryStart, queryEnd, thisScore, POSITIVE, elems[0], elems[0]);
                        alignmentMatchsMapT.push_back(orthologPair);
                    } else {
                        AlignmentMatch orthologPair(databaseChr, queryChr, databaseStart, databaseEnd, queryStart, queryEnd, thisScore, NEGATIVE, elems[0], elems[0]);
                        alignmentMatchsMapT.push_back(orthologPair);
                    }
                    if (geneScores.find(elems[0]) == geneScores.end()) {
                        std::map<std::string, std::vector<double>> a;
                        geneScores[elems[0]] = a;
                    }
                    if (geneScores[elems[0]].find(queryChr) == geneScores[elems[0]].end()) {
                        std::vector<double> a;
                        geneScores[elems[0]][queryChr] = a;
                    }
                    geneScores[elems[0]][queryChr].push_back(thisScore);
                }
            }
        }
    }

    std::vector<int> teRemoveIndexes;
    if (expectCopy > 0) {
        for (std::map<std::string, std::map<std::string, std::vector<double>>>::iterator it0 = geneScores.begin(); it0 != geneScores.end(); ++it0) {
            std::string geneName = it0->first;
            for (std::map<std::string, std::vector<double>>::iterator it1 = geneScores[geneName].begin(); it1 != geneScores[geneName].end(); ++it1) {
                std::vector<double> scores = it1->second;
                std::sort(scores.begin(), scores.end());
                std::reverse(scores.begin(), scores.end());
                if (scores.size() > expectCopy && scores[expectCopy] / scores[0] > secondarySimilarity) {
                    blackGeneList.insert(geneName);
                }
            }
        }
    }

    for (size_t i = 0; i < alignmentMatchsMapT.size(); ++i) {
        std::string geneName = alignmentMatchsMapT[i].getReferenceGeneName();
        if (blackGeneList.find(geneName) != blackGeneList.end()) {
            teRemoveIndexes.push_back(i);
        }
    }

    for (int j = teRemoveIndexes.size() - 1; j >= 0; --j) {
        alignmentMatchsMapT.erase(alignmentMatchsMapT.begin() + teRemoveIndexes[j]);
    }

    if (alignmentMatchsMapT.size() == 0) {
        std::cout << "there is no match anchor found in the input sam file" << std::endl;
        std::exit(1);
    }
}

void readSam(std::vector<AlignmentMatch> &alignmentMatchsMapT, std::ifstream &infile, std::map<std::string, Transcript> &transcriptHashMap, int &expectCopy, const double &minimumSimilarity,
             double &secondarySimilarity, std::set<std::string> &blackGeneList,
             int32_t &matchingScore, int32_t &mismatchingPenalty, int32_t &openGapPenalty1,
             int32_t &extendGapPenalty1, int &k, bool &H, int &w) {

    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;

    std::string databaseChr;
    int32_t databaseStart;
    int32_t databaseEnd;
    std::string queryChr;
    int32_t queryStart;
    int32_t queryEnd;

    std::map<std::string, std::string> lastChr;
    std::map<std::string, int32_t> lastPosition;

    std::map<std::string, std::map<std::string, std::vector<double>>> geneScores; //first key is gene name, second key is chr value is a vector of similarity

    while (std::getline(infile, line)) { // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if (line.substr(0, 3) == "@PG") {
            std::vector<std::string> elements;
            std::vector<std::string> elements2;
            char seperator = ' ';
            char seperator2 = ',';
            split(line, seperator, elements);

            if (elements[0].find("ID:minimap2") != std::string::npos) {
//                std::cout << "using parameters detected from the input SAM file for novel anchors identification" << std::endl;
                for (size_t i = 0; i < elements.size(); ++i) {
                    std::string element = elements[i];
                    if (element == "-k" and (i + 1) < elements.size()) {
                        k = std::stoi(elements[i + 1]);
                        w = 0.666 * k;
                    }
                    if (element == "-A" and (i + 1) < elements.size()) {
                        matchingScore = std::stoi(elements[i + 1]);
                    }
                    if (element == "-B" and (i + 1) < elements.size()) {
                        mismatchingPenalty = std::stoi(elements[i + 1]);
                    }
                    if (element == "-O" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        openGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-E" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        extendGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-H") {
                        H = true;
                    }
                }
            }
        }

        if (line.size()>1 && line[0] != '@') { //ignore the header
            elems.clear();
            split(line, delim, elems);
            queryStart = stoi(elems[3]);
            queryChr = elems[2];
            if (queryChr.compare("*") != 0 && transcriptHashMap.find(elems[0]) != transcriptHashMap.end()) { // ignore those none mapping records
                databaseChr = transcriptHashMap[elems[0]].getChromeSomeName();

                databaseStart = transcriptHashMap[elems[0]].getPStart();
                databaseEnd = transcriptHashMap[elems[0]].getPEnd();
                queryEnd = queryStart;
                std::vector<std::string> cigarElems;
                splitCIGAR(elems[5], cigarElems);
                int headClipping = 0;
                int tailClipping = 0;
                int numberofMatch = 0;
                for (size_t i = 0; i < cigarElems.size(); ++i) {
                    std::string cVal = cigarElems[i];
                    char cLetter = cVal[cVal.length() - 1];
                    int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                    if (i == cigarElems.size() - 1 && (cLetter == 'H' || cLetter == 'S')) { // ignore the last soft/hard clipping
                        tailClipping += cLen;
                        continue;
                    }
                    switch (cLetter) {
                        case 'H':
                            headClipping += cLen;
                            break;
                        case 'S':
                            headClipping += cLen;
                            break;
                        case 'M':
                            queryEnd += cLen;
                            numberofMatch += cLen;
                            break;
                        case '=':
                            queryEnd += cLen;
                            numberofMatch += cLen;
                            break;
                        case 'X':
                            numberofMatch += cLen;
                            queryEnd += cLen;
                            break;
                        case 'I':
                            break;
                        case 'D':
                            queryEnd += cLen;
                            break;
                        case 'N':
                            queryEnd += cLen;
                            break;
                        case 'P':
                            break;
                        default:
                            std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                            break;
                    }
                }

                --queryEnd;
                int samFlag = stoi(elems[1]);

                std::map<int32_t, int32_t> positionsMap;
                size_t cdsSequenceLength = 0;
                if (transcriptHashMap[elems[0]].getStrand() == POSITIVE) {
                    int32_t cdsPosition = 0;
                    int32_t chromosomePosition = transcriptHashMap[elems[0]].getPStart() - 1;
                    for (size_t cdsIndex = 0; cdsIndex < transcriptHashMap[elems[0]].getCdsVector().size(); ++cdsIndex) {
                        for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart(); i <= transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd(); ++i) {
                            cdsPosition++;
                            chromosomePosition++;
                            cdsSequenceLength++;
                            positionsMap[cdsPosition] = chromosomePosition;
                        }
                        if (cdsIndex < transcriptHashMap[elems[0]].getCdsVector().size() - 1) { // for intron
                            for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd() + 1; i < transcriptHashMap[elems[0]].getCdsVector()[cdsIndex + 1].getStart(); ++i) {
                                chromosomePosition++;
                            }
                        }
                    }

                    assert(chromosomePosition == transcriptHashMap[elems[0]].getPEnd());
                } else {
                    int32_t cdsPosition = 0;
                    int32_t chromosomePosition = transcriptHashMap[elems[0]].getPEnd() + 1;
                    for (int32_t cdsIndex = transcriptHashMap[elems[0]].getCdsVector().size() - 1; cdsIndex >= 0; --cdsIndex) {
                        for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd(); i >= transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart(); --i) {
                            cdsPosition++;
                            chromosomePosition--;
                            cdsSequenceLength++;
                            positionsMap[cdsPosition] = chromosomePosition;
                        }
                        if (cdsIndex > 0) { // for intron
                            for (int32_t i = transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart() - 1; i > transcriptHashMap[elems[0]].getCdsVector()[cdsIndex - 1].getEnd(); --i) {
                                chromosomePosition--;
                            }
                        }
                    }

                    assert(transcriptHashMap[elems[0]].getPStart() == chromosomePosition);
                }

                if (0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == POSITIVE) {
                    databaseStart = positionsMap[headClipping + 1];
                    databaseEnd = positionsMap[cdsSequenceLength - tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if (0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == NEGATIVE) {
                    databaseEnd = positionsMap[headClipping + 1];
                    databaseStart = positionsMap[cdsSequenceLength - tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == POSITIVE) {
                    databaseEnd = positionsMap[cdsSequenceLength - headClipping];
                    databaseStart = positionsMap[1 + tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == NEGATIVE) {
                    databaseStart = positionsMap[cdsSequenceLength - headClipping];
                    databaseEnd = positionsMap[1 + tailClipping];
                    assert(databaseStart < databaseEnd);
                }

                if (lastChr.find(elems[0]) != lastChr.end() && lastChr[elems[0]] == queryChr &&
                    std::min((std::abs(lastPosition[elems[0]] - queryEnd)), std::abs(lastPosition[elems[0]] - queryStart)) < std::abs(transcriptHashMap[elems[0]].getPStart() - transcriptHashMap[elems[0]].getPEnd())) {
                    blackGeneList.insert(elems[0]);
                } // remove those genes generated weired alignment

                lastChr[elems[0]] = queryChr;
                lastPosition[elems[0]] = queryEnd;

                //double thisScore = 1.0 - (tailClipping+headClipping)/(double)cdsSequenceLength;
                double thisScore = (double) numberofMatch / (double) cdsSequenceLength;
                if (thisScore > minimumSimilarity) {
                    if ((0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == POSITIVE) || (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand() == NEGATIVE)) {
                        AlignmentMatch orthologPair(databaseChr, queryChr, databaseStart, databaseEnd, queryStart, queryEnd, thisScore, POSITIVE, elems[0], elems[0]);
                        alignmentMatchsMapT.push_back(orthologPair);
                    } else {
                        AlignmentMatch orthologPair(databaseChr, queryChr, databaseStart, databaseEnd, queryStart, queryEnd, thisScore, NEGATIVE, elems[0], elems[0]);
                        alignmentMatchsMapT.push_back(orthologPair);
                    }

                    if (geneScores.find(elems[0]) == geneScores.end()) {
                        std::map<std::string, std::vector<double>> a;
                        geneScores[elems[0]] = a;
                    }

                    if (geneScores[elems[0]].find(queryChr) == geneScores[elems[0]].end()) {
                        std::vector<double> a;
                        geneScores[elems[0]][queryChr] = a;
                    }
                    geneScores[elems[0]][queryChr].push_back(thisScore);
                }
            }
        }
    }

    std::vector<int> teRemoveIndexes;
    if (expectCopy > 0) {
        for (std::map<std::string, std::map<std::string, std::vector<double>>>::iterator it0 = geneScores.begin(); it0 != geneScores.end(); ++it0) {  // remove those genes with too much copies
            std::string geneName = it0->first;
            for (std::map<std::string, std::vector<double>>::iterator it1 = geneScores[geneName].begin(); it1 != geneScores[geneName].end(); ++it1) {
                std::vector<double> scores = it1->second;
                std::sort(scores.begin(), scores.end());
                std::reverse(scores.begin(), scores.end());
                if (scores.size() > expectCopy && scores[expectCopy] / scores[0] > secondarySimilarity) {
                    blackGeneList.insert(geneName);
//                    std::cout << "removing " << geneName << " due to too much copies. " << it1->first << "\t" << scores.size() << "\t" << scores[0] << "\t" << scores[expectCopy] << std::endl;
                }
            }
        }
    }

    for (size_t i = 0; i < alignmentMatchsMapT.size(); ++i) {
        std::string geneName = alignmentMatchsMapT[i].getReferenceGeneName();
        if (blackGeneList.find(geneName) != blackGeneList.end()) {
            teRemoveIndexes.push_back(i);
        }
    }

    for (int j = teRemoveIndexes.size() - 1; j >= 0; --j) {
        alignmentMatchsMapT.erase(alignmentMatchsMapT.begin() + teRemoveIndexes[j]);
    }

    if (alignmentMatchsMapT.size() == 0) {
        std::cout << "there is no match anchor found in the input sam file" << std::endl;
        std::exit(1);
    }
}

void setupAnchorsWithSpliceAlignmentResult(const std::string &gffFilePath, const std::string &cdsSequenceFile, const std::string &samFile, std::map<std::string, std::vector<AlignmentMatch>> &map_v_am,
                                           double &inversion_PENALTY, double &MIN_ALIGNMENT_SCORE, bool &considerInversion, const int &minExon, const int64_t &windowWidth, const double &minimumSimilarity, const double &minimumSimilarity2,
                                           std::map<std::string, std::tuple<std::string, long, long, int> > &map_ref,
                                           std::map<std::string, std::tuple<std::string, long, long, int> > &map_qry,
                                           int &expectedCopies, double &maximumSimilarity,
                                           const std::string &referenceSamFilePath, const int32_t &wfaSize3, const bool &searchForNewAnchors, const bool &exonModel) {
    std::ifstream infile(samFile);
    if (!infile.good()) {
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit(1);
    }

    // those are default parameter from minimap2
    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = 4;
    int32_t openGapPenalty1 = 4;
    int32_t extendGapPenalty1 = 2;
    int k = 15;
    int w = 0.666 * k;
    bool H = false;

    //read genome and gff file begin
    std::map<std::string, std::vector<Transcript> > map_v_ts;
//    if (exonModel) {
//        readGffFile(gffFilePath, map_v_ts, "exon", minExon);
//    }
//    else {
//        readGffFile(gffFilePath, map_v_ts, "CDS", minExon);
//    }

    std::string regex = "([\\s\\S]*)Parent=([a-zA-Z0-9.:_-]+)";
    if (exonModel) {
        readGffFile_exon(gffFilePath, map_v_ts, regex, minExon);
    }
    else {
        readGffFile1(gffFilePath, map_v_ts, regex, minExon);
    }


    std::set<std::string> set_rm_chr;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = map_v_ts.begin(); it != map_v_ts.end(); ++it) {
        if (map_ref.find(it->first) == map_ref.end()) {
            set_rm_chr.insert(it->first);
        }
        if (map_qry.find(it->first) == map_qry.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_ref.begin(); it != map_ref.end(); ++it) {
        if (map_v_ts.find(it->first) == map_v_ts.end()) {
            set_rm_chr.insert(it->first);
        }
        if (map_qry.find(it->first) == map_qry.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_qry.begin(); it != map_qry.end(); ++it) {
        if (map_v_ts.find(it->first) == map_v_ts.end()) {
            set_rm_chr.insert(it->first);
        }
        if (map_ref.find(it->first) == map_ref.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::string chr: set_rm_chr) {
        std::cerr << "There is not enough anchors found on " << chr << std::endl;
        if (map_v_ts.find(chr) != map_v_ts.end()) {
            map_v_ts.erase(chr);
        }
        if (map_ref.find(chr) != map_ref.end()) {
            map_ref.erase(chr);
        }
        if (map_qry.find(chr) != map_qry.end()) {
            map_qry.erase(chr);
        }
    }

    std::map<std::string, Transcript> transcriptHashMap; // key is transcript name, value is a transcript structure
    for (std::map<std::string, std::vector<Transcript> >::iterator it = map_v_ts.begin(); it != map_v_ts.end(); ++it) {
        for (Transcript transcript: it->second) {
            transcriptHashMap[transcript.getName()] = transcript;
        }
    }

    //read genome and gff file end

    // set gene black list by reading reference gff file begin
    std::set<std::string> blackGeneList;
    if (referenceSamFilePath.size() > 0) {
        std::ifstream infileReferencSam(referenceSamFilePath);
        if (!infileReferencSam.good()) {
            std::cerr << "error in opening sam file " << referenceSamFilePath << std::endl;
            exit(1);
        }

        std::vector<AlignmentMatch> alignmentMatchsMapT0;
        std::cout << "reading reference sam begin" << std::endl;
        readSam(alignmentMatchsMapT0, infileReferencSam, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w);
        std::cout << "reading reference sam done" << std::endl;
        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;

        for (AlignmentMatch orthologPair2: alignmentMatchsMapT0) {
            if (alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end()) {
                alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
            }

            if (orthologPair2.getRefChr() == orthologPair2.getQueryChr() && orthologPair2.getStrand() == POSITIVE) {
                alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(orthologPair2);
            }
        }

        bool keepTandemDuplication = false;
        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMapT.begin(); it != alignmentMatchsMapT.end(); ++it) {
            myAlignmentMatchSort(it->second, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, false);
            std::vector<AlignmentMatch> sortedAlignmentMatchs;

            if (it->second.size() > 1) {
                longestPath(it->second, sortedAlignmentMatchs, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
            }
            else {
                sortedAlignmentMatchs = it->second;
            }

            map_v_am[it->first] = std::vector<AlignmentMatch>();

            for (unsigned long i = 0; i < sortedAlignmentMatchs.size(); ++i) {
                map_v_am[it->first].push_back(sortedAlignmentMatchs[i]);
            }
        }

        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
            for (size_t rangeIndex = 0; rangeIndex < it->second.size(); ++rangeIndex) {
                if (it->second[rangeIndex].getRefStartPos() != it->second[rangeIndex].getQueryStartPos() || it->second[rangeIndex].getRefEndPos() != it->second[rangeIndex].getQueryEndPos()) {
                    blackGeneList.insert(it->second[rangeIndex].getReferenceGeneName());
                }
            }
        }

        infileReferencSam.close();
        map_v_am.clear();
    }

    // set gene black list by reading reference gff file end

    std::vector<AlignmentMatch> alignmentMatchsMapT0;
    std::cout << "reading qry sam begin" << std::endl;
    readSam(alignmentMatchsMapT0, infile, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, cdsSequenceFile, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w, map_qry);
    std::cout << "reading qry sam end" << std::endl;

    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;
    for (AlignmentMatch orthologPair2: alignmentMatchsMapT0) {
        if (alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end()) {
            alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
        }

        if (orthologPair2.getRefChr() == orthologPair2.getQueryChr() && map_v_ts.find(orthologPair2.getRefChr()) != map_v_ts.end() &&
                map_ref.find(orthologPair2.getRefChr()) != map_ref.end() && map_qry.find(orthologPair2.getRefChr()) != map_qry.end()) {
            if (!considerInversion && orthologPair2.getStrand() == NEGATIVE) {

            } else {
                alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(orthologPair2);
            }
        }
    }

    bool keepTandemDuplication = false;
    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_ref.begin(); it != map_ref.end(); ++it) {
        if ( alignmentMatchsMapT.find( it->first ) == alignmentMatchsMapT.end() ){
            std::cerr << "There is not enough anchors found on " << it->first << std::endl;
        }
    }
    for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMapT.begin(); it != alignmentMatchsMapT.end(); ++it) {
        if (it->second.size() < 3) {
            std::cerr << "There is not enough anchors found on " << it->first << std::endl;
            continue;
        }

        std::vector<AlignmentMatch> temp;

        myAlignmentMatchSort(it->second, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);
        if (it->second.size() > 1) {
            longestPath(it->second, temp, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
        } else {
            temp = it->second;
        }

        int is_hpc = 0; // no, do not use  homopolymer-compressed (HPC) minimizers.
        if (H) {
            is_hpc = 1;
        }

        int bucket_bits = 2;
        int n = 1;
        bool changed = false;
        if (searchForNewAnchors) {
            changed = true;
        }

        std::set<int32_t> blackList;
        while (changed) {
            std::vector<AlignmentMatch> temp2;
            changed = false;

            size_t startRef = 1;
            size_t startQuery = 1;
            size_t endRef;
            size_t endQuery;

            std::string refChr = it->second[0].getRefChr();
            std::string queryChr = refChr;

            STRAND lastStrand = POSITIVE;

            bool hasInversion = false;
            int32_t temp_size = temp.size();
            myAlignmentMatchSort(temp, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, false);
            for (int32_t m = 0; m < temp_size; ++m) {
                AlignmentMatch alignmentMatch = temp[m];
                if (alignmentMatch.getStrand() == NEGATIVE) {
                    hasInversion = true;
                }

                if (lastStrand == POSITIVE && alignmentMatch.getStrand() == POSITIVE) {
                    if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() != startQuery) {
                        endQuery = alignmentMatch.getQueryStartPos() - 1;
                    } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryStartPos() == startQuery) {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                    } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() == startQuery) {

                    } else {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        endQuery = alignmentMatch.getQueryStartPos() - 1;

                        std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                        std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                        if ((refSeq.size() * querySeq.size() > wfaSize3 * wfaSize3) && (refSeq.size() > k && querySeq.size() > k) && blackList.find(startRef) == blackList.end()) {

                            mm_idxopt_t iopt;
                            mm_mapopt_t mopt;
                            mm_verbose = 2; // disable message output to stderr
                            mm_set_opt(0, &iopt, &mopt);
                            mopt.flag |= MM_F_CIGAR; // DO NOT perform alignment
                            mopt.flag |= MM_F_NO_PRINT_2ND;
                            mopt.bw = windowWidth / 5;
                            mopt.flag |= MM_F_NO_LJOIN; // together the last one, control the maximum gap length on the local alignment region (novel seed)
                            mopt.a = matchingScore;
                            mopt.b = mismatchingPenalty;
                            mopt.q = openGapPenalty1;
                            mopt.e = extendGapPenalty1;
                            mopt.q2 = openGapPenalty1;
                            mopt.e2 = extendGapPenalty1;
                            mopt.bw = k;
                            mopt.mid_occ = 2; // ignore seeds with occurrences above this threshold
                            mopt.min_cnt = 2;// min number of minimizers on each chain
                            int32_t referenceSeqLength = refSeq.length();
                            int32_t querySeqLength = querySeq.length();
                            char *reference_seq_array = new char[referenceSeqLength + 1];
                            char *query_seq_array = new char[querySeqLength + 1];
                            strcpy(reference_seq_array, refSeq.c_str());

                            const char *refseq[1] = {reference_seq_array};
                            strcpy(query_seq_array, querySeq.c_str());

                            std::string queryName = "temp";
                            int nameLength = queryName.length();
                            // declaring character array
                            char name_array[nameLength + 1];
                            // copying the contents of the string to char array
                            strcpy(name_array, queryName.c_str());
                            const char *name[1] = {"temp"};

                            mm_idx_t *mi = mm_idx_str(w, k, is_hpc, bucket_bits, n, refseq, name);
                            mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
                            mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
                            int j, n_reg;

                            mm_reg1_t *reg = mm_map(mi, querySeqLength, query_seq_array, &n_reg, tbuf, &mopt, name_array); // get all hits for the query

                            if (n_reg > 0 && (&reg[0])->rev == 0) {
                                mm_reg1_t *r = &reg[0];
//                                assert(r->p); // with MM_F_CIGAR, this should not be NULL
                                int32_t newAnchorRefEnd = r->re - 1;
//                                int32_t newAnchorRefStart = startRef + r->rs;
                                int32_t newAnchorQueryEnd = r->qe - 1;
                                std::string alignmentName = "localAlignment_" + refChr + "_" + std::to_string(startRef + r->rs) + "_" + std::to_string(startRef + newAnchorRefEnd);
                                double length = r->re - r->rs + 1.0;
                                double numberofMs = 0.0;
                                for (uint32_t i = 0; i < r->p->n_cigar; ++i) {
                                    if ("MIDNSH"[r->p->cigar[i] & 0xf] == 'M') {
                                        int32_t thisLength = r->p->cigar[i] >> 4;
                                        numberofMs = numberofMs + thisLength;
                                    }
                                }

                                double similarity = numberofMs / length;
                                if (similarity > minimumSimilarity2) {
                                    changed = true;
                                    double thisScore = r->score / 2;
                                    thisScore = thisScore / ((r->re - r->rs) * matchingScore);
                                    AlignmentMatch orthologPair(refChr, queryChr,
                                                                startRef + r->rs, startRef + newAnchorRefEnd, startQuery + r->qs,
                                                                startQuery + newAnchorQueryEnd, thisScore, POSITIVE, alignmentName,
                                                                alignmentName);
                                    temp2.push_back(orthologPair);
                                }
                            } else {
                                blackList.insert(startRef);
                            }

                            for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                                mm_reg1_t *r = &reg[j];
                                free(r->p);
                            }

                            free(reg);
                            mm_tbuf_destroy(tbuf);
                            mm_idx_destroy(mi);
                            delete [] reference_seq_array;
                            delete [] query_seq_array;
                        }
                    }
                } else if (lastStrand == NEGATIVE && alignmentMatch.getStrand() == NEGATIVE) {
                    if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() != startQuery) {

                    } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryEndPos() == startQuery) {

                    } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() == startQuery) {

                    } else {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        endQuery = alignmentMatch.getQueryEndPos() + 1;

                        std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                        std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                        if ((refSeq.size() * querySeq.size() > wfaSize3 * wfaSize3) && (refSeq.size() > k && querySeq.size() > k) && blackList.find(startRef) == blackList.end()) {
                            mm_idxopt_t iopt;
                            mm_mapopt_t mopt;

                            mm_verbose = 2; // disable message output to stderr
                            mm_set_opt(0, &iopt, &mopt);
                            //mopt.flag &= ~ MM_F_CIGAR; // DO NOT perform alignment
                            mopt.flag |= MM_F_CIGAR; // DO NOT perform alignment
                            mopt.flag |= MM_F_NO_PRINT_2ND;
                            mopt.bw = windowWidth / 5;
                            mopt.flag |= MM_F_NO_LJOIN; // together the last one, control the maximum gap length on the local alignment region (novel seed)
                            mopt.a = matchingScore;
                            mopt.b = mismatchingPenalty;
                            mopt.q = openGapPenalty1;
                            mopt.e = extendGapPenalty1;
                            mopt.q2 = openGapPenalty1;
                            mopt.e2 = extendGapPenalty1;
                            mopt.bw = k;
                            mopt.mid_occ = 2; // ignore seeds with occurrences above this threshold
                            mopt.min_cnt = 2;// min number of minimizers on each chain
                            int32_t referenceSeqLength = refSeq.length();
                            int32_t querySeqLength = querySeq.length();
                            char *reference_seq_array = new char[referenceSeqLength + 1];
                            char *query_seq_array = new char[querySeqLength + 1];

                            strcpy(reference_seq_array, refSeq.c_str());
                            const char *refseq[1] = {reference_seq_array};
                            strcpy(query_seq_array, querySeq.c_str());

                            std::string queryName = "temp";
                            int nameLength = queryName.length();
                            // declaring character array
                            char name_array[nameLength + 1];
                            // copying the contents of the string to char array
                            strcpy(name_array, queryName.c_str());
                            const char *name[1] = {"temp"};

                            mm_idx_t *mi = mm_idx_str(w, k, is_hpc, bucket_bits, n, refseq, name);

                            mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
                            mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
                            int j, n_reg;

                            mm_reg1_t *reg = mm_map(mi, querySeqLength, query_seq_array, &n_reg, tbuf, &mopt, name_array); // get all hits for the query
                            if (n_reg > 0 && (&reg[0])->rev == 0) {
                                mm_reg1_t *r = &reg[0];
//                                assert(r->p); // with MM_F_CIGAR, this should not be NULL
                                int32_t newAnchorRefEnd = r->re - 1;
//                                int32_t newAnchorRefStart = startRef + r->rs;
                                double length = r->re - r->rs + 1;
                                double numberofMs = 0;
                                for (uint32_t i = 0; i < r->p->n_cigar; ++i) {
                                    if ("MIDNSH"[r->p->cigar[i] & 0xf] == 'M') {
                                        int32_t thisLength = r->p->cigar[i] >> 4;
                                        numberofMs = numberofMs + thisLength;
                                    }
                                }

                                double similarity = numberofMs / length;
                                if (similarity > minimumSimilarity2) {
                                    changed = true;
                                    int32_t newAnchorQueryEnd = r->qe - 1;
                                    std::string alignmentName = "localAlignment_" + refChr + "_" + std::to_string(startRef + r->rs) + "_" + std::to_string(startRef + newAnchorRefEnd);
                                    double thisScore = r->score / 2;
                                    thisScore = thisScore / ((r->re - r->rs) * matchingScore);
                                    AlignmentMatch orthologPair(refChr, queryChr,
                                                                startRef + r->rs, startRef + newAnchorRefEnd,
                                                                startQuery - newAnchorQueryEnd, startQuery - r->qs,
                                                                thisScore, NEGATIVE,
                                                                alignmentName, alignmentName);
                                    temp2.push_back(orthologPair);
                                }
                            } else {
                                blackList.insert(startRef);
                            }

                            for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                                mm_reg1_t *r = &reg[j];
                                free(r->p);
                            }

                            free(reg);
                            mm_tbuf_destroy(tbuf);
                            mm_idx_destroy(mi);
                            delete [] reference_seq_array;
                            delete [] query_seq_array;
                        }
                    }
                }

                startRef = alignmentMatch.getRefEndPos() + 1;
                startQuery = alignmentMatch.getQueryEndPos() + 1;
                if (alignmentMatch.getStrand() == NEGATIVE) {
                    startQuery = alignmentMatch.getQueryStartPos() - 1;
                }
                lastStrand = alignmentMatch.getStrand();
            }

            if (!hasInversion) {
                endRef = getSequenceSizeFromPath2(map_ref[refChr]);
                endQuery = getSequenceSizeFromPath2(map_qry[queryChr]);

                if (startRef > endRef && startQuery <= endQuery) {

                } else if (startRef <= endRef && startQuery > endQuery) {

                } else if (startRef > endRef && startQuery > endQuery) {

                } else { // have bug by aligning b73 with mo17
                    std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                    if ((refSeq.size() * querySeq.size() > wfaSize3 * wfaSize3) &&
                        (refSeq.size() > k && querySeq.size() > k) && blackList.find(startRef) == blackList.end()) {
                        mm_idxopt_t iopt;
                        mm_mapopt_t mopt;
                        mm_verbose = 2; // disable message output to stderr
                        mm_set_opt(0, &iopt, &mopt);
                        //mopt.flag &= ~ MM_F_CIGAR; // DO NOT perform alignment
                        mopt.flag |= MM_F_CIGAR; // DO NOT perform alignment
                        mopt.flag |= MM_F_NO_PRINT_2ND;
                        mopt.bw = windowWidth / 5;
                        mopt.flag |= MM_F_NO_LJOIN; // together the last one, control the maximum gap length on the local alignment region (novel seed)
                        mopt.a = matchingScore;
                        mopt.b = mismatchingPenalty;
                        mopt.q = openGapPenalty1;
                        mopt.e = extendGapPenalty1;
                        mopt.q2 = openGapPenalty1;
                        mopt.e2 = extendGapPenalty1;
                        mopt.bw = k;
                        mopt.mid_occ = 2; // ignore seeds with occurrences above this threshold
                        mopt.min_cnt = 2;// min number of minimizers on each chain
                        int referenceSeqLength = refSeq.length();
                        int querySeqLength = querySeq.length();
                        char *reference_seq_array = new char[referenceSeqLength + 1];
                        char *query_seq_array = new char[querySeqLength + 1];

                        strcpy(reference_seq_array, refSeq.c_str());
                        const char *refseq[1] = {reference_seq_array};
                        strcpy(query_seq_array, querySeq.c_str());

                        std::string queryName = "temp";
                        int nameLength = queryName.length();
                        // declaring character array
                        char name_array[nameLength + 1];
                        // copying the contents of the string to char array
                        strcpy(name_array, queryName.c_str());
                        const char *name[1] = {"temp"};

                        mm_idx_t *mi = mm_idx_str(w, k, is_hpc, bucket_bits, n, refseq, name);

                        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
                        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
                        int j, n_reg;
                        mm_reg1_t *reg = mm_map(mi, querySeqLength, query_seq_array, &n_reg, tbuf, &mopt, name_array); // get all hits for the query

                        if (n_reg > 0 && (&reg[0])->rev == 0) {
                            mm_reg1_t *r = &reg[0];
//                            assert(r->p); // with MM_F_CIGAR, this should not be NULL
                            int32_t newAnchorRefEnd = r->re - 1;
//                            int32_t newAnchorRefStart = startRef + r->rs;
                            double length = r->re - r->rs + 1;
                            double numberofMs = 0;
                            for (uint32_t i = 0; i < r->p->n_cigar; ++i) {
                                if ("MIDNSH"[r->p->cigar[i] & 0xf] == 'M') {
                                    int32_t thisLength = r->p->cigar[i] >> 4;
                                    numberofMs = numberofMs + thisLength;
                                }
                            }

                            double similarity = numberofMs / length;
                            if (similarity > minimumSimilarity2) {
                                changed = true;
                                int32_t newAnchorQueryEnd = r->qe - 1;
                                std::string alignmentName = "localAlignment_" + refChr + "_" + std::to_string(startRef + r->rs) + "_" + std::to_string(startRef + newAnchorRefEnd);
                                double thisScore = r->score / 2;
                                thisScore = thisScore / ((r->re - r->rs) * matchingScore);
                                AlignmentMatch orthologPair(refChr, queryChr, startRef + r->rs,
                                                            startRef + newAnchorRefEnd, startQuery + r->qs,
                                                            startQuery + newAnchorQueryEnd, thisScore, POSITIVE,
                                                            alignmentName, alignmentName);
                                temp2.push_back(orthologPair);
                            }
                        } else {
                            blackList.insert(startRef);
                        }

                        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                            mm_reg1_t *r = &reg[j];
                            free(r->p);
                        }

                        free(reg);
                        mm_tbuf_destroy(tbuf);
                        mm_idx_destroy(mi);
                        delete [] reference_seq_array;
                        delete [] query_seq_array;
                    }
                }
            }

            for (AlignmentMatch alignmentMatch: temp2) {
                temp.push_back(alignmentMatch);
            }
        }

        if (temp.size() > 0) {
            myAlignmentMatchSort(temp, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);

            std::vector<AlignmentMatch> sortedAlignmentMatchs;

            if (temp.size() > 1) {
                longestPath(temp, sortedAlignmentMatchs, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
            } else {
                sortedAlignmentMatchs = temp;
            }

            for (unsigned long i = 0; i < sortedAlignmentMatchs.size(); ++i) {
                map_v_am[it->first].push_back(sortedAlignmentMatchs[i]);
            }
        }
    }

    infile.close();
}

void setupAnchorsWithSpliceAlignmentResultQuota(const std::string &gffFilePath, const std::string &samFile, const std::string &cdsSequenceFile, std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                                double &INDEL_SCORE, double &GAP_OPEN_PENALTY, double &MIN_ALIGNMENT_SCORE, int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                                                double &calculateIndelDistance, const int &minExon, const int64_t &windowWidth, const double &minimumSimilarity,
                                                const double &minimumSimilarity2,
                                                std::map<std::string, std::tuple<std::string, long, long, int> > &map_ref,
                                                std::map<std::string, std::tuple<std::string, long, long, int> > &map_qry,
                                                int &expectedCopies, const int32_t &wfaSize3,
                                                double &maximumSimilarity, const std::string &referenceSamFilePath, bool &searchForNewAnchors, const bool &exonModel) {

    // they are default parameter from minimap2
    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = 4;
    int32_t openGapPenalty1 = 4;
    int32_t extendGapPenalty1 = 2;
    int k = 15;
    int w = 0.666 * k;
    bool H = false;

    std::ifstream infile(samFile);
    if (!infile.good()) {
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit(1);
    }

    // read reference genome and gff begin
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
//    if (exonModel) {
//        readGffFile(gffFilePath, transcriptHashSet, "exon", minExon);
//    } else {
//        readGffFile(gffFilePath, transcriptHashSet, "CDS", minExon);
//    }

    std::string regex = "([\\s\\S]*)Parent=([a-zA-Z0-9.:_-]+)";
    if (exonModel) {
        readGffFile_exon(gffFilePath, transcriptHashSet, regex, minExon);
    }
    else {
        readGffFile1(gffFilePath, transcriptHashSet, regex, minExon);
    }

    // clean reference genome and annotation chromosomes
    std::set<std::string> toRemoveChrs;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = transcriptHashSet.begin();
         it != transcriptHashSet.end(); ++it) {
        if (map_ref.find(it->first) == map_ref.end()) {
            toRemoveChrs.insert(it->first);
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_ref.begin(); it != map_ref.end(); ++it) {
        if (transcriptHashSet.find(it->first) == transcriptHashSet.end()) {
            toRemoveChrs.insert(it->first);
        }
    }

    for (std::string chr: toRemoveChrs) {
        if (transcriptHashSet.find(chr) != transcriptHashSet.end()) {
            transcriptHashSet.erase(chr);
        }

        if (map_ref.find(chr) != map_ref.end()) {
            map_ref.erase(chr);
        }
    }

    std::map<std::string, Transcript> transcriptHashMap; // key is transcript name, value is a transcript structure
    for (std::map<std::string, std::vector<Transcript> >::iterator it = transcriptHashSet.begin();
         it != transcriptHashSet.end(); ++it) {
        for (Transcript transcript: it->second) {
            transcriptHashMap[transcript.getName()] = transcript;
        }
    }

    std::set<std::string> blackGeneList;
    if (referenceSamFilePath.size() > 0) {
        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap000;
        std::ifstream infileReferencSam(referenceSamFilePath);
        if (!infileReferencSam.good()) {
            std::cerr << "error in opening sam file " << referenceSamFilePath << std::endl;
            exit(1);
        }

        std::vector<AlignmentMatch> alignmentMatchsMapT0;
        readSam(alignmentMatchsMapT0, infileReferencSam, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w);
        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;

        for (AlignmentMatch orthologPair2: alignmentMatchsMapT0) {
            if (alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end()) {
                alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
            }
            if (orthologPair2.getRefChr() == orthologPair2.getQueryChr() && orthologPair2.getStrand() == POSITIVE) {
                alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(orthologPair2);
            }
        }

        bool keepTandemDuplication = false;
        double inversion_PENALTY = -1;
        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMapT.begin(); it != alignmentMatchsMapT.end(); ++it) {
            myAlignmentMatchSort(it->second, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, false); // since here the considerInversion is false, so the parameters of inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, would not be used
            std::vector<AlignmentMatch> sortedAlignmentMatchs;
            if (it->second.size() > 1) {
                longestPath(it->second, sortedAlignmentMatchs, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
            } else {
                sortedAlignmentMatchs = it->second;
            }

            alignmentMatchsMap000[it->first] = std::vector<AlignmentMatch>();
            for (unsigned long i = 0; i < sortedAlignmentMatchs.size(); ++i) {
                alignmentMatchsMap000[it->first].push_back(sortedAlignmentMatchs[i]);
            }
        }

        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMap000.begin(); it != alignmentMatchsMap000.end(); ++it) {
            for (size_t rangeIndex = 0; rangeIndex < it->second.size(); ++rangeIndex) {
                if (it->second[rangeIndex].getRefStartPos() != it->second[rangeIndex].getQueryStartPos() ||
                    it->second[rangeIndex].getRefEndPos() != it->second[rangeIndex].getQueryEndPos()) {
                    blackGeneList.insert(it->second[rangeIndex].getReferenceGeneName());
                }
            }
        }

        infileReferencSam.close();
        alignmentMatchsMap.clear();
        alignmentMatchsMap000.clear();
    }

    {
        std::vector<AlignmentMatch> alignmentMatchsMapT;
        readSam(alignmentMatchsMapT, infile, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, cdsSequenceFile, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w, map_qry);

        // begin setting index, they are necessary in the longest path approach
        std::map<std::string, int64_t> queryIndex;
        std::map<std::string, std::map<int64_t, AlignmentMatch>> queryIndexMap; // chr, index, AlignmentMatch
        myOrthologPairsSortQueryQuota(alignmentMatchsMapT);
        for (size_t ii = 0; ii < alignmentMatchsMapT.size(); ++ii) {
            if (queryIndex.find(alignmentMatchsMapT[ii].getQueryChr()) == queryIndex.end()) {
                queryIndex[alignmentMatchsMapT[ii].getQueryChr()] = 0;
                queryIndexMap[alignmentMatchsMapT[ii].getQueryChr()] = std::map<int64_t, AlignmentMatch>();
            } else {
                queryIndex[alignmentMatchsMapT[ii].getQueryChr()] = queryIndex[alignmentMatchsMapT[ii].getQueryChr()] + 1;
            }

            alignmentMatchsMapT[ii].setQueryId(queryIndex[alignmentMatchsMapT[ii].getQueryChr()]);
            queryIndexMap[alignmentMatchsMapT[ii].getQueryChr()][queryIndex[alignmentMatchsMapT[ii].getQueryChr()]] = alignmentMatchsMapT[ii];
        }

//        std::cout << "query id setting done" << std::endl;

        myOrthologPairsSortQuota(alignmentMatchsMapT);
        std::map<std::string, int64_t> refIndex; // key is chr
        std::map<std::string, std::map<int64_t, AlignmentMatch>> refIndexMap; // chr, index, AlignmentMatch
        for (size_t ii = 0; ii < alignmentMatchsMapT.size(); ++ii) {
            if (refIndex.find(alignmentMatchsMapT[ii].getRefChr()) == refIndex.end()) {
                refIndex[alignmentMatchsMapT[ii].getRefChr()] = 0;
                refIndexMap[alignmentMatchsMapT[ii].getRefChr()] = std::map<int64_t, AlignmentMatch>();
            }
            else {
                refIndex[alignmentMatchsMapT[ii].getRefChr()] = refIndex[alignmentMatchsMapT[ii].getRefChr()] + 1;
            }
            alignmentMatchsMapT[ii].setRefId(refIndex[alignmentMatchsMapT[ii].getRefChr()]);
            refIndexMap[alignmentMatchsMapT[ii].getRefChr()][refIndex[alignmentMatchsMapT[ii].getRefChr()]] = alignmentMatchsMapT[ii];
        }

        // index setting end
//        std::cout << "reference id setting done" << std::endl;

        longestPathQuotav2(alignmentMatchsMapT, alignmentMatchsMap, refIndexMap, queryIndexMap, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES, refMaximumTimes, queryMaximumTimes, calculateIndelDistance, false);

        for (size_t i = 0; i < alignmentMatchsMap.size(); ++i) {
            std::vector<size_t> keepIndexs;
            std::set<size_t> keepIndexset;
            for (size_t j = 0; j < alignmentMatchsMap[i].size(); ++j) {
                assert(alignmentMatchsMap[i][j].getStrand() == alignmentMatchsMap[i][0].getStrand());
                if (alignmentMatchsMap[i][j].getStrand() == POSITIVE) {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryEndPos() <
                               alignmentMatchsMap[i][j].getQueryStartPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() <
                               alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        std::cout << "line 1329" << std::endl;
                    }
                } else {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryStartPos() >
                               alignmentMatchsMap[i][j].getQueryEndPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() <
                               alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        std::cout << "line 1342" << std::endl;
                    }
                }
            }

            int this_size = alignmentMatchsMap[i].size();
            for (int jj = this_size - 1; jj >= 0; jj--) {
                if (keepIndexset.find(jj) == keepIndexset.end()) {
                    alignmentMatchsMap[i].erase(alignmentMatchsMap[i].begin() + jj);
                }
            }
            bool keepTandemDuplication = false;
            double inversion_PENALTY = -1;
            bool considerInversion = false;
            myAlignmentMatchSort(alignmentMatchsMap[i], inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion); // since here the considerInversion is false, so the parameters of inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, would not be used
        }
    }

    int is_hpc = 0; // no, do not use  homopolymer-compressed (HPC) minimizers.
    if (H) {
        is_hpc = 1;
    }

    int bucket_bits = 2;
    int n = 1;
    bool changed = false;
    if (searchForNewAnchors) {
        for (size_t i = 0; i < alignmentMatchsMap.size(); ++i) {
            std::set<int32_t> blackList;
            changed = true;
            while (alignmentMatchsMap[i].size() > 1 && changed) {
                changed = false;
                size_t startRef = alignmentMatchsMap[i][0].getRefStartPos(); // to skip the first one
                size_t startQuery;
                STRAND strand = alignmentMatchsMap[i][0].getStrand();
                size_t endRef;
                size_t endQuery;

                std::string refChr = alignmentMatchsMap[i][0].getRefChr();
                std::string queryChr = alignmentMatchsMap[i][0].getQueryChr();

                if (POSITIVE == strand) {
                    startQuery = alignmentMatchsMap[i][0].getQueryStartPos(); // to skip the first one
                    int32_t alignmentMatchsMap_i_size = alignmentMatchsMap[i].size();
                    for (int32_t alignmentMatchsMap_i_index = 0;
                         alignmentMatchsMap_i_index < alignmentMatchsMap_i_size; alignmentMatchsMap_i_index++) {
                        AlignmentMatch orthologPair = alignmentMatchsMap[i][alignmentMatchsMap_i_index];
                        if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() != startQuery) {
                            endQuery = orthologPair.getQueryStartPos() - 1;
                        } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryStartPos() == startQuery) {
                            endRef = orthologPair.getRefStartPos() - 1;
                        } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() == startQuery) {

                        } else {
                            endRef = orthologPair.getRefStartPos() - 1;
                            endQuery = orthologPair.getQueryStartPos() - 1;

                            std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                            std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                            if ((refSeq.size() * querySeq.size() > wfaSize3 * wfaSize3) && (refSeq.size() > k && querySeq.size() > k) && blackList.find(startRef) == blackList.end()) {
                                mm_idxopt_t iopt;
                                mm_mapopt_t mopt;
                                mm_verbose = 2; // disable message output to stderr
                                mm_set_opt(0, &iopt, &mopt);
                                mopt.flag |= MM_F_CIGAR; // DO NOT perform alignment
                                mopt.flag |= MM_F_NO_PRINT_2ND;
                                mopt.bw = windowWidth / 5;
                                mopt.flag |= MM_F_NO_LJOIN; // together the last one, control the maximum gap length on the local alignment region (novel seed)
                                mopt.a = matchingScore;
                                mopt.b = mismatchingPenalty;
                                mopt.q = openGapPenalty1;
                                mopt.e = extendGapPenalty1;
                                mopt.q2 = openGapPenalty1;
                                mopt.e2 = extendGapPenalty1;
                                mopt.bw = k;
                                mopt.mid_occ = 2; // ignore seeds with occurrences above this threshold
                                mopt.min_cnt = 2;// min number of minimizers on each chain
                                int32_t referenceSeqLength = refSeq.length();
                                int32_t querySeqLength = querySeq.length();
                                char *reference_seq_array = new char[referenceSeqLength + 1];
                                char *query_seq_array = new char[querySeqLength + 1];

                                strcpy(reference_seq_array, refSeq.c_str());
                                const char *refseq[1] = {reference_seq_array};
                                strcpy(query_seq_array, querySeq.c_str());

                                std::string queryName = "temp";
                                int nameLength = queryName.length();
                                // declaring character array
                                char name_array[nameLength + 1];
                                // copying the contents of the string to char array
                                strcpy(name_array, queryName.c_str());
                                const char *name[1] = {"temp"};

                                mm_idx_t *mi = mm_idx_str(w, k, is_hpc, bucket_bits, n, refseq, name);
                                mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
                                mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
                                int j, n_reg;
                                mm_reg1_t *reg = mm_map(mi, querySeqLength, query_seq_array, &n_reg, tbuf, &mopt, name_array); // get all hits for the query

                                if (n_reg > 0 && (&reg[0])->rev == 0) {
                                    mm_reg1_t *r = &reg[0];
                                    int32_t newAnchorRefEnd = r->re - 1;
                                    int32_t newAnchorRefStart = startRef + r->rs;
                                    if (blackList.find(newAnchorRefStart) == blackList.end()) {
                                        int32_t newAnchorQueryEnd = r->qe - 1;
                                        std::string alignmentName = "localAlignment_" + refChr + "_" + std::to_string(startRef + r->rs) + "_" + std::to_string(startRef + newAnchorRefEnd);
                                        double length = r->re - r->rs + 1;
                                        double numberofMs = 0;
                                        for (uint32_t i = 0; i < r->p->n_cigar; ++i) {
                                            if ("MIDNSH"[r->p->cigar[i] & 0xf] == 'M') {
                                                int32_t thisLength = r->p->cigar[i] >> 4;
                                                numberofMs = numberofMs + thisLength;
                                            }
                                        }
                                        double similarity = numberofMs / length;
                                        if (similarity > minimumSimilarity2) {
                                            double thisScore = r->score / 2;
                                            thisScore = thisScore / (abs(r->re - r->rs) * matchingScore);

                                            AlignmentMatch orthologPair(refChr, queryChr,
                                                                        startRef + r->rs, startRef + newAnchorRefEnd,
                                                                        startQuery + r->qs,
                                                                        startQuery + newAnchorQueryEnd, thisScore,
                                                                        POSITIVE,
                                                                        alignmentName, alignmentName);

                                            alignmentMatchsMap[i].push_back(orthologPair);
                                            changed = true;
                                        }
                                        blackList.insert(newAnchorRefStart);
                                    }
                                } else {
                                    blackList.insert(startRef);
                                }

                                for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                                    mm_reg1_t *r = &reg[j];
                                    free(r->p);
                                }

                                free(reg);
                                mm_tbuf_destroy(tbuf);
                                mm_idx_destroy(mi);
                                delete [] reference_seq_array;
                                delete [] query_seq_array;
                            } else {
                                blackList.insert(startRef);
                            }
                        }
                        startRef = orthologPair.getRefEndPos() + 1;
                        startQuery = orthologPair.getQueryEndPos() + 1;
                    }
                } else {
                    endQuery = alignmentMatchsMap[i][0].getQueryEndPos(); // skip the first one
                    int32_t size_i = alignmentMatchsMap[i].size();
                    for (int32_t i_index = 0; i_index < size_i; i_index++) {
                        AlignmentMatch orthologPair = alignmentMatchsMap[i][i_index];

                        if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() != endQuery) {
                            startQuery = orthologPair.getQueryEndPos() + 1;
                        } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryEndPos() == endQuery) {
                            endRef = orthologPair.getRefStartPos() - 1;
                        } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() == endQuery) {

                        } else {
                            endRef = orthologPair.getRefStartPos() - 1;
                            startQuery = orthologPair.getQueryEndPos() + 1;

                            std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                            std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);

                            if ((refSeq.size() * querySeq.size() > wfaSize3 * wfaSize3) && (refSeq.size() > k && querySeq.size() > k) && blackList.find(startRef) == blackList.end()) {
                                mm_idxopt_t iopt;
                                mm_mapopt_t mopt;

                                mm_verbose = 2; // disable message output to stderr
                                mm_set_opt(0, &iopt, &mopt);
                                //mopt.flag &= ~ MM_F_CIGAR; // DO NOT perform alignment
                                mopt.flag |= MM_F_CIGAR; // DO NOT perform alignment
                                mopt.flag |= MM_F_NO_PRINT_2ND;
                                mopt.bw = windowWidth / 5;
                                mopt.flag |= MM_F_NO_LJOIN; // together the last one, control the maximum gap length on the local alignment region (novel seed)
                                mopt.a = matchingScore;
                                mopt.b = mismatchingPenalty;
                                mopt.q = openGapPenalty1;
                                mopt.e = extendGapPenalty1;
                                mopt.q2 = openGapPenalty1;
                                mopt.e2 = extendGapPenalty1;
                                mopt.bw = k;
                                mopt.mid_occ = 2; // ignore seeds with occurrences above this threshold
                                mopt.min_cnt = 2;// min number of minimizers on each chain
                                int32_t referenceSeqLength = refSeq.length();
                                int32_t querySeqLength = querySeq.length();
                                char *reference_seq_array = new char[referenceSeqLength + 1];
                                char *query_seq_array = new char[querySeqLength + 1];

                                strcpy(reference_seq_array, refSeq.c_str());
                                const char *refseq[1] = {reference_seq_array};

                                strcpy(query_seq_array, querySeq.c_str());

                                std::string queryName = "temp";
                                int nameLength = queryName.length();
                                // declaring character array
                                char name_array[nameLength + 1];
                                // copying the contents of the string to char array
                                strcpy(name_array, queryName.c_str());
                                const char *name[1] = {"temp"};

                                mm_idx_t *mi = mm_idx_str(w, k, is_hpc, bucket_bits, n, refseq, name);
                                mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
                                mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
                                mm_reg1_t *reg;
                                int j, n_reg;
                                reg = mm_map(mi, querySeqLength, query_seq_array, &n_reg, tbuf, &mopt, name_array); // get all hits for the query
                                if (n_reg > 0 && (&reg[0])->rev == 0) {
                                    mm_reg1_t *r = &reg[0];
                                    int32_t newAnchorRefEnd = r->re - 1;
                                    int32_t newAnchorRefStart = startRef + r->rs;
                                    if (blackList.find(newAnchorRefStart) == blackList.end()) {
                                        double length = r->re - r->rs + 1;
                                        double numberofMs = 0;
                                        for (uint32_t i = 0; i < r->p->n_cigar; ++i) {
                                            if ("MIDNSH"[r->p->cigar[i] & 0xf] == 'M') {
                                                int32_t thisLength = r->p->cigar[i] >> 4;
                                                numberofMs = numberofMs + thisLength;
                                            }
                                        }

                                        double similarity = numberofMs / length;
                                        if (similarity > minimumSimilarity2) {
                                            int32_t newAnchorQueryEnd = r->qe - 1;
                                            std::string alignmentName = "localAlignment_" + refChr + "_" + std::to_string(startRef + r->rs) + "_" + std::to_string(startRef + newAnchorRefEnd);
                                            double thisScore = r->score / 2;
                                            thisScore = thisScore / (abs(r->re - r->rs + 1) * matchingScore);
                                            AlignmentMatch orthologPair(refChr, queryChr,
                                                                        startRef + r->rs, startRef + newAnchorRefEnd,
                                                                        endQuery - newAnchorQueryEnd, endQuery - r->qs,
                                                                        thisScore, NEGATIVE, alignmentName, alignmentName);
                                            alignmentMatchsMap[i].push_back(orthologPair);
                                            changed = true;
                                        }
                                        blackList.insert(newAnchorRefStart);
                                    }
                                } else {
                                    blackList.insert(startRef);
                                }

                                for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                                    mm_reg1_t *r = &reg[j];
                                    free(r->p);
                                }

                                free(reg);
                                mm_tbuf_destroy(tbuf);
                                mm_idx_destroy(mi);
                                delete [] reference_seq_array;
                                delete [] query_seq_array;
                            } else {
                                blackList.insert(startRef);
                            }
                        }
                        startRef = orthologPair.getRefEndPos() + 1;
                        endQuery = orthologPair.getQueryStartPos() - 1;
                    }
                }
                bool keepTandemDuplication = false;
                double inversion_PENALTY = -1;
                bool considerInversion = false;
                myAlignmentMatchSort(alignmentMatchsMap[i], inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);
            }
        }
    }

    {
        for (size_t i = 0; i < alignmentMatchsMap.size(); ++i) {
            std::vector<size_t> keepIndexs;
            std::set<size_t> keepIndexset;
            for (size_t j = 0; j < alignmentMatchsMap[i].size(); ++j) {
                if (alignmentMatchsMap[i][j].getStrand() == POSITIVE) {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryEndPos() < alignmentMatchsMap[i][j].getQueryStartPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() < alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        // here should never run.
                        std::cout << "removed" << alignmentMatchsMap[i][j].getRefChr() << "\t"
                                  << alignmentMatchsMap[i][j].getRefStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getRefEndPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryChr() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryEndPos() << "\t"
                                  << "+" << "\t"
                                  << alignmentMatchsMap[i][j].getReferenceGeneName() << "\t"
                                  << alignmentMatchsMap[i][j].getScore() << std::endl;
                    }
                } else {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryStartPos() > alignmentMatchsMap[i][j].getQueryEndPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() < alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        std::cout << "removed" << alignmentMatchsMap[i][j].getRefChr() << "\t"
                                  << alignmentMatchsMap[i][j].getRefStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getRefEndPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryChr() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryEndPos() << "\t"
                                  << "-" << "\t"
                                  << alignmentMatchsMap[i][j].getReferenceGeneName() << "\t"
                                  << alignmentMatchsMap[i][j].getScore() << std::endl;
                    }
                }
            }

            int this_size = alignmentMatchsMap[i].size();
            for (int jj = this_size - 1; jj >= 0; jj--) {
                if (keepIndexset.find(jj) == keepIndexset.end()) {
                    alignmentMatchsMap[i].erase(alignmentMatchsMap[i].begin() + jj);
                }
            }

            bool keepTandemDuplication = false;
            double inversion_PENALTY = -1;
            bool considerInversion = false;

            // since here the considerInversion is false, so the parameters of inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, would not be used
            myAlignmentMatchSort(alignmentMatchsMap[i], inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);
        }
    }
}
