//
// Created by Baoxing song on 20.10.18.
//

#include "deNovoGenomeVariantCalling.h"

std::mutex g_num_mutex;

void genomeAlignmentSingleThread(std::vector<AlignmentMatch> alignmentMatchs,
                                 const bool outPutMaf, const bool outPutFraged,
                                 std::ofstream &omaffile, std::ofstream &ofragfile,
                                 std::string refChr, std::string queryChr,
                                 std::string refSequence, std::string targetSequence,
                                 const int chrWidth, const std::string refFileName, const std::string queryFileName,
                                 const int32_t widownWidth, const int32_t wfaSize, const int32_t wfaSize2,
                                 const int32_t matchingScore, const int32_t mismatchingPenalty,
                                 const int32_t openGapPenalty1, const int32_t extendGapPenalty1,
                                 const int32_t openGapPenalty2, const int32_t extendGapPenalty2,
                                 const int32_t min_wavefront_length, const int32_t max_distance_threshold,
                                 std::atomic_int &number_of_runing_threads, std::map<std::string, std::string> parameters) {

    Scorei m(matchingScore, mismatchingPenalty);


    bool checkResult = false;

    STRAND strand = alignmentMatchs[0].getStrand();

//    std::cout << "line 39" << std::endl;
    if (POSITIVE == strand) {
        size_t startRef = alignmentMatchs[0].getRefStartPos();
        size_t startQuery;
        size_t endRef;
        size_t endQuery;
        std::stringstream refAlign;
        std::stringstream queryAlign;
        int64_t alignmentScore = 0;
        startQuery = alignmentMatchs[0].getQueryStartPos();
        for (AlignmentMatch orthologPair: alignmentMatchs) {
            if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() != startQuery) {
                endQuery = orthologPair.getQueryStartPos() - 1;
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);

                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.length(), '-');
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * querySeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * querySeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q
                              << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryStartPos() == startQuery) {
                endRef = orthologPair.getRefStartPos() - 1;
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.length(), '-');
                std::string _alignment_d = refSeq;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * refSeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * refSeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() == startQuery) {

            } else {
                endRef = orthologPair.getRefStartPos() - 1;
                endQuery = orthologPair.getQueryStartPos() - 1;
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);
                {

                    std::string _alignment_q;
                    std::string _alignment_d;
//                    std::cout << "line 78" << std::endl;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                          m, parameters);
                    //                  std::cout << " line 80" << std::endl;
                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;

                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());

                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                            m);
                            tempd = _alignment_d;
                            tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                            tempq = _alignment_q;
                            tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                        }
                        assert(tempd.compare(refSeq) == 0);
                        assert(tempq.compare(querySeq) == 0);
                    }
                    alignmentScore += thiScore;

                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q
                                  << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            {
                startRef = orthologPair.getRefStartPos();
                startQuery = orthologPair.getQueryStartPos();
                endRef = orthologPair.getRefEndPos();
                endQuery = orthologPair.getQueryEndPos();
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;
//                    std::cout << "line 126" << std::endl;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                          m, parameters);
                    //                  std::cout << "line 128" << std::endl;
                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                            m);
                            tempd = _alignment_d;
                            tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                            tempq = _alignment_q;
                            tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        }
                        assert(tempd.compare(refSeq) == 0);
                        assert(tempq.compare(querySeq) == 0);
                    }
                    alignmentScore += thiScore;

                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q
                                  << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            startRef = orthologPair.getRefEndPos() + 1;
            startQuery = orthologPair.getQueryEndPos() + 1;
        }
//        std::cout << "line 153" << std::endl;
        {
            std::string refGenomerSequence = getSubsequence(refSequence, alignmentMatchs[0].getRefStartPos(), alignmentMatchs[alignmentMatchs.size() - 1].getRefEndPos());
            std::string queryGenomerSequence = getSubsequence(targetSequence, alignmentMatchs[0].getQueryStartPos(), alignmentMatchs[alignmentMatchs.size() - 1].getQueryEndPos());

            //if (checkResult) {
            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            assert(temp.compare(refGenomerSequence) == 0);

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            assert(temp.compare(queryGenomerSequence) == 0);
            //}
            if (outPutMaf) {
                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getRefStartPos() - 1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << refSequence.size() << "\t"
                         << refAlign.str() << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getQueryStartPos() - 1 << "\t" << std::setw(9) << queryGenomerSequence.size() << "\t+\t"
                         << targetSequence.size() << "\t" << queryAlign.str() << std::endl
                         << std::endl;
                g_num_mutex.unlock();
            }
        }
    } else {
//        std::cout << "line 183" << std::endl;
        size_t startRef = alignmentMatchs[0].getRefStartPos();
        size_t startQuery;
        size_t endRef;
        size_t endQuery = alignmentMatchs[0].getQueryEndPos();
        std::stringstream refAlign;
        std::stringstream queryAlign;
        std::string refChr = alignmentMatchs[0].getRefChr();
        std::string queryChr = alignmentMatchs[0].getQueryChr();

        int64_t alignmentScore = 0;
        for (AlignmentMatch orthologPair: alignmentMatchs) {

            if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() != endQuery) {
                startQuery = orthologPair.getQueryEndPos() + 1;
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery, strand);
                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.length(), '-');
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * querySeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * querySeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequence.size() << "\t" << _alignment_q
                              << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryEndPos() == endQuery) {
                endRef = orthologPair.getRefStartPos() - 1;
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.length(), '-');
                std::string _alignment_d = refSeq;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * refSeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * refSeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t-\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() == endQuery) {

            } else {
                endRef = orthologPair.getRefStartPos() - 1;
                startQuery = orthologPair.getQueryEndPos() + 1;
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery, strand);
                {

                    std::string _alignment_q;
                    std::string _alignment_d;
//                    std::cout << "line 288" << std::endl;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                          m, parameters);
//                    std::cout << "line 230" << std::endl;

                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                            m);
                            tempd = _alignment_d;
                            tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                            tempq = _alignment_q;
                            tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                        }
                        assert(tempd.compare(refSeq) == 0);
                        assert(tempq.compare(querySeq) == 0);
                    }
                    alignmentScore += thiScore;

                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequence.size() << "\t" << _alignment_q
                                  << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            {
                startRef = orthologPair.getRefStartPos();
                endQuery = orthologPair.getQueryEndPos();
                endRef = orthologPair.getRefEndPos();
                startQuery = orthologPair.getQueryStartPos();
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery, strand);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;
//                    std::cout << refSeq << std::endl << querySeq << std::endl << "line 276" << std::endl;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                          m, parameters);

                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold,
                                                            m);
                            tempd = _alignment_d;
                            tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                            tempq = _alignment_q;
                            tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                        }
                        assert(tempd.compare(refSeq) == 0);
                        assert(tempq.compare(querySeq) == 0);
                    }
                    alignmentScore += thiScore;

                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequence.size() << "\t" << _alignment_q
                                  << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            startRef = orthologPair.getRefEndPos() + 1;
            endQuery = orthologPair.getQueryStartPos() - 1;
        }
        {
            std::string refGenomerSequence = getSubsequence(refSequence, alignmentMatchs[0].getRefStartPos(), alignmentMatchs[alignmentMatchs.size() - 1].getRefEndPos());
            std::string queryGenomerSequence = getSubsequence(targetSequence, alignmentMatchs[0].getQueryEndPos(), alignmentMatchs[alignmentMatchs.size() - 1].getQueryStartPos(), strand);
            //if (checkResult) {
            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            assert(temp.compare(refGenomerSequence) == 0);

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            assert(temp.compare(queryGenomerSequence) == 0);
            //}
            if (outPutMaf) {
                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getRefStartPos() - 1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << refSequence.size() << "\t"
                         << refAlign.str() << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) <</* queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << alignmentMatchs[alignmentMatchs.size() - 1].getQueryStartPos() - 1 << "\t" << std::setw(9) << queryGenomerSequence.size()
                         << "\t-\t" << targetSequence.size() << "\t" << queryAlign.str() << std::endl
                         << std::endl;
                g_num_mutex.unlock();
            }
        }
    }
    number_of_runing_threads = number_of_runing_threads - 1;
//    std::cout << "line 336" << std::endl;
}


void genomeAlignment(std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                     const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                     const int32_t &widownWidth, const int32_t &wfaSize, const int32_t &wfaSize2,
                     const std::string &outPutMafFile, const std::string &outPutFragedFile, /*std::string & outPutLocalalignmentFile,*/
                     const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                     const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1,
                     const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                     int32_t &seed_window_size, const int32_t &mini_cns_score, const int32_t &step_size,
                     const int32_t &matrix_boundary_distance, const int32_t &scoreThreshold, const int32_t &w, const int32_t &xDrop,
                     const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, const int &maxThread, std::map<std::string, std::string> &parameters) {


    bool outPutMaf = false;
    bool outPutFraged = false;
//    bool outPutLocalalignment = false;

    if (outPutMafFile.size() > 0) {
        outPutMaf = true;
    }
    if (outPutFragedFile.size() > 0) {
        outPutFraged = true;
    }
//    if ( outPutLocalalignmentFile.size() > 0 ){
//        outPutLocalalignment = true;
//    }
//    std::cout << "line 354" << std::endl;
    std::map<std::string, std::string> refSequences;
    readFastaFile(refFastaFilePath, refSequences);

    std::map<std::string, std::string> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);
    //  std::cout << "line 360" << std::endl;
    int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    queryFileName = elems.back();
//    std::cout << "line 370" << std::endl;
    for (std::map<std::string, std::string>::iterator itchr = refSequences.begin(); itchr != refSequences.end(); ++itchr) {
        if ((refFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }
    for (std::map<std::string, std::string>::iterator itchr = targetSequences.begin(); itchr != targetSequences.end(); ++itchr) {
        if ((queryFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }
//    std::cout << "line 381" << std::endl;


    std::ofstream omaffile;
    std::ofstream ofragfile;
    std::ofstream oLocalalignment;
    if (outPutMaf) {
        omaffile.open(outPutMafFile);
        omaffile << "##maf version=1" << std::endl;
    }

    if (outPutFraged) {
        ofragfile.open(outPutFragedFile);
        ofragfile << "##maf version=1" << std::endl;
    }

    SequenceCharToUInt8 sequenceCharToUInt8;
    int32_t sizei = alignmentMatchsMap.size();


    std::atomic_int number_of_runing_threads(0);



//    std::cout << "line 404" << std::endl;
    for (int32_t i = sizei - 1; i >= 0; --i) {
        std::vector<AlignmentMatch> alignmentMatchs = alignmentMatchsMap[i];
        //      std::cout << "line 407" << std::endl;
        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if (number_of_runing_threads < maxThread) {
                std::string refChr = alignmentMatchs[0].getRefChr();
                std::string queryChr = alignmentMatchs[0].getQueryChr();
                std::thread t(genomeAlignmentSingleThread, alignmentMatchs, outPutMaf, outPutFraged,
                              std::ref(omaffile), std::ref(ofragfile), refChr, queryChr, refSequences[refChr], targetSequences[queryChr],
                              chrWidth, refFileName, queryFileName,
                              widownWidth, wfaSize, wfaSize2, matchingScore, mismatchingPenalty,
                              openGapPenalty1, extendGapPenalty1,
                              openGapPenalty2, extendGapPenalty2,
                              min_wavefront_length, max_distance_threshold, std::ref(number_of_runing_threads), parameters);
                ++number_of_runing_threads;
                t.detach();
                isThisThreadUnrun = false;
                break;
            } else {
                usleep(1000);
            }
        }
    }
    while (number_of_runing_threads > 0) {// wait for all the thread
        usleep(1000);
    }
    if (outPutMaf) {
        omaffile.close();
    }
    if (outPutFraged) {
        ofragfile.close();
    }
}


void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::ofstream &ovcffile, std::string chr, std::string &refSequence) {
    alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, chr, refSequence, 0);
}

void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::ofstream &ovcffile, std::string chr, std::map<std::string, std::string> &refSequences) {
    alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, chr, refSequences, 0);
}

void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::ofstream &ovcffile, std::string chr, std::map<std::string, std::string> &refSequences, int32_t refLetterNumber) {
    alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, chr, refSequences[chr], refLetterNumber);
}

void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::ofstream &ovcffile, std::string chr, std::string &refSequence, int32_t refLetterNumber) {
    std::vector<Variant> sdiRecordsThisOne;
    alignmentToVcf(queryAlignSeq, refAlignSeq, sdiRecordsThisOne, chr, refSequence, refLetterNumber);
    g_num_mutex.lock();
    for (std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin(); itVariant != sdiRecordsThisOne.end(); ++itVariant) {
        ovcffile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChromosome() << "_" << itVariant->getPosition() << "\t" + itVariant->getReference() << "\t" << itVariant->getAlternative() << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1:50:1"
                 << std::endl;
    }
    g_num_mutex.unlock();
}


void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::ofstream &ovcffile, std::string chr, std::map<std::string, std::string> &refSequences, int32_t refLetterNumber, const bool &gvcf) {
    alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, chr, refSequences[chr], refLetterNumber, gvcf);
}

void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::ofstream &ovcffile, std::string chr, std::string &refSequence, int32_t refLetterNumber, const bool &gvcf) {
    std::vector<Variant> sdiRecordsThisOne;
    alignmentToVcf(queryAlignSeq, refAlignSeq, sdiRecordsThisOne, chr, refSequence, refLetterNumber);
    if (gvcf) {
        g_num_mutex.lock();
        int32_t lastPosition = refLetterNumber + 1;
        for (std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin(); itVariant != sdiRecordsThisOne.end(); ++itVariant) {
            if (lastPosition < itVariant->getPosition()) {
                ovcffile << itVariant->getChromosome() << "\t" << lastPosition << "\t" << itVariant->getChromosome() << "_" << lastPosition << "\t" << refSequence[lastPosition - 1] << "\t" << "<NON_REF>" << "\t.\t.\tEND=" << std::to_string(itVariant->getPosition() - 1)
                         << "\tGT:AD:DP\t0:1,0:1" << std::endl;
            }
            ovcffile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChromosome() << "_" << itVariant->getPosition() << "\t" + itVariant->getReference() << "\t" << itVariant->getAlternative() << ",<NON_REF>"
                     << "\t.\t.\t.\tGT:AD:DP\t1:0,1,0:1" << std::endl;
            lastPosition = itVariant->getPosition() + itVariant->getReference().length();
        }
        g_num_mutex.unlock();
    } else {
        g_num_mutex.lock();
        for (std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin(); itVariant != sdiRecordsThisOne.end(); ++itVariant) {
            ovcffile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChromosome() << "_" << itVariant->getPosition() << "\t" + itVariant->getReference() << "\t" << itVariant->getAlternative() << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1:50:1"
                     << std::endl;
        }
        g_num_mutex.unlock();
    }
}

void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::vector<Variant> &sdiRecordsThisOne, std::string chr, std::map<std::string, std::string> &refSequences, int32_t refLetterNumber) {
    alignmentToVcf(queryAlignSeq, refAlignSeq, sdiRecordsThisOne, chr, refSequences[chr], refLetterNumber);
}

bool to_T(char c) { return c == 'U'; }

bool to_N(char c) {
    return c == 'R' || c == 'Y' || c == 'S' || c == 'W' || c == 'K' || c == 'M' || c == 'B' || c == 'D' || c == 'H' || c == 'V';
}

void alignmentToVcf(std::string &queryAlignSeq, std::string &refAlignSeq, std::vector<Variant> &sdiRecordsThisOne, std::string chr, std::string &refSequence, int32_t refLetterNumber) {

    std::replace_if(refAlignSeq.begin(), refAlignSeq.end(), to_T, 'T');
    std::replace_if(refAlignSeq.begin(), refAlignSeq.end(), to_N, 'N');

    std::replace_if(queryAlignSeq.begin(), queryAlignSeq.end(), to_T, 'T');
    std::replace_if(queryAlignSeq.begin(), queryAlignSeq.end(), to_N, 'N');

    FirstLastList sdiRecords;
    std::cout << "reference length:" << std::to_string(refAlignSeq.length()) << std::endl;
    for (int ai = 0; ai < refAlignSeq.length(); ai++) {
        if (refAlignSeq[ai] != '-') {
            ++refLetterNumber;
        }
//        std::cout << "ai:" << std::to_string(ai) << std::endl;
        if (refAlignSeq[ai] != queryAlignSeq[ai]) {
            if (queryAlignSeq[ai] == '-') {
                if (NULL != sdiRecords.getLast() &&
                    sdiRecords.getLast()->getMapSingleRecord().getChanginglength() < 0 &&
                    sdiRecords.getLast()->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                    sdiRecords.getLast()->getMapSingleRecord().getPosition() - sdiRecords.getLast()->getMapSingleRecord().getChanginglength() == refLetterNumber) {
                    int position = sdiRecords.getLast()->getMapSingleRecord().getPosition();
                    std::string ori = sdiRecords.getLast()->getMapSingleRecord().getReference() + refAlignSeq[ai];
                    std::string result = "-";
//                    std::cout << "liner 581" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                    Variant mapSingleRecord = Variant(chr, position, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.deleteLast();
                    sdiRecords.insertLast(data);

                } else {
                    std::string ori(1, refAlignSeq[ai]);
                    std::string result = "-";
//                    std::cout << "liner 590" << chr <<  std::to_string(refLetterNumber) <<  ori << result << std::endl;
                    Variant mapSingleRecord = Variant(chr, refLetterNumber, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                }
            } else if (refAlignSeq[ai] == '-') {
                if (sdiRecords.getLast() == NULL) {
                    int position = refLetterNumber + 1;
                    std::string ori = "-";
                    std::string result(1, queryAlignSeq[ai]);
//                    std::cout << "liner 600" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                    Variant mapSingleRecord = Variant(chr, position, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                } else {
                    if (sdiRecords.getLast()->getMapSingleRecord().getPosition() ==
                        (refLetterNumber + 1)
                        && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                        sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                        int position = refLetterNumber + 1;
                        std::string ori = "-";
                        std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + queryAlignSeq[ai];
//                        std::cout << "liner 614" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                        Variant mapSingleRecord = Variant(chr, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.deleteLast();
                        sdiRecords.insertLast(data);
                    } else {
                        int position = refLetterNumber + 1;
                        std::string ori = "-";
                        std::string result(1, queryAlignSeq[ai]);
//                        std::cout << "liner 623" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                        Variant mapSingleRecord = Variant(chr, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
                    }
                }
            } else {
                int position = refLetterNumber;
                std::string ori(1, refAlignSeq[ai]);
                std::string result(1, queryAlignSeq[ai]);
//                std::cout << "liner 633" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                Variant mapSingleRecord = Variant(chr, position, ori, result);
                Data *data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            }
        }
    }

    std::cout << "initial variants calling done" << std::endl;

    for (int runingCound = 0; runingCound < 2; ++runingCound) {
        if ((sdiRecords.getFirst() != NULL)
            && (sdiRecords.getFirst()->getNext() != NULL)) {
            Data *prevOne = sdiRecords.getFirst();
            Data *currOne = (sdiRecords.getFirst()->getNext());
            while ((currOne != NULL) && (NULL != currOne->getNext())) {
//                std::cout << "647" << std::endl;
                if (sdiRecords.getFirst() == currOne) {
                    prevOne = currOne;
                    currOne = prevOne->getNext();
                }
                if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                    prevOne->getMapSingleRecord().getChanginglength() < 0 &&
                    (prevOne->getMapSingleRecord().getPosition() + abs(prevOne->getMapSingleRecord().getChanginglength())) ==
                    currOne->getMapSingleRecord().getPosition()
                    && prevOne->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                    currOne->getMapSingleRecord().getAlternative().compare("-") == 0) { // merge two deletions

                    int position = prevOne->getMapSingleRecord().getPosition();
                    std::string ori = prevOne->getMapSingleRecord().getReference() + currOne->getMapSingleRecord().getReference();
                    std::string result = "-";
//                    std::cout << "liner 662" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                    Variant mapSingleRecord2(chr, position, ori, result);
                    //delete prev one begin
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                    } else {
                        prevOne->getPrev()->setNext(currOne);
                        currOne->setPrev(prevOne->getPrev());
                        delete (prevOne);
                    } //delete prev one end
                    currOne->setMapSingleRecord(mapSingleRecord2);
                    prevOne = currOne->getPrev();
                    if (prevOne == NULL) {
                        prevOne = currOne;
                        currOne = prevOne->getNext();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() == 0 &&
                           0 == currOne->getMapSingleRecord().getReference().compare(currOne->getMapSingleRecord().getAlternative())) { // nonsense records
                    //delete current one
//                    std::cout << "820 delete prev" << std::endl;
                    prevOne->setNext(currOne->getNext());
                    currOne->getNext()->setPrev(prevOne);
                    delete (currOne);
                    currOne = prevOne->getNext();
                } else if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                           prevOne->getMapSingleRecord().getChanginglength() > 0
                           && currOne->getMapSingleRecord().getReference().compare(
                        prevOne->getMapSingleRecord().getAlternative()) == 0 &&
                           currOne->getMapSingleRecord().getPosition() ==
                           prevOne->getMapSingleRecord().getPosition()) { //delete one insertion and next reverse sence deletion
                    //delete current one and prev
//                    std::cout << "830 delete current one and prev" << std::endl;
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                        sdiRecords.deleteFirst();
                        prevOne = sdiRecords.getFirst();
                        currOne = (sdiRecords.getFirst()->getNext());
                    } else if (currOne == sdiRecords.getLast()) {
                        sdiRecords.deleteLast();
                        sdiRecords.deleteLast();
                        currOne = sdiRecords.getLast();
                        prevOne = currOne->getPrev();
                    } else {
                        currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                        currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                        Data *temp = currOne->getNext();
                        delete (currOne->getPrev());
                        delete (currOne);
                        currOne = temp;
                        prevOne = temp->getPrev();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() > 0 &&
                           prevOne->getMapSingleRecord().getChanginglength() < 0
                           && currOne->getMapSingleRecord().getAlternative().compare(
                        prevOne->getMapSingleRecord().getReference()) == 0 &&
                           (currOne->getMapSingleRecord().getPosition() - 1) ==
                           prevOne->getMapSingleRecord().getPosition()) {
                    //delete current one and prev
//                    std::cout << "850 delete current one and prev" << std::endl;
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                        sdiRecords.deleteFirst();
                        prevOne = sdiRecords.getFirst();
                        currOne = (sdiRecords.getFirst()->getNext());
                    } else {
                        currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                        currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                        Data *temp = currOne->getNext();
                        delete (currOne->getPrev());
                        delete (currOne);
                        currOne = temp;
                        prevOne = temp->getPrev();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() < 0 && currOne->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                           ( //(sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition()-1)  && sdiRecordsThisOne[j - 1].getReference() != "-" ) ||
                                   (currOne->getPrev()->getMapSingleRecord().getReference()[0] != '-' &&
                                    currOne->getPrev()->getMapSingleRecord().getPosition() + currOne->getPrev()->getMapSingleRecord().getReference().size() == currOne->getMapSingleRecord().getPosition()))) {
                    int position = currOne->getPrev()->getMapSingleRecord().getPosition();
                    std::string ori = currOne->getPrev()->getMapSingleRecord().getReference() + currOne->getMapSingleRecord().getReference();
                    std::string result = currOne->getPrev()->getMapSingleRecord().getAlternative() + currOne->getMapSingleRecord().getAlternative();

                    ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                    result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 846" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                    Variant mapSingleRecord2(chr, position, ori, result);

                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                    } else {
                        prevOne->getPrev()->setNext(currOne);
                        currOne->setPrev(prevOne->getPrev());
                        delete (prevOne);
                    } //delete prev one end
                    currOne->setMapSingleRecord(mapSingleRecord2);
                    prevOne = currOne->getPrev();
                    if (prevOne == NULL) {
                        prevOne = currOne;
                        currOne = prevOne->getNext();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                           currOne->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                           (currOne->getPrev()->getMapSingleRecord().getReference() == "-" && currOne->getPrev()->getMapSingleRecord().getPosition() == (currOne->getMapSingleRecord().getPosition()))
                        ) {
                    int position = currOne->getPrev()->getMapSingleRecord().getPosition();
                    std::string ori = currOne->getPrev()->getMapSingleRecord().getReference() + currOne->getMapSingleRecord().getReference();
                    std::string result = currOne->getPrev()->getMapSingleRecord().getAlternative() + currOne->getMapSingleRecord().getAlternative();

                    ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                    result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 873" << chr <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                    Variant mapSingleRecord2(chr, position, ori, result);
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                    } else {
                        prevOne->getPrev()->setNext(currOne);
                        currOne->setPrev(prevOne->getPrev());
                        delete (prevOne);
                    } //delete prev one end
                    currOne->setMapSingleRecord(mapSingleRecord2);
                    prevOne = currOne->getPrev();
                    if (prevOne == NULL) {
                        prevOne = currOne;
                        currOne = prevOne->getNext();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() > 0 &&
                           currOne->getMapSingleRecord().getReference().compare("-") == 0 &&
                           currOne->getPrev()->getMapSingleRecord().getReference().compare("-") != 0 &&
                           (currOne->getPrev()->getMapSingleRecord().getPosition() + currOne->getPrev()->getMapSingleRecord().getReference().size()) == (currOne->getMapSingleRecord().getPosition())) {
                    int position = currOne->getPrev()->getMapSingleRecord().getPosition();
                    std::string ori = currOne->getPrev()->getMapSingleRecord().getReference() + currOne->getMapSingleRecord().getReference();
                    std::string result = currOne->getPrev()->getMapSingleRecord().getAlternative() + currOne->getMapSingleRecord().getAlternative();
                    ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                    result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 888" << chr <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                    Variant mapSingleRecord2(chr, position, ori, result);
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                    } else {
                        prevOne->getPrev()->setNext(currOne);
                        currOne->setPrev(prevOne->getPrev());
                        delete (prevOne);
                    } //delete prev one end
                    currOne->setMapSingleRecord(mapSingleRecord2);
                    prevOne = currOne->getPrev();
                    if (prevOne == NULL) {
                        prevOne = currOne;
                        currOne = prevOne->getNext();
                    }
                } else {
                    prevOne = currOne;
                    currOne = prevOne->getNext();
                }//std::cout <<  (*itName) << ": link data structure end " << currOne->getMapSingleRecord().getPosition() << std::endl;
//                std::cout << "820" << std::endl;
            }
        }
    }
    std::cout << "variants calling merging" << std::endl;

    if (sdiRecords.getFirst() != NULL) {
        Data *thisone = sdiRecords.getFirst();
        while (thisone != NULL) {
            sdiRecordsThisOne.push_back(thisone->getMapSingleRecord());
            thisone = (thisone->getNext());
        }
    }

    // clear RAM assigned by new Data() begin
    if (sdiRecords.getFirst() != NULL) {
        Data *thisone = sdiRecords.getFirst();
        while (thisone != NULL) {
            Data *tempData = thisone;
            thisone = (thisone->getNext());
            delete (tempData);
        }
    }// clear RAM assigned by new Data() end
    std::cout << sdiRecordsThisOne.size() << std::endl;
    std::cout << "variants calling to vector" << std::endl;


    // transform link to vector and sort and merge nearby records begin
    bool ifChanged = true;
    while (ifChanged) { // link deletions together
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                sdiRecordsThisOne[j - 1].getChanginglength() < 0 &&
                (sdiRecordsThisOne[j - 1].getPosition() + abs(sdiRecordsThisOne[j - 1].getChanginglength())) == sdiRecordsThisOne[j].getPosition() &&
                sdiRecordsThisOne[j - 1].getAlternative().compare("-") == 0 &&
                sdiRecordsThisOne[j].getAlternative().compare("-") == 0) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = "-";
//                std::cout << "liner 778" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }
    // transform link to vector and sort and merge nearby records end
    std::cout << sdiRecordsThisOne.size() << std::endl;
    std::cout << "merge deletions together" << std::endl;

    ifChanged = true;
    while (ifChanged && sdiRecordsThisOne.size() > 1) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j - 1].getChanginglength() < 0 && sdiRecordsThisOne[j - 1].getAlternative().compare("-") == 0 &&
                sdiRecordsThisOne[j].getChanginglength() > 0 && sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                ((sdiRecordsThisOne[j - 1].getPosition() + sdiRecordsThisOne[j - 1].getReference().size()) == sdiRecordsThisOne[j].getPosition())) {
                int position = sdiRecordsThisOne[j - 1].getPosition();

                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 812" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }
    std::cout << sdiRecordsThisOne.size() << std::endl;
    std::cout << "merge deletion insertions together" << std::endl;

    // merge nearby indels begin
    ifChanged = true;
    while (ifChanged && sdiRecordsThisOne.size() > 1) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j].getChanginglength() < 0 && sdiRecordsThisOne[j].getAlternative().compare("-") == 0 &&
                ( //(sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition()-1)  && sdiRecordsThisOne[j - 1].getReference() != "-" ) ||
                        (sdiRecordsThisOne[j - 1].getReference()[0] != '-' &&
                         sdiRecordsThisOne[j - 1].getPosition() + sdiRecordsThisOne[j - 1].getReference().size() == sdiRecordsThisOne[j].getPosition()))) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 846" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                       sdiRecordsThisOne[j].getAlternative().compare("-") == 0 && (
                               (sdiRecordsThisOne[j - 1].getReference() != "-" && (sdiRecordsThisOne[j - 1].getPosition() + (sdiRecordsThisOne[j - 1].getReference().size())) < (sdiRecordsThisOne[j].getPosition()))
                               || (sdiRecordsThisOne[j - 1].getReference() == "-" && sdiRecordsThisOne[j - 1].getPosition() != (sdiRecordsThisOne[j].getPosition())))
                    ) {
                int position = sdiRecordsThisOne[j].getPosition() - 1;
//                std::cout << std::to_string(position-1) << chr << "\t" <<  std::to_string( refSequence[position-1]) <<  refSequence.substr(position-1, 20) << std::endl;
                std::string ori(1, refSequence[position - 1]);
                std::string result = ori;
//                std::cout << ori << " liner 860" << chr << "\t" <<  std::to_string(sdiRecordsThisOne[j].getPosition()) << "\t" <<  sdiRecordsThisOne[j].getReference() << "\t" << sdiRecordsThisOne[j].getAlternative() << std::endl;
                ori = ori + sdiRecordsThisOne[j].getReference();
//                std::cout << "liner 861" << chr << "\t" <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                       sdiRecordsThisOne[j].getAlternative().compare("-") == 0 && (sdiRecordsThisOne[j - 1].getReference() == "-" && sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition()))
                    ) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 873" << chr <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                       sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                       sdiRecordsThisOne[j - 1].getReference().compare("-") != 0 &&
                       (sdiRecordsThisOne[j - 1].getPosition() + sdiRecordsThisOne[j - 1].getReference().size()) == (sdiRecordsThisOne[j].getPosition())) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();
                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 888" << chr <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                       sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                       sdiRecordsThisOne[j - 1].getReference() != "-" &&
                       (sdiRecordsThisOne[j - 1].getPosition() + (sdiRecordsThisOne[j - 1].getReference().size())) < (sdiRecordsThisOne[j].getPosition())) {
                int position = sdiRecordsThisOne[j].getPosition() - 1;
                std::string ori(1, refSequence[position - 1]);
                std::string result = ori + sdiRecordsThisOne[j].getAlternative();
//                std::cout << "liner 901" << chr <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                       sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                       sdiRecordsThisOne[j - 1].getReference() == "-") {
                int position = sdiRecordsThisOne[j].getPosition() - 1;
                std::string ori(1, refSequence[position - 1]);
                std::string result = ori + sdiRecordsThisOne[j].getAlternative();
//                std::cout << "liner 910" << chr <<  std::to_string(position) << "\t" <<  ori << "\t" << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            } else if (sdiRecordsThisOne[j].getChanginglength() == 0 || (sdiRecordsThisOne[j].getReference().compare("-") != 0 && sdiRecordsThisOne[j].getAlternative().compare("-") != 0)) {

            } else {
                std::cout << "maybe something is missing here" << std::endl;
                std::cout << sdiRecordsThisOne[j - 1].getChromosome() << "\t" << std::to_string(sdiRecordsThisOne[j - 1].getPosition()) << "\t" << sdiRecordsThisOne[j - 1].getReference() << "\t" << sdiRecordsThisOne[j - 1].getAlternative() << std::endl;
                std::cout << sdiRecordsThisOne[j].getChromosome() << "\t" << std::to_string(sdiRecordsThisOne[j].getPosition()) << "\t" << sdiRecordsThisOne[j].getReference() << "\t" << sdiRecordsThisOne[j].getAlternative() << std::endl;
            }
        }
        std::cout << sdiRecordsThisOne.size() << std::endl;
//        return;
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
            sdiRecordsThisOne.shrink_to_fit();
        }
        std::cout << sdiRecordsThisOne.size() << std::endl << std::endl;
//        sdiRecordsThisOne.shrink_to_fit();
    }
    std::cout << sdiRecordsThisOne.size() << std::endl;
    std::cout << "merge nearby records together" << std::endl;

    ifChanged = true;
    while (ifChanged && sdiRecordsThisOne.size() > 1) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j - 1].getReference().compare("-") != 0 &&
                sdiRecordsThisOne[j].getChanginglength() > 0 && sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                sdiRecordsThisOne[j - 1].getPosition() == sdiRecordsThisOne[j].getPosition()) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j].getReference() + sdiRecordsThisOne[j - 1].getReference();
                std::string result = sdiRecordsThisOne[j].getAlternative() + sdiRecordsThisOne[j - 1].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//                std::cout << "liner 931" << chr <<  std::to_string(position) <<  ori << result << std::endl;
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }
    std::cout << sdiRecordsThisOne.size() << std::endl;
    std::cout << "merge same location records together" << std::endl;

    if (sdiRecordsThisOne.size() > 0 && sdiRecordsThisOne[0].getChanginglength() > 0 &&
        sdiRecordsThisOne[0].getReference().compare("-") == 0) {
        if (sdiRecordsThisOne[0].getPosition() != 1) {
            int position = sdiRecordsThisOne[0].getPosition() - 1;
            std::string ori(1, refSequence[position - 1]);
            std::string result = ori + sdiRecordsThisOne[0].getAlternative();
//            std::cout << "liner 953" << chr <<  std::to_string(position) <<  ori << result << std::endl;
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        } else if (sdiRecordsThisOne[0].getPosition() == 1 && sdiRecordsThisOne.size() > 1 && sdiRecordsThisOne[1].getPosition() == 1) {
            int position = 1;
            std::string ori = sdiRecordsThisOne[0].getReference() + sdiRecordsThisOne[1].getReference();
            std::string result = sdiRecordsThisOne[0].getAlternative() + sdiRecordsThisOne[1].getAlternative();

            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//            std::cout << "liner 963" << chr <<  std::to_string(position) <<  ori << result << std::endl;
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + 1);
            sdiRecordsThisOne.shrink_to_fit();
            sdiRecordsThisOne[0] = mapSingleRecord2;
        } else {
            int position = 1;
            std::string ori(1, refSequence[0]);
            std::string result = sdiRecordsThisOne[0].getAlternative() + ori;
            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//            std::cout << "liner 974" << chr <<  std::to_string(position) <<  ori << result << std::endl;
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }
    }
    std::cout << "first insertion" << std::endl;
    if (sdiRecordsThisOne.size() > 0 && sdiRecordsThisOne[0].getChanginglength() < 0 &&
        sdiRecordsThisOne[0].getAlternative().compare("-") == 0) {
        if (sdiRecordsThisOne[0].getPosition() != 1) {
            int position = sdiRecordsThisOne[0].getPosition() - 1;
            std::string ori(1, refSequence[position - 1]);
            std::string result = ori;
            ori += sdiRecordsThisOne[0].getReference();
//            std::cout << "liner 987" << chr <<  std::to_string(position) <<  ori << result << std::endl;
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        } else if (sdiRecordsThisOne[0].getPosition() == 1 && sdiRecordsThisOne.size() > 1 && sdiRecordsThisOne[1].getPosition() == sdiRecordsThisOne[0].getPosition() + sdiRecordsThisOne[0].getReference().size()) {
            int position = 1;
            std::string ori = sdiRecordsThisOne[0].getReference() + sdiRecordsThisOne[1].getReference();
            std::string result = sdiRecordsThisOne[0].getAlternative() + sdiRecordsThisOne[1].getAlternative();

            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//            std::cout << "liner 997" << chr <<  std::to_string(position) <<  ori << result << std::endl;
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + 1);
            sdiRecordsThisOne.shrink_to_fit();
            sdiRecordsThisOne[0] = mapSingleRecord2;
        } else {
            int position = 1;
            std::string result(1, refSequence[sdiRecordsThisOne[0].getReference().length()]);
            std::string ori = sdiRecordsThisOne[0].getReference() + result;
            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
//            std::cout << "liner 1008" << chr <<  std::to_string(position) <<  ori << result << std::endl;
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }
    }
    std::cout << "first deletion" << std::endl;
}


void genomeAlignmentAndVariantCallingSingleThread(
        std::string refChr, std::string queryChr,
        std::string refSequence, std::string targetSequence,
        const std::vector<AlignmentMatch> it0,
        const int chrWidth, const std::string refFileName, const std::string queryFileName,
        const bool outPutMaf, const bool outPutVcf, const bool outPutFraged, std::ofstream &omaffile,
        std::ofstream &ovcffile, std::ofstream &ofragfile,
        const int32_t widownWidth, const int32_t wfaSize, const int32_t wfaSize2,
        const int32_t matchingScore, const int32_t mismatchingPenalty,
        const int32_t openGapPenalty1, const int32_t extendGapPenalty1, const int32_t openGapPenalty2, const int32_t extendGapPenalty2,
        const int32_t min_wavefront_length, const int32_t max_distance_threshold, std::atomic_int &number_of_runing_threads,
        std::map<std::string, std::string> parameters) {
//    std::cout << "line 877" << std::endl;

    Scorei m(matchingScore, mismatchingPenalty);

//    std::cout << "line 893" << std::endl;
    size_t startRef = 1;
    size_t startQuery = 1;
    size_t endRef;
    size_t endQuery;
    std::stringstream refAlign;
    std::stringstream queryAlign;
//    std::cout << "line 900" << std::endl;
//    std::cout << refChr << " align start" << std::endl;

    int64_t alignmentScore = 0;
    STRAND lastStrand = POSITIVE;
    AlignmentMatch lastAlignmentMatch;

    int mafRefStart = 0;
    int mafQueryStart = 0;
    std::string mafStrand = "+";
    bool hasInversion = false;
//    std::cout << "line 913" << std::endl;
    bool checkResult = true;
    for (AlignmentMatch alignmentMatch: it0) {
        if (alignmentMatch.getStrand() == NEGATIVE) {
            hasInversion = true;
        }

        if (lastStrand == POSITIVE && alignmentMatch.getStrand() == POSITIVE) {
//                    std::cout << "line 1625" << std::endl;
            if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() != startQuery) {
                endQuery = alignmentMatch.getQueryStartPos() - 1;
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);

                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.length(), '-');
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * querySeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * querySeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q
                              << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryStartPos() == startQuery) {
                endRef = alignmentMatch.getRefStartPos() - 1;
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.length(), '-');
                std::string _alignment_d = refSeq;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * refSeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * refSeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }

            } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() == startQuery) {

            } else {
                endRef = alignmentMatch.getRefStartPos() - 1;
                endQuery = alignmentMatch.getQueryStartPos() - 1;

                std::string refSeq = getSubsequence(refSequence, startRef, endRef);
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);

                std::string _alignment_q;
                std::string _alignment_d;

                int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m,
                                                      parameters);
                if (checkResult) {
                    std::string tempd;
                    std::string tempq;
                    tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                    tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                    if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                        std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                        thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                    }
                    assert(tempd.compare(refSeq) == 0);
                    assert(tempq.compare(querySeq) == 0);
                }


                alignmentScore += thiScore;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;

                //assert(temp.compare(querySeq) == 0);

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q
                              << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
            mafStrand = "+";
        } else if (lastStrand == NEGATIVE && alignmentMatch.getStrand() == NEGATIVE
                   && lastAlignmentMatch.getRefEndPos() < alignmentMatch.getRefStartPos()
                   && lastAlignmentMatch.getQueryStartPos() > alignmentMatch.getQueryEndPos()
                ) {
            if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() != startQuery) {
                endQuery = alignmentMatch.getQueryEndPos() + 1;
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery, alignmentMatch.getStrand());

                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.length(), '-');
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * querySeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * querySeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequence.size() << "\t" << _alignment_q
                              << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }

            } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryEndPos() == startQuery) {
                endRef = alignmentMatch.getRefStartPos() - 1;
                std::string refSeq = getSubsequence(refSequence, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.length(), '-');
                std::string _alignment_d = refSeq;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.length();
                if (thiScore < openGapPenalty2 + extendGapPenalty2 * refSeq.length()) {
                    thiScore = openGapPenalty2 + extendGapPenalty2 * refSeq.length();
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t-\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() == startQuery) {

            } else {
                endRef = alignmentMatch.getRefStartPos() - 1;
                endQuery = alignmentMatch.getQueryEndPos() + 1;

                std::string refSeq = getSubsequence(refSequence, startRef, endRef);
                std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery, alignmentMatch.getStrand());

                std::string _alignment_q;
                std::string _alignment_d;

                int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m,
                                                      parameters);
                if (checkResult) {
                    std::string tempd;
                    std::string tempq;
                    tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                    tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                    if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                        std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                        thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore,
                                                        mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                        openGapPenalty2, extendGapPenalty2, min_wavefront_length,
                                                        max_distance_threshold, m);
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                    }
                    assert(tempd.compare(refSeq) == 0);
                    assert(tempq.compare(querySeq) == 0);
                }

                alignmentScore += thiScore;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << endQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequence.size() << "\t" << _alignment_q
                              << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
//            if ( mafQueryStart > alignmentMatch.getQueryStartPos() ){
//                mafQueryStart = alignmentMatch.getQueryStartPos();
//            }
        } else {

            std::string temp1 = refAlign.str();
            std::string temp2 = queryAlign.str();
            if (outPutMaf && temp1.size() > 0) {
                temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
                temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());
                std::string refGenomerSequence;
                std::string queryGenomerSequence;
                if (lastStrand == POSITIVE) {
                    refGenomerSequence = getSubsequence(refSequence, mafRefStart + 1, mafRefStart + temp1.size());
                    queryGenomerSequence = getSubsequence(targetSequence, mafQueryStart + 1, mafQueryStart + temp2.size(), lastStrand);
                } else {
                    refGenomerSequence = getSubsequence(refSequence, mafRefStart + 1, mafRefStart + temp1.size());
                    queryGenomerSequence = getSubsequence(targetSequence, mafQueryStart - temp2.size() + 2, mafQueryStart + 1, lastStrand);
                }

//                std::cout << temp1 << std::endl;
//                std::cout << refGenomerSequence << std::endl;
//                std::cout << temp2 << std::endl;
//                std::cout << queryGenomerSequence << std::endl;

                assert(temp1.compare(refGenomerSequence) == 0);
                assert(temp2.compare(queryGenomerSequence) == 0);

                int32_t tm = mafQueryStart - temp2.size() + 1;
                if (lastStrand == POSITIVE) {
                    tm = mafQueryStart;
                }
                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right
                         << std::setw(9) << mafRefStart << "\t" << std::setw(9)
                         << temp1.size() << "\t+\t" << refSequence.size() << "\t"
                         << refAlign.str() << std::endl

                         << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right
                         << std::setw(9) << tm << "\t" << std::setw(9)
                         << temp2.size() << "\t" << mafStrand << "\t" << targetSequence.size() << "\t"
                         << queryAlign.str() << std::endl
                         << std::endl;
                g_num_mutex.unlock();
            }
            alignmentScore = 0;
            refAlign.str(std::string());
            queryAlign.str(std::string());

            mafRefStart = alignmentMatch.getRefStartPos() - 1;
            mafQueryStart = alignmentMatch.getQueryStartPos() - 1;
            if (NEGATIVE == alignmentMatch.getStrand()) {
                mafQueryStart = alignmentMatch.getQueryEndPos() - 1;
            }
        }
        {
            if (POSITIVE == alignmentMatch.getStrand()) {
                mafStrand = "+";
            } else {
                mafStrand = "-";
            }
            startRef = alignmentMatch.getRefStartPos();
            startQuery = alignmentMatch.getQueryStartPos();
            endRef = alignmentMatch.getRefEndPos();
            endQuery = alignmentMatch.getQueryEndPos();

            std::string refSeq = getSubsequence(refSequence, startRef, endRef);
            std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery, alignmentMatch.getStrand());

            std::string _alignment_q;
            std::string _alignment_d;

            int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m,
                                                  parameters);
            if (checkResult) {
                std::string tempd;
                std::string tempq;
                tempd = _alignment_d;
                tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                tempq = _alignment_q;
                tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                    std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                    thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore,
                                                    mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                    openGapPenalty2, extendGapPenalty2, min_wavefront_length,
                                                    max_distance_threshold, m);
                    tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                    tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                }
                assert(tempd.compare(refSeq) == 0);
                assert(tempq.compare(querySeq) == 0);
            }
            alignmentScore += thiScore;
            refAlign << _alignment_d;
            queryAlign << _alignment_q;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t" << mafStrand << "\t" << targetSequence.size() << "\t"
                          << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
        }
        startRef = alignmentMatch.getRefEndPos() + 1;
        startQuery = alignmentMatch.getQueryEndPos() + 1;
        if (alignmentMatch.getStrand() == NEGATIVE) {
            startQuery = alignmentMatch.getQueryStartPos() - 1;
        }
        lastStrand = alignmentMatch.getStrand();
        lastAlignmentMatch = alignmentMatch;
    }
//    std::cout << "line 1126" << std::endl;
    if (!hasInversion) {

//        std::cout << refChr << " last align" << std::endl;
        endRef = refSequence.length();
        endQuery = targetSequence.length();
        if (startRef > endRef && startQuery <= endQuery) {
            std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);
            std::string _alignment_q = querySeq;
            std::string _alignment_d = std::string(querySeq.length(), '-');
            refAlign << _alignment_d;
            queryAlign << _alignment_q;
            int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.length();
            if (thiScore < openGapPenalty2 + extendGapPenalty2 * querySeq.length()) {
                thiScore = openGapPenalty2 + extendGapPenalty2 * querySeq.length();
            }
            alignmentScore += thiScore;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }

        } else if (startRef <= endRef && startQuery > endQuery) {
            std::string refSeq = getSubsequence(refSequence, startRef, endRef);

            std::string _alignment_q = std::string(refSeq.length(), '-');
            std::string _alignment_d = refSeq;
            refAlign << _alignment_d;
            queryAlign << _alignment_q;
            int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.length();
            if (thiScore < openGapPenalty2 + extendGapPenalty2 * refSeq.length()) {
                thiScore = openGapPenalty2 + extendGapPenalty2 * refSeq.length();
            }
            alignmentScore += thiScore;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
        } else if (startRef > endRef && startQuery > endQuery) {

        } else {
            std::string refSeq = getSubsequence(refSequence, startRef, endRef);
            std::string querySeq = getSubsequence(targetSequence, startQuery, endQuery);

            std::string _alignment_q;
            std::string _alignment_d;

            int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m,
                                                  parameters);
            if (checkResult) {
                std::string tempd;
                std::string tempq;
                tempd = _alignment_d;
                tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                tempq = _alignment_q;
                tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                    std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                    thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, wfaSize, matchingScore,
                                                    mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                    openGapPenalty2, extendGapPenalty2, min_wavefront_length,
                                                    max_distance_threshold, m);
                    tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                    tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                }
                assert(tempd.compare(refSeq) == 0);
                assert(tempq.compare(querySeq) == 0);
            }
            alignmentScore += thiScore;
            refAlign << _alignment_d;
            queryAlign << _alignment_q;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequence.size() << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequence.size() << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
        }

        std::string temp = refAlign.str();
        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
        assert(temp.compare(refSequence) == 0);

        temp = queryAlign.str();
        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
//                temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
        assert(temp.compare(targetSequence) == 0);
    }

    if (outPutMaf) {

        std::string temp1 = refAlign.str();
        std::string temp2 = queryAlign.str();
        temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
        temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());

        std::string refGenomerSequence;
        std::string queryGenomerSequence;
        if (lastStrand == POSITIVE) {
            refGenomerSequence = getSubsequence(refSequence, mafRefStart + 1, mafRefStart + temp1.size());
            queryGenomerSequence = getSubsequence(targetSequence, mafQueryStart + 1, mafQueryStart + temp2.size(), lastStrand);
        } else {
            refGenomerSequence = getSubsequence(refSequence, mafRefStart + 1, mafRefStart + temp1.size());
            queryGenomerSequence = getSubsequence(targetSequence, mafQueryStart + 1, mafQueryStart - temp2.size() + 2, lastStrand);
        }


//        std::cout << temp1 << std::endl;
//        std::cout << refGenomerSequence << std::endl;
//        std::cout << temp2 << std::endl;
//        std::cout << queryGenomerSequence << std::endl;
        assert(temp1.compare(refGenomerSequence) == 0);
        assert(temp2.compare(queryGenomerSequence) == 0);

        int32_t tm = mafQueryStart - temp2.size() + 1;
        if (lastStrand == POSITIVE) {
            tm = mafQueryStart;
        }

        g_num_mutex.lock();
        omaffile << "a\tscore=" << alignmentScore << std::endl
                 << "s\t" << std::left << std::setw(chrWidth) << /*refFileName + "." + */refChr << "\t" << std::right
                 << std::setw(9) << mafRefStart << "\t" << std::setw(9)
                 << temp1.size() << "\t+\t" << refSequence.size() << "\t"
                 << refAlign.str() << std::endl
                 << "s\t" << std::left << std::setw(chrWidth) << /*queryFileName + "." + */queryChr << "\t" << std::right
                 << std::setw(9) << tm << "\t" << std::setw(9)
                 << temp2.size() << "\t" + mafStrand + "\t" << targetSequence.size() << "\t"
                 << queryAlign.str() << std::endl
                 << std::endl;
        g_num_mutex.unlock();
    }

    if (outPutVcf) {
        std::string queryAlignSeq = queryAlign.str();
        std::string refAlignSeq = refAlign.str();
        alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, refChr, refSequence);
    }
//    std::cout << refChr << " align done" << std::endl;
//    std::cout << "line 1223" << std::endl;
    --number_of_runing_threads;
}


void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                      const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                                      const int32_t &widownWidth, const int32_t &wfaSize, const int32_t &wfaSize2, const std::string &outPutMafFile, const std::string &outPutVcfFile,
                                      const std::string &outPutFragedFile, /*std::string & outPutLocalalignmentFile, */const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                                      const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                      const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, int32_t &seed_window_size, const int32_t &mini_cns_score, const int32_t &step_size,
                                      const int32_t &matrix_boundary_distance, const int32_t &scoreThreshold, const int32_t &w, const int32_t &xDrop, const int &maxThread, std::map<std::string, std::string> &parameters) {

//    Scorei m(matchingScore, mismatchingPenalty);


    bool outPutMaf = false;
    bool outPutVcf = false;
    bool outPutFraged = false;
//    bool outPutLocalalignment = false;

    if (outPutMafFile.size() > 0) {
        outPutMaf = true;
    }
    if (outPutVcfFile.size() > 0) {
        outPutVcf = true;
    }
    if (outPutFragedFile.size() > 0) {
        outPutFraged = true;
    }
//    if ( outPutLocalalignmentFile.size() > 0 ){
//        outPutLocalalignment = true;
//    }

    std::map<std::string, std::string> refSequences;
    readFastaFile(refFastaFilePath, refSequences);

    std::map<std::string, std::string> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);

    int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    queryFileName = elems.back();

    for (std::map<std::string, std::string>::iterator itchr = refSequences.begin(); itchr != refSequences.end(); ++itchr) {
        if ((refFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }
    for (std::map<std::string, std::string>::iterator itchr = targetSequences.begin(); itchr != targetSequences.end(); ++itchr) {
        if ((queryFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ovcffile;
    std::ofstream ofragfile;
    if (outPutMaf) {
        omaffile.open(outPutMafFile);
        omaffile << "##maf version=1" << std::endl;
    }

    SequenceCharToUInt8 sequenceCharToUInt8;

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon));
    if (ltm->tm_mday < 10) {
        filedate = filedate + "0" + std::to_string((ltm->tm_mday));
    } else {
        filedate = filedate + std::to_string((ltm->tm_mday));
    }

    if (outPutVcf) {
        ovcffile.open(outPutVcfFile);
        ovcffile << "##fileformat=VCFv4.3" << std::endl;
        ovcffile << "##fileDate=" << filedate << std::endl;
        ovcffile << "##source=proali" << std::endl;
        ovcffile << "##reference=" << refFileName << std::endl;
        ovcffile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
        ovcffile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
        ovcffile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
        ovcffile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
        ovcffile << "##FILTER=<ID=q30,Description=\"Quality below 30\">" << std::endl;
        std::string accession = queryFileName;
        accession = songStrReplaceAll(accession, ".fasta", "");
        accession = songStrReplaceAll(accession, ".fa", "");
        accession.erase(std::remove(accession.begin(), accession.end(), ' '), accession.end());
        ovcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << accession << std::endl;
    }

    if (outPutFraged) {
        ofragfile.open(outPutFragedFile);
        ofragfile << "##maf version=1" << std::endl;
    }
    std::atomic_int number_of_runing_threads(0);
    for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it0 = alignmentMatchsMap.begin(); it0 != alignmentMatchsMap.end(); ++it0) {
        if (it0->second.size() > 0) {
            bool isThisThreadUnrun = true;
            std::string refChr = it0->second[0].getRefChr();
            std::string queryChr = it0->second[0].getQueryChr();

            while (isThisThreadUnrun) {
                if (number_of_runing_threads < maxThread) {
                    std::thread t(genomeAlignmentAndVariantCallingSingleThread, refChr, queryChr, refSequences[refChr], targetSequences[queryChr], it0->second, chrWidth,
                                  refFileName, queryFileName, outPutMaf, outPutVcf, outPutFraged, std::ref(omaffile),
                                  std::ref(ovcffile), std::ref(ofragfile), widownWidth, wfaSize, wfaSize2, matchingScore, mismatchingPenalty,
                                  openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                                  min_wavefront_length, max_distance_threshold, std::ref(number_of_runing_threads), parameters);
                    ++number_of_runing_threads;
                    t.detach();
                    isThisThreadUnrun = false;
                    break;
                } else {
                    usleep(1000);
                }
            }
        }
    }
    while (number_of_runing_threads > 0) {// wait for all the thread
        usleep(1000);
    }

    if (outPutMaf) {
        omaffile.close();
    }
    if (outPutVcf) {
        ovcffile.close();
    }
    if (outPutFraged) {
        ofragfile.close();
    }
}
