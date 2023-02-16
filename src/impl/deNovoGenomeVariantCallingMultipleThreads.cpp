//
// Created by Baoxing song on 20.10.18.
//

#include "deNovoGenomeVariantCalling.h"

std::mutex g_num_mutex;

void genomeAlignmentSingleThread(std::vector<AlignmentMatch> alignmentMatchs,
                                 const bool outPutMaf, const bool outPutFraged,
                                 std::ofstream &omaffile, std::ofstream &ofragfile,
                                 std::map<std::string, std::tuple<std::string, long, long, int> > &map_ref,
                                 std::map<std::string, std::tuple<std::string, long, long, int> > &map_qry,
                                 const int chrWidth, const std::string refFileName, const std::string queryFileName,
                                 const int32_t windowWidth,
                                 const int32_t matchingScore, const int32_t mismatchingPenalty,
                                 const int32_t openGapPenalty1, const int32_t extendGapPenalty1,
                                 const int32_t openGapPenalty2, const int32_t extendGapPenalty2,
                                 std::atomic_int &number_of_runing_threads) {

    std::string refChr = alignmentMatchs[0].getRefChr();
    std::string queryChr = alignmentMatchs[0].getQueryChr();

    bool checkResult = false;

    STRAND strand = alignmentMatchs[0].getStrand();

    size_t size_refSequence = getSequenceSizeFromPath2(map_ref[refChr]);
    size_t size_targetSequence = getSequenceSizeFromPath2(map_qry[queryChr]);

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
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.size(), '-');
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * querySeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryStartPos() == startQuery) {
                endRef = orthologPair.getRefStartPos() - 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.size(), '-');
                std::string _alignment_d = refSeq;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * refSeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() == startQuery) {

            } else {
                endRef = orthologPair.getRefStartPos() - 1;
                endQuery = orthologPair.getQueryStartPos() - 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;

                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;

                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());

                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
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
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
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
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                {
                    std::string _alignment_q;
                    std::string _alignment_d;

                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    if (checkResult) {
                        std::string tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        std::string tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
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
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            startRef = orthologPair.getRefEndPos() + 1;
            startQuery = orthologPair.getQueryEndPos() + 1;
        }

        {
            std::string refGenomeSequence = getSubsequence2(map_ref, refChr, alignmentMatchs[0].getRefStartPos(), alignmentMatchs[alignmentMatchs.size() - 1].getRefEndPos());
            std::string queryGenomeSequence = getSubsequence2(map_qry, queryChr, alignmentMatchs[0].getQueryStartPos(), alignmentMatchs[alignmentMatchs.size() - 1].getQueryEndPos());

            //if (checkResult) {
            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            assert(temp.compare(refGenomeSequence) == 0);

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            assert(temp.compare(queryGenomeSequence) == 0);
            //}
            if (outPutMaf) {
                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getRefStartPos() - 1 << "\t" << std::setw(9) << refGenomeSequence.size() << "\t+\t" << size_refSequence << "\t"
                         << refAlign.str() << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getQueryStartPos() - 1 << "\t" << std::setw(9) << queryGenomeSequence.size() << "\t+\t"
                         << size_targetSequence << "\t" << queryAlign.str() << std::endl
                         << std::endl;
                g_num_mutex.unlock();
            }
        }
    } else {
        size_t startRef = alignmentMatchs[0].getRefStartPos();
        size_t startQuery=0;
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
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);
                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.size(), '-');
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * querySeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryEndPos() == endQuery) {
                endRef = orthologPair.getRefStartPos() - 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.size(), '-');
                std::string _alignment_d = refSeq;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * refSeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t-\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() == endQuery) {

            } else {
                endRef = orthologPair.getRefStartPos() - 1;
                startQuery = orthologPair.getQueryEndPos() + 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;

                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
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
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << size_targetSequence << "\t" << _alignment_q
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
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;
//                    std::cout << refSeq << std::endl << querySeq << std::endl << "line 276" << std::endl;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

                    if (checkResult) {
                        std::string tempd;
                        std::string tempq;
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                        if (tempd.compare(refSeq) != 0 || tempq.compare(querySeq) != 0) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
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
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            startRef = orthologPair.getRefEndPos() + 1;
            endQuery = orthologPair.getQueryStartPos() - 1;
        }
        {
            std::string refGenomerSequence = getSubsequence2(map_ref, refChr, alignmentMatchs[0].getRefStartPos(), alignmentMatchs[alignmentMatchs.size() - 1].getRefEndPos());
            std::string queryGenomerSequence = getSubsequence2(map_qry, queryChr, alignmentMatchs[0].getQueryEndPos(), alignmentMatchs[alignmentMatchs.size() - 1].getQueryStartPos(), strand);
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
                         << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getRefStartPos() - 1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << size_refSequence << "\t"
                         << refAlign.str() << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << alignmentMatchs[alignmentMatchs.size() - 1].getQueryStartPos() - 1 << "\t" << std::setw(9) << queryGenomerSequence.size()
                         << "\t-\t" << size_targetSequence << "\t" << queryAlign.str() << std::endl
                         << std::endl;
                g_num_mutex.unlock();
            }
        }
    }

    number_of_runing_threads = number_of_runing_threads - 1;
}

void genomeAlignment(std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                     const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                     const int32_t &windowWidth,
                     const std::string &outPutMafFile, const std::string &outPutFragedFile,
                     const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                     const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1,
                     const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                     const int &maxThread) {

    bool outPutMaf = false;
    bool outPutFraged = false;

    if (outPutMafFile.size() > 0) {
        outPutMaf = true;
    }
    if (outPutFragedFile.size() > 0) {
        outPutFraged = true;
    }

    std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
    readFastaFile(refFastaFilePath, map_ref);

    std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
    readFastaFile(targetFastaFilePath, map_qry);

    long unsigned int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    queryFileName = elems.back();

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_ref.begin(); itchr != map_ref.end(); ++itchr) {
        if ((refFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }
    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_qry.begin(); itchr != map_qry.end(); ++itchr) {
        if ((queryFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ofragfile;

    if (outPutMaf) {
        omaffile.open(outPutMafFile);
        omaffile << "##maf version=1" << std::endl;
    }

    if (outPutFraged) {
        ofragfile.open(outPutFragedFile);
        ofragfile << "##maf version=1" << std::endl;
    }

    int32_t size = alignmentMatchsMap.size();

    std::atomic_int number_of_runing_threads(0);

    for (int32_t i = size - 1; i >= 0; --i) {
        std::vector<AlignmentMatch> alignmentMatchs = alignmentMatchsMap[i];
        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if (number_of_runing_threads < maxThread) {
                std::thread t(genomeAlignmentSingleThread, alignmentMatchs, outPutMaf, outPutFraged,
                              std::ref(omaffile), std::ref(ofragfile),
                              std::ref(map_ref), std::ref(map_qry),
                              chrWidth, refFileName, queryFileName,
                              windowWidth, matchingScore, mismatchingPenalty,
                              openGapPenalty1, extendGapPenalty1,
                              openGapPenalty2, extendGapPenalty2,
                              std::ref(number_of_runing_threads));
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

void genomeAlignmentAndVariantCallingSingleThread(
        std::map<std::string, std::tuple<std::string, long, long, int> > & map_ref,
        std::map<std::string, std::tuple<std::string, long, long, int> > & map_qry,
        const std::vector<AlignmentMatch> & v_am,
        const int & chrWidth, const std::string & refFileName, const std::string & queryFileName,
        const bool & outPutMaf, const bool & outPutFraged, std::ofstream &omaffile,
        std::ofstream &ofragfile,
        const int32_t & windowWidth,
        const int32_t & matchingScore, const int32_t & mismatchingPenalty,
        const int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
        std::atomic_int &num_runing_threads) {

    std::string refChr = v_am[0].getRefChr();
    std::string queryChr = v_am[0].getQueryChr();

    size_t startRef = 1;
    size_t startQuery = 1;
    size_t endRef;
    size_t endQuery;
    std::stringstream refAlign;
    std::stringstream queryAlign;

    int64_t alignmentScore = 0;
    STRAND lastStrand = POSITIVE;
    AlignmentMatch lastAlignmentMatch;

    int mafRefStart = 0;
    int mafQueryStart = 0;
    std::string mafStrand = "+";
    bool hasInversion = false;
    bool checkResult = false;

    size_t size_ref_sq = getSequenceSizeFromPath2(map_ref[refChr]);
    size_t size_target_sq = getSequenceSizeFromPath2(map_qry[queryChr]);

    std::string path_ref = std::get<0>(map_ref[refChr]);
    std::string path_qry = std::get<0>(map_qry[queryChr]);

    int fd_ref = open(path_ref.c_str(), O_RDONLY);
    int fd_qry = open(path_qry.c_str(), O_RDONLY);

    for (AlignmentMatch alignmentMatch: v_am) {
        if (alignmentMatch.getStrand() == NEGATIVE) {
            hasInversion = true;
        }

        if (lastStrand == POSITIVE && alignmentMatch.getStrand() == POSITIVE) {
            if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() != startQuery) {
                endQuery = alignmentMatch.getQueryStartPos() - 1;
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

                std::string _alignment_q = seq_qry;
                std::string _alignment_d = std::string(seq_qry.size(), '-');
                if ( outPutMaf || checkResult ) {
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;
                }
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_qry.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_qry.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryStartPos() == startQuery) {
                endRef = alignmentMatch.getRefStartPos() - 1;
                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(seq_ref.size(), '-');
                std::string _alignment_d = seq_ref;
                if (outPutMaf || checkResult ) {
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;
                }
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_ref.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_ref.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() == startQuery) {

            }
            else {
                endRef = alignmentMatch.getRefStartPos() - 1;
                endQuery = alignmentMatch.getQueryStartPos() - 1;

                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

                std::string _alignment_q;
                std::string _alignment_d;

                int64_t thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                if (checkResult) {
                    std::string tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());

                    std::string tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                    if (tempd.compare(seq_ref) != 0 || tempq.compare(seq_qry) != 0) {
                        thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                    }

                    assert(tempd == seq_ref);
                    assert(tempq == seq_qry);
                }

                alignmentScore += thiScore;
                if (outPutMaf || checkResult  ) {
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;
                }

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
            mafStrand = "+";
        }
        else if (lastStrand == NEGATIVE && alignmentMatch.getStrand() == NEGATIVE && lastAlignmentMatch.getRefEndPos() < alignmentMatch.getRefStartPos() && lastAlignmentMatch.getQueryStartPos() > alignmentMatch.getQueryEndPos() ) {
            if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() != startQuery) {
                endQuery = alignmentMatch.getQueryEndPos() + 1;
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                std::string _alignment_q = seq_qry;
                std::string _alignment_d = std::string(seq_qry.size(), '-');
                if (outPutMaf  || checkResult ) {
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;
                }
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_qry.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_qry.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t-\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryEndPos() == startQuery) {
                endRef = alignmentMatch.getRefStartPos() - 1;
                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(seq_ref.size(), '-');
                std::string _alignment_d = seq_ref;
                if (outPutMaf  || checkResult ) {
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;
                }
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_ref.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_ref.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }

                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t-\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() == startQuery) {

            }
            else {
                endRef = alignmentMatch.getRefStartPos() - 1;
                endQuery = alignmentMatch.getQueryEndPos() + 1;

                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                std::string _alignment_q;
                std::string _alignment_d;

                int64_t thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                if (checkResult) {
                    std::string tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());

                    std::string tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                    if (tempd.compare(seq_ref) != 0 || tempq.compare(seq_qry) != 0) {
                        thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                        tempd = _alignment_d;
                        tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                        tempq = _alignment_q;
                        tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                    }
                    assert(tempd == seq_ref);
                    assert(tempq == seq_qry);
                }

                alignmentScore += thiScore;
                if (outPutMaf  || checkResult  ) {
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;
                }

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << endQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t-\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
            }
        }
        else {
            std::string temp1 = refAlign.str();
            std::string temp2 = queryAlign.str();
            if (outPutMaf && temp1.size() > 0) {
                temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
                temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());

                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, mafRefStart + 1, mafRefStart + temp1.size());
                std::string seq_qry;

                if (lastStrand == POSITIVE) {
                    seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, mafQueryStart + 1, mafQueryStart + temp2.size(), lastStrand);
                }
                else {
                    seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, mafQueryStart - temp2.size() + 2, mafQueryStart + 1, lastStrand);
                }

                assert(temp1 == seq_ref);
                assert(temp2 == seq_qry);

                int32_t tm = mafQueryStart - temp2.size() + 1;
                if (lastStrand == POSITIVE) {
                    tm = mafQueryStart;
                }

                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << mafRefStart << "\t" << std::setw(9) << temp1.size() << "\t+\t" << size_ref_sq << "\t" << refAlign.str() << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << tm << "\t" << std::setw(9) << temp2.size() << "\t" << mafStrand << "\t" << size_target_sq << "\t" << queryAlign.str() << std::endl
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

            std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
            std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

            std::string _alignment_q;
            std::string _alignment_d;

            int64_t thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
            if (checkResult) {
                std::string tempd = _alignment_d;
                tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());

                std::string tempq = _alignment_q;
                tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                if (tempd.compare(seq_ref) != 0 || tempq.compare(seq_qry) != 0) {
                    thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                    tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                }

                assert(tempd == seq_ref);
                assert(tempq == seq_qry);
            }

            alignmentScore += thiScore;
            if (outPutMaf  || checkResult ) {
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
            }

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t" << mafStrand << "\t" << size_target_sq << "\t" << _alignment_q << std::endl
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


    // last align
    if (!hasInversion) {
        endRef = getSequenceSizeFromPath2(map_ref[refChr]);
        endQuery = getSequenceSizeFromPath2(map_qry[queryChr]);

        if (startRef > endRef && startQuery <= endQuery) {
            std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

            std::string _alignment_q = seq_qry;
            std::string _alignment_d = std::string(seq_qry.size(), '-');
            if (outPutMaf  || checkResult  ) {
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
            }
            int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_qry.size();
            int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_qry.size();
            if (thiScore < thiScore2) {
                thiScore = thiScore2;
            }

            alignmentScore += thiScore;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
        }
        else if (startRef <= endRef && startQuery > endQuery) {
            std::string refSeq = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);

            std::string _alignment_q = std::string(refSeq.size(), '-');
            std::string _alignment_d = refSeq;
            if (outPutMaf  || checkResult ) {
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
            }

            int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.size();
            int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * refSeq.size();
            if (thiScore < thiScore2) {
                thiScore = thiScore2;
            }
            alignmentScore += thiScore;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
        }
        else if (startRef > endRef && startQuery > endQuery) {

        }
        else {
            std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
            std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

            std::string _alignment_q;
            std::string _alignment_d;

            int64_t thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
            if (checkResult) {
                std::string tempd = _alignment_d;
                tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());

                std::string tempq = _alignment_q;
                tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());

                if (tempd.compare(seq_ref) != 0 || tempq.compare(seq_qry) != 0) {
                    thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    tempd = _alignment_d;
                    tempd.erase(std::remove(tempd.begin(), tempd.end(), '-'), tempd.end());
                    tempq = _alignment_q;
                    tempq.erase(std::remove(tempq.begin(), tempq.end(), '-'), tempq.end());
                }
                assert(tempd == seq_ref);
                assert(tempq == seq_qry);
            }

            alignmentScore += thiScore;
            if (outPutMaf  || checkResult  ) {
                refAlign << _alignment_d;
                queryAlign << _alignment_q;
            }

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
        }
        if (outPutMaf  || checkResult  ) {
            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());

            std::string r2 = getSubsequence3(map_ref, fd_ref, refChr);
            assert(temp == r2);

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());

            std::string t2 = getSubsequence3(map_qry, fd_qry, queryChr);
            assert(temp == t2);
        }
    }

    if (outPutMaf) {
        std::string temp1 = refAlign.str();
        std::string temp2 = queryAlign.str();
        temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
        temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());

        std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, mafRefStart + 1, mafRefStart + temp1.size());
        std::string seq_qry;

        if (lastStrand == POSITIVE) {
            seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, mafQueryStart + 1, mafQueryStart + temp2.size(), lastStrand);
        }
        else {
            seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, mafQueryStart + 1, mafQueryStart - temp2.size() + 2, lastStrand);
        }

        assert(temp1 == seq_ref);
        assert(temp2 == seq_qry);

        int32_t tm = mafQueryStart - temp2.size() + 1;
        if (lastStrand == POSITIVE) {
            tm = mafQueryStart;
        }

        g_num_mutex.lock();
        omaffile << "a\tscore=" << alignmentScore << std::endl
                 << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << mafRefStart << "\t" << std::setw(9) << temp1.size() << "\t+\t" << size_ref_sq << "\t" << refAlign.str() << std::endl
                 << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << tm << "\t" << std::setw(9) << temp2.size() << "\t" + mafStrand + "\t" << size_target_sq << "\t" << queryAlign.str() << std::endl
                 << std::endl;
        g_num_mutex.unlock();
    }

    close(fd_ref);
    close(fd_qry);

    --num_runing_threads;
}


void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> &map_v_am,
                                      const std::string &path_ref_GenomeSequence, const std::string &path_target_GenomeSequence,
                                      const int32_t &windowWidth, const std::string &outPutMafFile,
                                      const std::string &outPutFragedFile, const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                                      const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                      const int &maxThread) {

    bool outPutMaf = false;
    bool outPutFraged = false;

    if (outPutMafFile.size() > 0) {
        outPutMaf = true;
    }

    if (outPutFragedFile.size() > 0) {
        outPutFraged = true;
    }

    std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
    readFastaFile(path_ref_GenomeSequence, map_ref);

    std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
    readFastaFile(path_target_GenomeSequence, map_qry);

    long unsigned int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(path_ref_GenomeSequence, delim, elems);
    refFileName = elems.back();
    split(path_target_GenomeSequence, delim, elems);
    queryFileName = elems.back();

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_ref.begin(); itchr != map_ref.end(); ++itchr) {
        if ((refFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_qry.begin(); itchr != map_qry.end(); ++itchr) {
        if ((queryFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ofragfile;
    if (outPutMaf) {
        omaffile.open(outPutMafFile);
        omaffile << "##maf version=1" << std::endl;
    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon));
    if (ltm->tm_mday < 10) {
        filedate = filedate + "0" + std::to_string((ltm->tm_mday));
    } else {
        filedate = filedate + std::to_string((ltm->tm_mday));
    }

    if (outPutFraged) {
        ofragfile.open(outPutFragedFile);
        ofragfile << "##maf version=1" << std::endl;
    }

    std::atomic_int num_runing_threads(0);

    for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
        if (it->second.size() > 0) {
            bool isThisThreadUnrun = true;

            while (isThisThreadUnrun) {
                if (num_runing_threads < maxThread) {

//                    std::thread t(genomeAlignmentAndVariantCallingSingleThread, std::ref(refChr), std::ref(queryChr), std::ref(map_ref), std::ref(map_qry), std::ref(it->second), std::ref(chrWidth),
//                                  std::ref(refFileName), std::ref(queryFileName), std::ref(outPutMaf), std::ref(outPutFraged), std::ref(omaffile),
//                                  std::ref(ofragfile), std::ref(windowWidth), std::ref(matchingScore), std::ref(mismatchingPenalty),
//                                  std::ref(openGapPenalty1), std::ref(extendGapPenalty1), std::ref(openGapPenalty2), std::ref(extendGapPenalty2),
//                                  std::ref(num_runing_threads));
                    std::thread t(genomeAlignmentAndVariantCallingSingleThread, std::ref(map_ref), std::ref(map_qry), std::ref(it->second), std::ref(chrWidth),
                                  std::ref(refFileName), std::ref(queryFileName), std::ref(outPutMaf), std::ref(outPutFraged), std::ref(omaffile),
                                  std::ref(ofragfile), std::ref(windowWidth), std::ref(matchingScore), std::ref(mismatchingPenalty),
                                  std::ref(openGapPenalty1), std::ref(extendGapPenalty1), std::ref(openGapPenalty2), std::ref(extendGapPenalty2),
                                  std::ref(num_runing_threads));

                    ++num_runing_threads;
                    t.detach();
                    isThisThreadUnrun = false;
                    break;
                } else {
                    usleep(1000);
                }
            }
        }
    }

    while (num_runing_threads > 0) {// wait for all the thread
        usleep(1000);
    }

    if (outPutMaf) {
        omaffile.close();
    }

    if (outPutFraged) {
        ofragfile.close();
    }
}
