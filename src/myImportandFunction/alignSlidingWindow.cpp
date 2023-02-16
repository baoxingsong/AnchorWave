//
// Created by song on 8/5/18.
//

#include "alignSlidingWindow.h"

int32_t minimap2_alignment(const std::string &_dna_q, const std::string &_dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int32_t &matchingScore, int32_t mismatchingPenalty,
                           int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2) {
    int8_t a = 0, b = mismatchingPenalty < 0 ? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = _dna_d.length(), ql = _dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = _dna_d.c_str();
    const char *qseq = _dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);

    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table

    ts = (uint8_t *) malloc(tl);
    qs = (uint8_t *) malloc(ql);

    int i;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t) qseq[i]];
    }

    //flag  0x01 score only

#ifdef __AVX512BW__
    //    std::cout << "using AVX512" << std::endl;
        ksw_extd2_avx512(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0, & ez);
#elif __AVX2__
    //    std::cout << "using AVX2" << std::endl;
    ksw_extd2_avx2(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0, & ez);
#else
    //    std::cout << "using SSE" << std::endl;
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0, &ez);
#endif
    std::string cigarstring = "";
    for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
        //        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
        cigarstring = cigarstring + std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];
    }
    //    putchar('\n');

    std::vector <std::string> cigarElems;
    splitCIGAR(cigarstring, cigarElems);

    _alignment_q = "";
    _alignment_d = "";
    int pattern_pos = 0, text_pos = 0;

    for (i = 0; i < cigarElems.size(); ++i) {
        std::string cVal = cigarElems[i];
        char cLetter = cVal[cVal.length() - 1];
        int cLen = stoi(cVal.substr(0, cVal.length() - 1));

        if (cLetter == 'M') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += _dna_q[text_pos];
                _alignment_d += _dna_d[pattern_pos];
                pattern_pos++;
                text_pos++;
            }
        } else if (cLetter == 'X') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += _dna_q[text_pos];
                _alignment_d += _dna_d[pattern_pos];
                pattern_pos++;
                text_pos++;
            }
        } else if (cLetter == 'I') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += _dna_q[text_pos];
                _alignment_d += '-';
                text_pos++;
            }
        } else if (cLetter == 'D') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += '-';
                _alignment_d += _dna_d[pattern_pos];
                pattern_pos++;
            }
        }
    }
    int32_t totalScore = ez.score;
    free(ez.cigar);
    free(ts);
    free(qs);
    return totalScore;
}


int32_t alignment_minimap2(const std::string &_dna_q, const std::string &_dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int32_t &matchingScore, int32_t mismatchingPenalty, int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2, int32_t &endPositionq, int32_t &endPositiont) {

    int8_t a = 0, b = mismatchingPenalty < 0 ? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = _dna_d.length(), ql = _dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = _dna_d.c_str();
    const char *qseq = _dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);

    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table

    ts = (uint8_t *) malloc(tl);
    qs = (uint8_t *) malloc(ql);

    int i;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t) qseq[i]];
    }

//flag  0x01 score only

#ifdef __AVX512BW__
    //    std::cout << "using AVX512" << std::endl;
    ksw_extd2_avx512(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0x01, & ez);
#elif __AVX2__
    //    std::cout << "using AVX2" << std::endl;
        ksw_extd2_avx2(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0x01, & ez);
#else
    //    std::cout << "using SSE" << std::endl;
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0x01, &ez);
#endif

    endPositionq = 0;
    endPositiont = 0;
    int32_t maxScore = -1000000;

    if (ez.mqe > ez.mte) {
        endPositiont = ez.mqe_t + 1;
        endPositionq = ql;
        maxScore = ez.mqe;
    } else {
        endPositiont = tl;
        endPositionq = ez.mte_q + 1;
        maxScore = ez.mte;
    }
    free(ez.cigar);
    free(ts);
    free(qs);
    std::string _dna_q_sub = _dna_q.substr(0, endPositionq);
    std::string _dna_d_sub = _dna_d.substr(0, endPositiont);

    int32_t thisScore = minimap2_alignment(_dna_q_sub, _dna_d_sub, _alignment_q, _alignment_d,
                                           matchingScore, mismatchingPenalty,
                                           _open_gap_penalty1, _extend_gap_penalty1,
                                           _open_gap_penalty2, _extend_gap_penalty2);
    assert(thisScore == maxScore);
    return maxScore;
}

int64_t alignSlidingWindow(const std::string &dna_q, const std::string &dna_d, int64_t _length_of_q, int64_t _length_of_d,
                           std::string &_alignment_q, std::string &_alignment_d, const int64_t &slidingWindowSize,
                           const int32_t &matchingScore, const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                           const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    //2^15 = 32768
    //of the maximum length of the windowSize of is about 32000/2 = 16000
    int32_t databaseStart = 1;
    int32_t databaseEnd = 0;
    int32_t queryStart = 1;
    int32_t queryEnd = 0;
    int64_t totalScore = 0;
    int32_t endPositionq;
    int32_t endPositiont;

    //if (sqrt(_length_of_d) * sqrt(_length_of_q) <= slidingWindowSize) {
    if ((_length_of_d * 1.0 / slidingWindowSize) * (_length_of_q * 1.0 / slidingWindowSize) <= 1) {
        totalScore += minimap2_alignment(dna_q, dna_d, _alignment_q, _alignment_d,
                                         matchingScore, mismatchingPenalty,
                                         openGapPenalty1, extendGapPenalty1,
                                         openGapPenalty2, extendGapPenalty2);
        databaseStart = _length_of_d + 1;
        queryStart = _length_of_q + 1;
    } else {
        while (databaseStart <= _length_of_d && queryStart <= _length_of_q) {
            databaseEnd = databaseStart + slidingWindowSize;
            queryEnd = queryStart + slidingWindowSize;
            if (databaseEnd > _length_of_d) {
                databaseEnd = _length_of_d;
            }

            if (queryEnd > _length_of_q) {
                queryEnd = _length_of_q;
            }

            std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd);
            std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd);
            std::string alignment_q = "";
            std::string alignment_d = "";

            if (slidingWindowSize > 1073741824) {
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            } else {
                if (extendGapPenalty2 + matchingScore < 0) {
                    totalScore += alignment_minimap2(qSeq, dSeq, alignment_q, alignment_d, matchingScore, mismatchingPenalty + matchingScore, openGapPenalty1 + matchingScore, extendGapPenalty1 + matchingScore, openGapPenalty2 + matchingScore, extendGapPenalty2 + matchingScore, endPositionq, endPositiont);
                } else {
                    totalScore += alignment_minimap2(qSeq, dSeq, alignment_q, alignment_d, matchingScore, mismatchingPenalty + matchingScore - 1, openGapPenalty1 + matchingScore - 1, extendGapPenalty1 + matchingScore - 1, openGapPenalty2 + matchingScore - 1, -1, endPositionq, endPositiont);
                }
            }

            std::string temp_q = alignment_q;
            temp_q.erase(std::remove(temp_q.begin(), temp_q.end(), '-'), temp_q.end());
            assert(endPositionq == temp_q.size());

            std::string temp_d = alignment_d;
            temp_d.erase(std::remove(temp_d.begin(), temp_d.end(), '-'), temp_d.end());
            assert(endPositiont == temp_d.size());

            _alignment_q += alignment_q;
            _alignment_d += alignment_d;
            queryStart += endPositionq;
            databaseStart += endPositiont;
        }
    }

    int32_t final_indel_length = 0;
    int32_t count_1 = _length_of_d - databaseStart;
    if (count_1 >= 0) {
        _alignment_q += std::string(count_1 + 1, '-');
        _alignment_d += dna_d.substr(databaseStart - 1, count_1 + 1);
        final_indel_length += count_1 + 1;
    }
    assert(_alignment_d.size() == _alignment_q.size());

    int32_t count_2 = _length_of_q - queryStart;
    if (count_2 >= 0) {
        _alignment_q += dna_q.substr(queryStart - 1, count_2 + 1);
        _alignment_d += std::string(count_2 + 1, '-');
        final_indel_length += count_2 + 1;
    }
    assert(_alignment_d.size() == _alignment_q.size());

    if (final_indel_length > 0) {
        totalScore += std::max(openGapPenalty1 + extendGapPenalty1 * final_indel_length, openGapPenalty2 + extendGapPenalty2 * final_indel_length);
    }

    return totalScore;
}


int64_t alignSlidingWindow_minimap2(const std::string &dna_q, const std::string &dna_d, int64_t _length_of_q, int64_t _length_of_d,
                                    std::string &_alignment_q, std::string &_alignment_d, const int64_t &slidingWindowSize,
                                    const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                                    const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    int8_t a = 0, b = mismatchingPenalty < 0 ? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = dna_d.length(), ql = dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = dna_d.c_str();
    const char *qseq = dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);

    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table

    ts = (uint8_t *) malloc(tl);
    qs = (uint8_t *) malloc(ql);

    int i;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t) qseq[i]];
    }


#ifdef __AVX512BW__
    ksw_extd2_avx512(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, & ez);
#elif __AVX2__
    ksw_extd2_avx2(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, &ez);
#else
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, &ez);
#endif

    std::string cigarstring = "";
    for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
        cigarstring = cigarstring + std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];
    }

    std::vector <std::string> cigarElems;
    splitCIGAR(cigarstring, cigarElems);

    _alignment_q = "";
    _alignment_d = "";
    int pattern_pos = 0, text_pos = 0;

    for (i = 0; i < cigarElems.size(); ++i) {
        std::string cVal = cigarElems[i];
        char cLetter = cVal[cVal.length() - 1];
        int cLen = stoi(cVal.substr(0, cVal.length() - 1));

        if (cLetter == 'M' || cLetter == 'X') {
            if (cLen >= 1) {
                _alignment_q += dna_q.substr(text_pos, cLen);
                _alignment_d += dna_d.substr(pattern_pos, cLen);
                pattern_pos += cLen;
                text_pos += cLen;
            }
        } else if (cLetter == 'I') {
            if (cLen >= 1) {
                _alignment_q += dna_q.substr(text_pos, cLen);
                _alignment_d += std::string(cLen, '-');
                text_pos += cLen;
            }
        } else if (cLetter == 'D') {
            if (cLen >= 1) {
                _alignment_q += std::string(cLen, '-');
                _alignment_d += dna_d.substr(pattern_pos, cLen);
                pattern_pos += cLen;
            }
        }
    }

    int32_t totalScore = ez.score;
    free(ez.cigar);
    free(ts);
    free(qs);

    return totalScore;
}

int64_t alignSlidingWindow_minimap2_or_NW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                                          const int64_t &slidingWindowSize, const int32_t &matchingScore,
                                          const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {

    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns end
    int32_t longerSeqLength = std::max(_length_of_d, _length_of_q);

    if ((_length_of_d * 1.0 / slidingWindowSize) * (_length_of_q * 1.0 / slidingWindowSize) <= 1) { //  _length_of_d*_length_of_q <= (slidingWindowSize*slidingWindowSize) this calculated via RAM cost
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = -1;
        totalScore = alignSlidingWindow_minimap2(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                 mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
    } else if (longerSeqLength * 2.0 < slidingWindowSize) { // this calculated with RAM cost longerSeqLength*slidingWindowSize*2 <= (slidingWindowSize*slidingWindowSize
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = (slidingWindowSize * 0.5 / longerSeqLength) * slidingWindowSize;
        totalScore = alignSlidingWindow_minimap2(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                 mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
    } else {
        totalScore = alignSlidingWindow(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                        mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
    }

    return totalScore;
}

int64_t alignSlidingWindow_local_wfa2_v2(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                                         const int64_t &slidingWindowSize, const int32_t &matchingScore,
                                         const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";

    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns begin
    bool flag_q_all_N = std::all_of(dna_q.begin(), dna_q.end(), [](char c) { return c == 'N'; });
    bool flag_d_all_N = std::all_of(dna_d.begin(), dna_d.end(), [](char c) { return c == 'N'; });

    if (flag_q_all_N || flag_d_all_N) {
        _alignment_q = dna_q;
        _alignment_d = dna_d;

        int32_t count_ = abs(_length_of_q - _length_of_d);
        if (_length_of_q < _length_of_d) {
            _alignment_q += std::string(count_, '-');
        }

        if (_length_of_d < _length_of_q) {
            _alignment_d += std::string(count_, '-');
        }

        return totalScore;
    }

    double ratio = _length_of_q * 1.0 / _length_of_d;
    if (_length_of_q < _length_of_d)
        ratio = _length_of_d * 1.0 / _length_of_q;

    if (_length_of_q < 5 || _length_of_d < 5 || ratio > 3.0) {
//        std::cout << " minimap2 1 " << std::endl;
        totalScore = alignSlidingWindow_minimap2_or_NW(dna_q, dna_d, _alignment_q, _alignment_d,
                                                       slidingWindowSize, matchingScore,
                                                       mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                       openGapPenalty2,
                                                       extendGapPenalty2);
    } else {
        //check all Ns end
        std::string qSeq = dna_q;
        std::string dSeq = dna_d;
        const char *pattern = dna_d.c_str();
        const char *text = dna_q.c_str();
        wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
        attributes.distance_metric = gap_affine_2p;
        attributes.affine2p_penalties.mismatch = -mismatchingPenalty;       // X > 0
        attributes.affine2p_penalties.gap_opening1 = -openGapPenalty1;      // O1 >= 0
        attributes.affine2p_penalties.gap_extension1 = -extendGapPenalty1;  // E1 > 0
        attributes.affine2p_penalties.gap_opening2 = -openGapPenalty2;      // O2 >= 0
        attributes.affine2p_penalties.gap_extension2 = -extendGapPenalty2;  // E2 > 0
        attributes.alignment_form.span = alignment_end2end;
        attributes.alignment_scope = compute_alignment;
        attributes.memory_mode = wavefront_memory_high;
        attributes.heuristic.strategy = wf_heuristic_none;
        attributes.system.max_memory_abort = slidingWindowSize * slidingWindowSize;

        wavefront_aligner_t *const wf_aligner = wavefront_aligner_new(&attributes);

        if (WF_STATUS_SUCCESSFUL == wavefront_align(wf_aligner, pattern, strlen(pattern), text, strlen(text))) { // Align
            totalScore = wf_aligner->cigar->score;
            cigar_t *const edit_cigar = wf_aligner->cigar;
            char *const operations = edit_cigar->operations;

            int pattern_pos = 0, text_pos = 0;
            for (int i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i) {
                if (operations[i] == 'M' || operations[i] == 'X') {
                    _alignment_q += qSeq[text_pos];
                    _alignment_d += dSeq[pattern_pos];
                    pattern_pos++;
                    text_pos++;
                } else if (operations[i] == 'I') {
                    _alignment_q += qSeq[text_pos];
                    _alignment_d += '-';
                    text_pos++;
                } else if (operations[i] == 'D') {
                    _alignment_q += '-';
                    _alignment_d += dSeq[pattern_pos];
                    pattern_pos++;
                }
            }
            wavefront_aligner_delete(wf_aligner);
        } else {
            wavefront_aligner_delete(wf_aligner);
            totalScore = alignSlidingWindow_minimap2_or_NW(dna_q, dna_d, _alignment_q, _alignment_d,
                                                           slidingWindowSize, matchingScore,
                                                           mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                           openGapPenalty2,
                                                           extendGapPenalty2);
        }

    }
    return totalScore;
}

int64_t alignSlidingWindow(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int64_t &slidingWindowSize,  const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {
    return alignSlidingWindow_local_wfa2_v2(dna_q, dna_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                            mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
}

int64_t alignSlidingWindowNW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                             const int64_t &slidingWindowSize, const int32_t &matchingScore,
                             const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                             const Scorei &m) {
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns begin
    bool flag_q_all_N = std::all_of(dna_q.begin(), dna_q.end(), [](char c) { return c == 'N'; });
    bool flag_d_all_N = std::all_of(dna_d.begin(), dna_d.end(), [](char c) { return c == 'N'; });

    if (flag_q_all_N || flag_d_all_N) {
        _alignment_q = dna_q;
        _alignment_d = dna_d;

        int32_t count_ = abs(_length_of_q - _length_of_d);
        if (_length_of_q < _length_of_d) {
            _alignment_q += std::string(count_, '-');
        }

        if (_length_of_d < _length_of_q) {
            _alignment_d += std::string(count_, '-');
        }

        return totalScore;
    }

    totalScore = alignSlidingWindow(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                    mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

    return totalScore;
}
