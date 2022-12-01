//
// Created by song on 8/5/18.
//

#include "alignSlidingWindow.h"


int max(const int32_t &a, const int32_t &b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

int max(const int32_t &m, const int32_t &e1, const int32_t &f1, const int32_t &e2, const int32_t &f2, int8_t &figure) {
    int32_t max = m;
    figure = -1;
    if (f1 > max) {
        max = f1;
        figure = 3;
    }
    if (f2 > max) {
        max = f2;
        figure = 4;
    }
    if (e1 > max) {
        max = e1;
        figure = 1;
    }
    if (e2 > max) {
        max = e2;
        figure = 2;
    }
    return max;
}

int max(const int32_t &m, const int32_t &e1, const int32_t &f1, const int32_t &e2, const int32_t &f2) {
    int32_t max = m;
    if (e1 > max) {
        max = e1;
    }
    if (e2 > max) {
        max = e2;
    }
    if (f1 > max) {
        max = f1;
    }
    if (f2 > max) {
        max = f2;
    }
    return max;
}


/**
 * The matrix operation implemented in this source file is very complex to understand for human.
 * But is very friendly for parallel RMA operation on HPC.
 * If I created several standard 2D matrix, the parallel performance on cornell bioHPC was very bad.
 * */
int32_t needleAlignment(const std::string &_dna_q, const std::string &_dna_d, std::stack<char> &SQ, std::stack<char> &SD,
                        const int32_t &mismatchingPenalty, const int32_t &_open_gap_penalty1, const int32_t &_extend_gap_penalty1, const int32_t &_open_gap_penalty2, const int32_t &_extend_gap_penalty2) {
    int32_t matchingScore = 0;
    int32_t length1 = _dna_d.length();
    int32_t length2 = _dna_q.length();
    int64_t matrixsize = (length1 + 1) * (length2 + 1);

    int32_t *MM = new int32_t[matrixsize * 6];
    memset(MM, 0, matrixsize * 6);

    int8_t *TM = new int8_t[matrixsize * 6];
    memset(TM, 0, matrixsize * 2);              //V and M
    memset(TM + matrixsize * 2, 1, matrixsize); //E1
    memset(TM + matrixsize * 3, 3, matrixsize); //F1
    memset(TM + matrixsize * 4, 2, matrixsize); //E2
    memset(TM + matrixsize * 5, 4, matrixsize); //F2

    int32_t i = 0, j;
    for (j = 0; j < (length2 + 1); ++j) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }

    j = 0;
    for (i = 0; i < (length1 + 1); ++i) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }

    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;

    for (i = 1; i <= length1; ++i) {
        for (j = 1; j <= length2; ++j) {
            mScore = _dna_d[i - 1] == _dna_q[j - 1] ? matchingScore : mismatchingPenalty;
            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i - 1) * (length2 + 1) + j - 1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if (_extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            } else {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            } else {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if (_extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            } else {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            } else {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j],
                                                             TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }

    int32_t endPosition1 = length1;
    int32_t endPosition2 = length2;
    int32_t maxScore = MM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SD.push(_dna_d[i - 1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_q[j - 1]);
            SD.push('-');
            --j;
        } else {
            if (0 == trackMatrix) { // come from V
                {
                    if (TM[matrixsize * 0 + i * (length2 + 1) + j] == -1) {
                        SQ.push(_dna_q[j - 1]);
                        SD.push(_dna_d[i - 1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i;
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 3) {
                        SD.push(_dna_d[i - 1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 4) {
                        SD.push(_dna_d[i - 1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 1) {
                        SD.push('-');
                        SQ.push(_dna_q[j - 1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 2) {
                        SD.push('-');
                        SQ.push(_dna_q[j - 1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    } else {

                    }
                }
            } else if (trackMatrix == -1) {
                SQ.push(_dna_q[j - 1]);
                SD.push(_dna_d[i - 1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i;
                --j;
            } else if (1 == trackMatrix) {
                SD.push('-');
                SQ.push(_dna_q[j - 1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            } else if (2 == trackMatrix) {
                SD.push('-');
                SQ.push(_dna_q[j - 1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            } else if (3 == trackMatrix) {
                SD.push(_dna_d[i - 1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            } else if (4 == trackMatrix) {
                SD.push(_dna_d[i - 1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            } else {

            }
        }
    }

    delete[] MM;
    delete[] TM;

    return maxScore;
}

int32_t needleAlignment(const std::string &_dna_ref2, const std::string &_dna_query, const std::string &_dna_ref1, std::stack<char> &SQ, std::stack<char> &SR1, std::stack<char> &SR2,
                        const int32_t &mismatchingPenalty, const int32_t &_open_gap_penalty1, const int32_t &_extend_gap_penalty1, const int32_t &_open_gap_penalty2, const int32_t &_extend_gap_penalty2) {
    int32_t matchingScore = 0;
    int32_t length1 = _dna_ref1.length();
    int32_t length2 = _dna_ref2.length();
    int64_t matrixsize = (length1 + 1) * (length2 + 1);

    int32_t *MM = new int32_t[matrixsize * 6];
    memset(MM, 0, matrixsize * 6);

    int8_t *TM = new int8_t[matrixsize * 6];
    memset(TM, 0, matrixsize * 2);              //V and M
    memset(TM + matrixsize * 2, 1, matrixsize); //E1
    memset(TM + matrixsize * 3, 3, matrixsize); //F1
    memset(TM + matrixsize * 4, 2, matrixsize); //E2
    memset(TM + matrixsize * 5, 4, matrixsize); //F2

    int32_t i = 0, j;
    for (j = 0; j < (length2 + 1); ++j) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }

    j = 0;
    for (i = 0; i < (length1 + 1); ++i) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }

    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;

    for (i = 1; i <= length1; ++i) {
        for (j = 1; j <= length2; ++j) {
            mScore = 0;
            mScore += _dna_ref2[j - 1] == '-' ? 0 : (_dna_ref1[i - 1] == _dna_ref2[j - 1] ? matchingScore : mismatchingPenalty);
            mScore += _dna_query[j - 1] == '-' ? 0 : (_dna_ref1[i - 1] == _dna_query[j - 1] ? matchingScore : mismatchingPenalty);
            mScore = mScore / 2;

            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i - 1) * (length2 + 1) + j - 1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if (_extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            } else {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            } else {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if (_extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            } else {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            } else {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j],
                                                             TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }

    int32_t endPosition1 = length1;
    int32_t endPosition2 = length2;
    int32_t maxScore = MM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SR2.push('-');
            SR1.push(_dna_ref1[i - 1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_query[j - 1]);
            SR2.push(_dna_ref2[j - 1]);
            SR1.push('-');
            --j;
        } else {
            if (0 == trackMatrix) { // come from V
                {
                    if (TM[matrixsize * 0 + i * (length2 + 1) + j] == -1) {
                        SQ.push(_dna_query[j - 1]);
                        SR2.push(_dna_ref2[j - 1]);
                        SR1.push(_dna_ref1[i - 1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i;
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 3) {
                        SR1.push(_dna_ref1[i - 1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 4) {
                        SR1.push(_dna_ref1[i - 1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 1) {
                        SR1.push('-');
                        SQ.push(_dna_query[j - 1]);
                        SR2.push(_dna_ref2[j - 1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 2) {
                        SR1.push('-');
                        SQ.push(_dna_query[j - 1]);
                        SR2.push(_dna_ref2[j - 1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    } else {

                    }
                }
            } else if (trackMatrix == -1) {
                SQ.push(_dna_query[j - 1]);
                SR2.push(_dna_ref2[j - 1]);
                SR1.push(_dna_ref1[i - 1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i;
                --j;
            } else if (1 == trackMatrix) {
                SR1.push('-');
                SQ.push(_dna_query[j - 1]);
                SR2.push(_dna_ref2[j - 1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            } else if (2 == trackMatrix) {
                SR1.push('-');
                SQ.push(_dna_query[j - 1]);
                SR2.push(_dna_ref2[j - 1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            } else if (3 == trackMatrix) {
                SR1.push(_dna_ref1[i - 1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            } else if (4 == trackMatrix) {
                SR1.push(_dna_ref1[i - 1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            } else {

            }
        }
    }

    delete[] MM;
    delete[] TM;

    return maxScore;
}

int32_t needleAlignment(const std::vector<std::string> &_dna_refs, const std::vector<std::string> &_dna_queries, std::vector<std::stack<char>> &SQs, std::vector<std::stack<char>> &SRs,
                        const int32_t &mismatchingPenalty, const int32_t &_open_gap_penalty1, const int32_t &_extend_gap_penalty1, const int32_t &_open_gap_penalty2, const int32_t &_extend_gap_penalty2) {
    int32_t matchingScore = 0;
    int32_t length1 = _dna_refs[0].length();
    int32_t length2 = _dna_queries[0].length();
    int64_t matrixsize = (length1 + 1) * (length2 + 1);

    int32_t *MM = new int32_t[matrixsize * 6];
    memset(MM, 0, matrixsize * 6);

    int8_t *TM = new int8_t[matrixsize * 6];
    memset(TM,                  0, matrixsize * 2); //V and M
    memset(TM + matrixsize * 2, 1, matrixsize); //E1
    memset(TM + matrixsize * 3, 3, matrixsize); //F1
    memset(TM + matrixsize * 4, 2, matrixsize); //E2
    memset(TM + matrixsize * 5, 4, matrixsize); //F2

    int32_t i = 0, j;
    for (j = 0; j < (length2 + 1); ++j) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }

    j = 0;
    for (i = 0; i < (length1 + 1); ++i) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }

    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;

    for (i = 1; i <= length1; ++i) {
        for (j = 1; j <= length2; ++j) {
            mScore = 0;
            for (int k = 0; k < _dna_refs.size(); ++k) {
                for (int l = 0; l < _dna_queries.size(); ++l) {
                    mScore += _dna_refs[k][i - 1] == '-' ? 0 : _dna_queries[l][j - 1] == '-' ? 0 : (_dna_refs[k][i - 1] == _dna_queries[l][j - 1] ? matchingScore : mismatchingPenalty);
                }
            }
            int32_t sizes = (_dna_refs.size() * _dna_queries.size());
            if (i == j) {

            }

            mScore = mScore / sizes;

            if (i == j) {

            }

            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i - 1) * (length2 + 1) + j - 1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if (_extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            } else {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            } else {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if (_extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            } else {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            } else {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j],
                                                             TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }

    int32_t endPosition1 = length1;
    int32_t endPosition2 = length2;
    int32_t maxScore = MM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    while (i != 0 || j != 0) {
        if (j == 0) {
            for (int l = 0; l < _dna_queries.size(); ++l) {
                SQs[l].push('-');
            }
            for (int k = 0; k < _dna_refs.size(); ++k) {
                SRs[k].push(_dna_refs[k][i - 1]);
            }
            --i;
        } else if (i == 0) {
            for (int l = 0; l < _dna_queries.size(); ++l) {
                SQs[l].push(_dna_queries[l][j - 1]);
            }
            for (int k = 0; k < _dna_refs.size(); ++k) {
                SRs[k].push('-');
            }
            --j;
        } else {
            if (0 == trackMatrix) { // come from V
                {
                    if (TM[matrixsize * 0 + i * (length2 + 1) + j] == -1) {
                        for (int l = 0; l < _dna_queries.size(); ++l) {
                            SQs[l].push(_dna_queries[l][j - 1]);
                        }
                        for (int k = 0; k < _dna_refs.size(); ++k) {
                            SRs[k].push(_dna_refs[k][i - 1]);
                        }
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i;
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 3) {
                        for (int k = 0; k < _dna_refs.size(); ++k) {
                            SRs[k].push(_dna_refs[k][i - 1]);
                        }
                        for (int l = 0; l < _dna_queries.size(); ++l) {
                            SQs[l].push('-');
                        }
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 4) {
                        for (int k = 0; k < _dna_refs.size(); ++k) {
                            SRs[k].push(_dna_refs[k][i - 1]);
                        }
                        for (int l = 0; l < _dna_queries.size(); ++l) {
                            SQs[l].push('-');
                        }
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 1) {
                        for (int k = 0; k < _dna_refs.size(); ++k) {
                            SRs[k].push('-');
                        }
                        for (int l = 0; l < _dna_queries.size(); ++l) {
                            SQs[l].push(_dna_queries[l][j - 1]);
                        }
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 2) {
                        for (int k = 0; k < _dna_refs.size(); ++k) {
                            SRs[k].push('-');
                        }
                        for (int l = 0; l < _dna_queries.size(); ++l) {
                            SQs[l].push(_dna_queries[l][j - 1]);
                        }
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    } else {

                    }
                }
            } else if (trackMatrix == -1) {
                for (int l = 0; l < _dna_queries.size(); ++l) {
                    SQs[l].push(_dna_queries[l][j - 1]);
                }
                for (int k = 0; k < _dna_refs.size(); ++k) {
                    SRs[k].push(_dna_refs[k][i - 1]);
                }
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i;
                --j;
            } else if (1 == trackMatrix) {
                for (int k = 0; k < _dna_refs.size(); ++k) {
                    SRs[k].push('-');
                }
                for (int l = 0; l < _dna_queries.size(); ++l) {
                    SQs[l].push(_dna_queries[l][j - 1]);
                }
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            } else if (2 == trackMatrix) {
                for (int k = 0; k < _dna_refs.size(); ++k) {
                    SRs[k].push('-');
                }
                for (int l = 0; l < _dna_queries.size(); ++l) {
                    SQs[l].push(_dna_queries[l][j - 1]);
                }
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            } else if (3 == trackMatrix) {
                for (int k = 0; k < _dna_refs.size(); ++k) {
                    SRs[k].push(_dna_refs[k][i - 1]);
                }
                for (int l = 0; l < _dna_queries.size(); ++l) {
                    SQs[l].push('-');
                }
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            } else if (4 == trackMatrix) {
                for (int k = 0; k < _dna_refs.size(); ++k) {
                    SRs[k].push(_dna_refs[k][i - 1]);
                }
                for (int l = 0; l < _dna_queries.size(); ++l) {
                    SQs[l].push('-');
                }
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            } else {

            }
        }
    }

    delete[] MM;
    delete[] TM;

    return maxScore;
}

int32_t alignment(const std::string &_dna_q, const std::string &_dna_d, std::stack<char> &SQ, std::stack<char> &SD,
                  const int32_t &matchingScore, int32_t mismatchingPenalty,
                  int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                  int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2) {
    int32_t length1 = _dna_d.length();
    int32_t length2 = _dna_q.length();

    int64_t matrixsize = (length1 + 1) * (length2 + 1);

    int32_t *MM = new int32_t[matrixsize * 6];
    memset(MM, 0, matrixsize * 6);

    int8_t *TM = new int8_t[matrixsize * 6];
    memset(TM, 0, matrixsize * 2);              //V and M
    memset(TM + matrixsize * 2, 1, matrixsize); //E1
    memset(TM + matrixsize * 3, 3, matrixsize); //F1
    memset(TM + matrixsize * 4, 2, matrixsize); //E2
    memset(TM + matrixsize * 5, 4, matrixsize); //F2

    int32_t i = 0, j;
    for (j = 0; j < (length2 + 1); ++j) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }

    j = 0;
    for (i = 0; i < (length1 + 1); ++i) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }

    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;

    for (i = 1; i <= length1; ++i) {
        for (j = 1; j <= length2; ++j) {
            mScore = _dna_d[i - 1] == _dna_q[j - 1] ? matchingScore : mismatchingPenalty;
            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i - 1) * (length2 + 1) + j - 1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if (_extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            } else {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            } else {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if (_extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            } else {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            } else {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j],
                                                             TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }

    int32_t endPosition1 = 0;
    int32_t endPosition2 = 0;
    int32_t maxScore = -1000000;

    i = length1;
    for (j = 1; j <= length2; ++j) {
        if (MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore) { // please do not change >= to >, since we are doing global alignment
            // >= will omit the first similar fragments
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2 = j;
        }
    }

    j = length2;
    for (i = 1; i <= length1; ++i) {
        if (MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore) { // please do not change >= to >, since we are doing global alignment
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2 = j;
        }
    }

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SD.push(_dna_d[i - 1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_q[j - 1]);
            SD.push('-');
            --j;
        } else {
            if (0 == trackMatrix) { // come from V
                {
                    if (TM[matrixsize * 0 + i * (length2 + 1) + j] == -1) {
                        SQ.push(_dna_q[j - 1]);
                        SD.push(_dna_d[i - 1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i;
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 3) {
                        SD.push(_dna_d[i - 1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 4) {
                        SD.push(_dna_d[i - 1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 1) {
                        SD.push('-');
                        SQ.push(_dna_q[j - 1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 2) {
                        SD.push('-');
                        SQ.push(_dna_q[j - 1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    } else {

                    }
                }
            } else if (trackMatrix == -1) {
                SQ.push(_dna_q[j - 1]);
                SD.push(_dna_d[i - 1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i;
                --j;
            } else if (1 == trackMatrix) {
                SD.push('-');
                SQ.push(_dna_q[j - 1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            } else if (2 == trackMatrix) {
                SD.push('-');
                SQ.push(_dna_q[j - 1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            } else if (3 == trackMatrix) {
                SD.push(_dna_d[i - 1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            } else if (4 == trackMatrix) {
                SD.push(_dna_d[i - 1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            } else {

            }
        }
    }

    delete[] MM;
    delete[] TM;

    return maxScore;
}

int32_t alignment(const std::string &_dna_ref2, const std::string &_dna_query, const std::string &_dna_ref1, std::stack<char> &SQ, std::stack<char> &SR1, std::stack<char> &SR2, const int32_t &matchingScore,
                  const int32_t &mismatchingPenalty, const int32_t &_open_gap_penalty1, const int32_t &_extend_gap_penalty1, const int32_t &_open_gap_penalty2, const int32_t &_extend_gap_penalty2) {
    int32_t length1 = _dna_ref1.length();
    int32_t length2 = _dna_ref2.length();

    int64_t matrixsize = (length1 + 1) * (length2 + 1);
    int32_t *MM = new int32_t[matrixsize * 6];
    memset(MM, 0, matrixsize * 6);

    int8_t *TM = new int8_t[matrixsize * 6];
    memset(TM, 0, matrixsize * 2);              //V and M
    memset(TM + matrixsize * 2, 1, matrixsize); //E1
    memset(TM + matrixsize * 3, 3, matrixsize); //F1
    memset(TM + matrixsize * 4, 2, matrixsize); //E2
    memset(TM + matrixsize * 5, 4, matrixsize); //F2

    int32_t i = 0, j;
    for (j = 0; j < (length2 + 1); ++j) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j * _extend_gap_penalty1, _open_gap_penalty2 + j * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }

    j = 0;
    for (i = 0; i < (length1 + 1); ++i) {
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i * _extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i * _extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }

    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;

    for (i = 1; i <= length1; ++i) {
        for (j = 1; j <= length2; ++j) {
            mScore = 0;
            mScore += _dna_ref2[j - 1] == '-' ? 0 : (_dna_ref1[i - 1] == _dna_ref2[j - 1] ? matchingScore : mismatchingPenalty);
            mScore += _dna_query[j - 1] == '-' ? 0 : (_dna_ref1[i - 1] == _dna_query[j - 1] ? matchingScore : mismatchingPenalty);
            mScore = mScore / 2;

            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i - 1) * (length2 + 1) + j - 1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if (_extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            } else {
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1]) {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            } else {
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j - 1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if (_extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            } else {
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if (_extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j]) {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            } else {
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i - 1) * (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j],
                                                             TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }

    int32_t endPosition1 = 0;
    int32_t endPosition2 = 0;
    int32_t maxScore = -1000000;

    i = length1;
    for (j = 1; j <= length2; ++j) {
        if (MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore) { // please do not change >= to >, since we are doing global alignment
            // >= will omit the first similar fragments
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2 = j;
        }
    }

    j = length2;
    for (i = 1; i <= length1; ++i) {
        if (MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore) { // please do not change >= to >, since we are doing global alignment
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2 = j;
        }
    }

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SR2.push('-');
            SR1.push(_dna_ref1[i - 1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_query[j - 1]);
            SR2.push(_dna_ref2[j - 1]);
            SR1.push('-');
            --j;
        } else {
            if (0 == trackMatrix) { // come from V
                {
                    if (TM[matrixsize * 0 + i * (length2 + 1) + j] == -1) {
                        SQ.push(_dna_query[j - 1]);
                        SR2.push(_dna_ref2[j - 1]);
                        SR1.push(_dna_ref1[i - 1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i;
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 3) {
                        SR1.push(_dna_ref1[i - 1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 4) {
                        SR1.push(_dna_ref1[i - 1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 1) {
                        SR1.push('-');
                        SQ.push(_dna_query[j - 1]);
                        SR2.push(_dna_ref2[j - 1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (TM[matrixsize * 0 + i * (length2 + 1) + j] == 2) {
                        SR1.push('-');
                        SQ.push(_dna_query[j - 1]);
                        SR2.push(_dna_ref2[j - 1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    } else {

                    }
                }
            } else if (trackMatrix == -1) {
                SQ.push(_dna_query[j - 1]);
                SR2.push(_dna_ref2[j - 1]);
                SR1.push(_dna_ref1[i - 1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i;
                --j;
            } else if (1 == trackMatrix) {
                SR1.push('-');
                SQ.push(_dna_query[j - 1]);
                SR2.push(_dna_ref2[j - 1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            } else if (2 == trackMatrix) {
                SR1.push('-');
                SQ.push(_dna_query[j - 1]);
                SR2.push(_dna_ref2[j - 1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            } else if (3 == trackMatrix) {
                SR1.push(_dna_ref1[i - 1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            } else if (4 == trackMatrix) {
                SR1.push(_dna_ref1[i - 1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            } else {

            }
        }
    }

    delete[] MM;
    delete[] TM;

    return maxScore;
}

int32_t alignment_position(const std::string &_dna_q, const std::string &_dna_d, std::stack<char> &SQ, std::stack<char> &SD,
                           const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                           const int32_t &_open_gap_penalty1, const int32_t &_extend_gap_penalty1,
                           const int32_t &_open_gap_penalty2, const int32_t &_extend_gap_penalty2) { // far from done


    int32_t length1 = _dna_d.length(), length2 = _dna_q.length(), i, j, f1, f2, mScore, endPosition1 = 0, endPosition2 = 0, maxScore = INT32_MIN;

    int32_t *MP = new int32_t[length1 + 1];
    int32_t *MC = new int32_t[length1 + 1];
    int32_t *E1 = new int32_t[length1 + 1];
    int32_t *E2 = new int32_t[length1 + 1];

    memset(MP, 0, length1 + 1);
    memset(MC, 0, length1 + 1);
    memset(E1, 0, length1 + 1);
    memset(E2, 0, length1 + 1);


    for (i = 0; i < (length1 + 1); ++i) {
        MP[i] = max(_open_gap_penalty1 + i * _extend_gap_penalty1, _open_gap_penalty2 + i * _extend_gap_penalty2);
        E1[i] = _open_gap_penalty1;
        E2[i] = _open_gap_penalty2;
    }

    MP[0] = 0;
    MC[0] = 0;
    for (j = 1; j <= length2; ++j) {
        f1 = _open_gap_penalty1;
        f2 = _open_gap_penalty2;
        for (i = 1; i <= length1; ++i) {

            f1 = (_open_gap_penalty1 + _extend_gap_penalty1 + MP[j - 1]) > (_extend_gap_penalty1 + f1) ? (_open_gap_penalty1 + _extend_gap_penalty1 + MP[j - 1]) : (_extend_gap_penalty1 + f1);
            f2 = (_open_gap_penalty2 + _extend_gap_penalty2 + MP[j - 1]) > (_extend_gap_penalty2 + f2) ? (_open_gap_penalty2 + _extend_gap_penalty2 + MP[j - 1]) : (_extend_gap_penalty2 + f2);

            mScore = _dna_d[i - 1] == _dna_q[j - 1] ? matchingScore : mismatchingPenalty;
            MC[i] = mScore + MP[i - 1];

            if (_extend_gap_penalty1 + E1[i] >= _open_gap_penalty1 + _extend_gap_penalty1 + MP[i]) {
                E1[i] = _extend_gap_penalty1 + E1[i];
            } else {
                E1[i] = _open_gap_penalty1 + _extend_gap_penalty1 + MP[i];
            }
            if (_extend_gap_penalty2 + E2[i] >= _open_gap_penalty2 + _extend_gap_penalty2 + MP[i]) {
                E2[i] = _extend_gap_penalty2 + E2[i];
            } else {
                E2[i] = _open_gap_penalty2 + _extend_gap_penalty2 + MP[i];
            }
            MC[i] = max(MC[i], E1[i], f1, E2[i], f2);
        }
        if (MP[length1] < maxScore) {
            endPosition1 = length1;
            endPosition2 = j;
        }
        int32_t *TEMP = MP;
        MP = MC;
        MC = TEMP;
    }

    for (j = 1; j <= length2; ++j) {
        if (MP[length1] < maxScore) {
            endPosition1 = length1;
        }
    }

    return maxScore;
}

int32_t minimap2_alignment(const std::string& _dna_q, const std::string& _dna_d, std::string & _alignment_q, std::string & _alignment_d,
                           const int32_t & matchingScore, int32_t mismatchingPenalty,
                           int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2){
    int8_t a = 0, b = mismatchingPenalty < 0? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = _dna_d.length(), ql = _dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = _dna_d.c_str();
    const char *qseq = _dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    //    std::cout << "line 359" << std::endl;
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    //    std::cout << "line 362" << std::endl;
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    //    std::cout << "line 365" << std::endl;
    int i;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t)qseq[i]];
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
        ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0, & ez);
#endif
    std::string cigarstring = "";
    for (i = 0; i < ez.n_cigar; ++i){ // print CIGAR
        //        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
        cigarstring = cigarstring + std::to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];
    }
    //    putchar('\n');

    std::vector<std::string> cigarElems;
    splitCIGAR(cigarstring, cigarElems);

    _alignment_q = "";
    _alignment_d = "";
    int pattern_pos = 0, text_pos = 0;

    //    std::cout << cigarstring << std::endl;
    //    std::cout << std::to_string(cigarElems.size()) << std::endl;

    for(i=0; i<cigarElems.size(); ++i) {
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
    free(ez.cigar); free(ts); free(qs);
    return totalScore;
}



int32_t alignment_minimap2(const std::string& _dna_q, const std::string& _dna_d, std::string & _alignment_q, std::string & _alignment_d,
                           const int32_t & matchingScore, int32_t mismatchingPenalty, int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                           int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2, int32_t & endPositionq, int32_t & endPositiont){

    int8_t a = 0, b = mismatchingPenalty < 0? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = _dna_d.length(), ql = _dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = _dna_d.c_str();
    const char *qseq = _dna_q.c_str();

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    //    std::cout << "line 359" << std::endl;
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    //    std::cout << "line 362" << std::endl;
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    //    std::cout << "line 365" << std::endl;
    int i;
    for (i = 0; i < tl; ++i) {
        ts[i] = c[(uint8_t) tseq[i]]; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = c[(uint8_t)qseq[i]];
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
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -_open_gap_penalty1, -_extend_gap_penalty1, -_open_gap_penalty2, -_extend_gap_penalty2, -1, -1, 0, 0x01, & ez);
#endif


    endPositionq=0;
    endPositiont=0;
    int32_t maxScore = -1000000;

    if( ez.mqe > ez.mte ){
        endPositiont = ez.mqe_t + 1;
        endPositionq = ql;
        maxScore = ez.mqe;
    }else{
        endPositiont = tl;
        endPositionq = ez.mte_q + 1;
        maxScore = ez.mte;
    }
    free(ez.cigar); free(ts); free(qs);
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

    if (sqrt(_length_of_d) * sqrt(_length_of_q) <= slidingWindowSize) {
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
//            std::stack<char> SQ;
//            std::stack<char> SD;
            std::string alignment_q="";
            std::string alignment_d="";
            if (slidingWindowSize > 1073741824) {
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            } else {
                if( extendGapPenalty2+matchingScore < 0){
                    totalScore += alignment_minimap2(qSeq, dSeq, alignment_q, alignment_d, matchingScore, mismatchingPenalty+matchingScore, openGapPenalty1+matchingScore, extendGapPenalty1+matchingScore, openGapPenalty2+matchingScore, extendGapPenalty2+matchingScore, endPositionq, endPositiont);
                }else{
                    totalScore += alignment_minimap2(qSeq, dSeq, alignment_q, alignment_d, matchingScore, mismatchingPenalty+matchingScore-1, openGapPenalty1+matchingScore-1, extendGapPenalty1+matchingScore-1, openGapPenalty2+matchingScore-1, -1, endPositionq, endPositiont);
                }
            }

            std::string tempalignment_q = alignment_q;
            tempalignment_q = songStrReplaceAll(tempalignment_q, "-", "");
            assert(endPositionq == tempalignment_q.size());
            std::string tempalignment_d = alignment_d;
            tempalignment_d = songStrReplaceAll(tempalignment_d, "-", "");
            assert(endPositiont == tempalignment_d.size());

            _alignment_q += alignment_q;
            _alignment_d += alignment_d;
            queryStart += endPositionq;
            databaseStart += endPositiont;
        }
    }

    int32_t final_indel_length = 0;
    while( databaseStart<=_length_of_d ){
        _alignment_q += '-';
        _alignment_d += dna_d[databaseStart-1];
        ++databaseStart;
        ++final_indel_length;
    }
    assert(_alignment_d.size() == _alignment_q.size());
    while( queryStart<=_length_of_q ){
        _alignment_q += dna_q[queryStart-1];
        _alignment_d += '-';
        ++queryStart;
        ++final_indel_length;
    }
    assert(_alignment_d.size() == _alignment_q.size());
    if( final_indel_length > 0 ){
        totalScore += max(openGapPenalty1 + extendGapPenalty1*final_indel_length, openGapPenalty2+extendGapPenalty2*final_indel_length);
    }
    return totalScore;
}

int64_t alignSlidingWindow(std::string &align_ref2, std::string &align_query, const std::string &dna_ref1,
                           std::string &align_ref1, const int64_t &slidingWindowSize,
                           const int32_t &matchingScore, const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1,
                           const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2) {

    int32_t databaseStart = 1;
    int32_t databaseEnd = 0;
    int32_t queryStart = 1;
    int32_t queryEnd = 0;
    int64_t totalScore = 0;
    int64_t _length_of_q = align_ref2.size();
    int64_t _length_of_d = dna_ref1.size();
    std::string dna_ref2 = align_ref2;
    std::string dna_query = align_query;
    align_ref2 = "";
    align_query = "";

    if (sqrt(_length_of_d) * sqrt(_length_of_q) <= slidingWindowSize) {
        std::stack<char> SQ;
        std::stack<char> SR1;
        std::stack<char> SR2;
        totalScore += needleAlignment(dna_ref2, dna_query, dna_ref1, SQ, SR1, SR2, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        while (!SQ.empty()) {
            align_query += SQ.top();
            align_ref1 += SR1.top();
            align_ref2 += SR2.top();
            if (SQ.top() != '-' || SR2.top() != '-') {
                ++queryStart;
            }
            if (SR1.top() != '-') {
                ++databaseStart;
            }
            SQ.pop();
            SR2.pop();
            SR1.pop();
        }
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

            std::string qSeq = getSubsequence(dna_query, queryStart, queryEnd);
            std::string d2Seq = getSubsequence(dna_ref2, queryStart, queryEnd);
            std::string d1Seq = getSubsequence(dna_ref1, databaseStart, databaseEnd);

            std::stack<char> SQ;
            std::stack<char> SR1;
            std::stack<char> SR2;
            if (slidingWindowSize > 1073741824) {
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            } else {
                if (extendGapPenalty2 + matchingScore < 0) {
                    totalScore += alignment(d2Seq, qSeq, d1Seq, SQ, SR1, SR2, matchingScore, mismatchingPenalty + matchingScore, openGapPenalty1 + matchingScore, extendGapPenalty1 + matchingScore, openGapPenalty2 + matchingScore, extendGapPenalty2 + matchingScore);
                } else {
                    totalScore += alignment(d2Seq, qSeq, d1Seq, SQ, SR1, SR2, matchingScore, mismatchingPenalty + matchingScore - 1, openGapPenalty1 + matchingScore - 1, extendGapPenalty1 + matchingScore - 1, openGapPenalty2 + matchingScore - 1, -1);
                }
            }
            while (!SQ.empty()) {
                align_query += SQ.top();
                align_ref1 += SR1.top();
                align_ref2 += SR2.top();
                if (SQ.top() != '-' || SR2.top() != '-') {
                    ++queryStart;
                }
                if (SR1.top() != '-') {
                    ++databaseStart;
                }
                SQ.pop();
                SR2.pop();
                SR1.pop();
            }
        }
    }

    int32_t final_indel_length = 0;

    int32_t count_1 = _length_of_d - databaseStart;
    if (count_1 >= 0) {
        align_query += std::string(count_1, '-');
        align_ref2 += std::string(count_1, '-');
        align_ref1 += dna_ref1.substr(databaseStart - 1, count_1 + 1);
        final_indel_length += count_1;
    }

    int32_t count_2 = _length_of_q - queryStart;
    if (count_2 >= 0) {
        align_query += dna_query.substr(queryStart - 1, count_2 + 1);;
        align_ref2 += dna_ref2.substr(queryStart - 1, count_2 + 1);
        align_ref1 += std::string(count_2, '-');
        final_indel_length += count_2;
    }
    if (final_indel_length > 0) {
        totalScore += max(openGapPenalty1 + extendGapPenalty1 * final_indel_length, openGapPenalty2 + extendGapPenalty2 * final_indel_length);
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
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, & ez);
#endif

    std::string cigarstring = "";
    for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
        cigarstring = cigarstring + std::to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];
    }

    std::vector<std::string> cigarElems;
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

std::string GetStdoutFromCommand(std::string &cmd) {

    std::string data;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1");

    FILE *stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
                data.append(buffer);
        pclose(stream);
    }

    rtrim(data);
    return data;
}

int64_t alignSlidingWindow_minimap2_or_NW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                                          const int64_t &slidingWindowSize, const int32_t &wfaSize, const int32_t &matchingScore,
                                          const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                          const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, const Scorei &m) {
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int32_t _length_of_q = dna_q.size();
    int32_t _length_of_d = dna_d.size();

    //check all Ns end
    int32_t longerSeqLength = max(_length_of_d, _length_of_q);

    if ((_length_of_d*1.0/slidingWindowSize) * (_length_of_q*1.0/slidingWindowSize) <= 1) { //  _length_of_d*_length_of_q <= (slidingWindowSize*slidingWindowSize) this calculated via RAM cost
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = -1;
        totalScore = alignSlidingWindow_minimap2(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                 mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
    }
    else if (longerSeqLength*2.0 < slidingWindowSize) { // this calculated with RAM cost longerSeqLength*slidingWindowSize*2 <= (slidingWindowSize*slidingWindowSize
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = (slidingWindowSize * 0.5 / longerSeqLength) * slidingWindowSize;
        totalScore = alignSlidingWindow_minimap2(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                 mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
    }
    else {
        totalScore = alignSlidingWindow(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                        mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
    }

    return totalScore;
}

int64_t alignSlidingWindow_local_wfa2_v2(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                                         const int64_t &slidingWindowSize, const int32_t &wfaSize, const int32_t &matchingScore,
                                         const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                         const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, const Scorei &m) {
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
        if(_length_of_q < _length_of_d) {
            _alignment_q += std::string(count_, '-');
        }

        if(_length_of_d < _length_of_q) {
            _alignment_d += std::string(count_, '-');
        }

        return totalScore;
    }

    double ratio = _length_of_q * 1.0 / _length_of_d;
    if(_length_of_q < _length_of_d) 
        ratio = _length_of_d * 1.0 / _length_of_q;
    
    if (_length_of_q < 5 || _length_of_d < 5 || ratio > 3.0) {
//        std::cout << " minimap2 1 " << std::endl;
        totalScore = alignSlidingWindow_minimap2_or_NW(dna_q, dna_d, _alignment_q, _alignment_d,
                                                       slidingWindowSize, wfaSize, matchingScore,
                                                       mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                       openGapPenalty2,
                                                       extendGapPenalty2, min_wavefront_length,
                                                       max_distance_threshold, m);
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
        } else {
            totalScore = alignSlidingWindow_minimap2_or_NW(dna_q, dna_d, _alignment_q, _alignment_d,
                                                        slidingWindowSize, wfaSize, matchingScore,
                                                        mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                                        openGapPenalty2,
                                                        extendGapPenalty2, min_wavefront_length,
                                                        max_distance_threshold, m);
        }

        wavefront_aligner_delete(wf_aligner);
    }
    return totalScore;
}

int64_t alignSlidingWindow(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int64_t &slidingWindowSize, const int32_t &wfaSize, const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                           const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, const Scorei &m, std::map<std::string, std::string> &parameters) {
    return alignSlidingWindow_local_wfa2_v2(dna_q, dna_d, _alignment_q, _alignment_d, slidingWindowSize, wfaSize, matchingScore,
                                            mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                                            min_wavefront_length, max_distance_threshold, m);
}

int64_t alignSlidingWindowNW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                             const int64_t &slidingWindowSize, const int32_t &wfaSize, const int32_t &matchingScore,
                             const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                             const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, const Scorei &m) {
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
        if(_length_of_q < _length_of_d) {
            _alignment_q += std::string(count_, '-');
        }

        if(_length_of_d < _length_of_q) {
            _alignment_d += std::string(count_, '-');
        }

        return totalScore;
    }

    totalScore = alignSlidingWindow(dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                    mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

    return totalScore;
}

int64_t alignSlidingWindow(std::string &dna_ref1, std::string &dna_ref2, std::string &dna_query, std::string &_alignment_ref1, std::string &_alignment_ref2, std::string &_alignment_query,
                           const int64_t &slidingWindowSize, const int32_t &wfaSize, const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                           const int32_t &min_wavefront_length, const int32_t &max_distance_threshold, const Scorei &m, std::map<std::string, std::string> &parameters) {

    return alignSlidingWindow_local_wfa2_v2(dna_query, dna_ref2, _alignment_query, _alignment_ref2, slidingWindowSize, wfaSize, matchingScore,
                                            mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                                            min_wavefront_length, max_distance_threshold, m);
}
