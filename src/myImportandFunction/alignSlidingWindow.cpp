//
// Created by song on 8/5/18.
//

#include "alignSlidingWindow.h"


int max( const int32_t & a, const int32_t & b){
    if( a > b){
        return a;
    }else{
        return b;
    }
}
int max( const int32_t & m, const int32_t & e1, const int32_t & f1, const int32_t & e2, const int32_t & f2, int8_t & figure) {
    int32_t max = m;
    figure = -1;
    if ( f1 > max ){
        max = f1;
        figure = 3;
    }
    if (f2 > max ) {
        max = f2;
        figure = 4;
    }
    if (e1 > max) {
        max = e1;
        figure = 1;
    }
    if( e2>max ){
        max = e2;
        figure = 2;
    }
    return  max ;
}

int max( const int32_t & m, const int32_t & e1, const int32_t & f1, const int32_t & e2, const int32_t & f2) {
    int32_t max = m;
    if( e1 > max ){
        max = e1;
    }
    if( e2 > max ){
        max = e2;
    }
    if( f1 > max ){
        max = f1;
    }
    if( f2 > max ){
        max = f2;
    }
    return  max ;
}


/**
 * The matrix operation implemented in this source file is very complex to understand for human.
 * But is very friendly for parallel RMA operation on HPC.
 * If I created several standard 2D matrix, the parallel performance on cornell bioHPC was very bad.
 * */







int32_t needleAlignment(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                        const int32_t & mismatchingPenalty, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2){
    int32_t matchingScore = 0;
    int32_t length1 = _dna_d.length();
    int32_t length2 = _dna_q.length();
//    std::cout << "line 51" << std::endl;
    int64_t matrixsize = (length1 +1) * (length2+1);
  //  std::cout << "line 53" << std::endl;
    int32_t * MM = new int32_t [matrixsize * 6];
    int8_t * TM = new int8_t [matrixsize * 6];
    //std::cout << "line 56" << std::endl;
    std::fill_n(MM, matrixsize * 6, 0);
    //std::cout << "line 58" << std::endl;
    std::fill_n(TM+matrixsize*0, matrixsize*2, 0); //V and M
    std::fill_n(TM+matrixsize*2, matrixsize, 1); //E1
    std::fill_n(TM+matrixsize*3, matrixsize, 3); //F1
    std::fill_n(TM+matrixsize*4, matrixsize, 2); //E2
    std::fill_n(TM+matrixsize*5, matrixsize, 4); //F2

//    std::cout << "line 65" << std::endl;
    int32_t i=0, j;
    for (j = 0; j < (length2 + 1); ++j  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }
//    std::cout << "line 81" << std::endl;
    j=0;
    for (i = 0; i < (length1 + 1); ++i  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }
//    std::cout << "line 97" << std::endl;
    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;
//    std::cout << "line 106" << std::endl;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            mScore = _dna_d[i-1] == _dna_q[j-1] ? matchingScore : mismatchingPenalty;
            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i-1) * (length2 + 1) + j-1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if( _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            }else{
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if( _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            }else{
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if(_extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            }else{
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if(_extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            }else{
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j], TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }
//    std::cout << "line 144" << std::endl;

    int32_t endPosition1 = length1;
    int32_t endPosition2 = length2;
    int32_t maxScore = MM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];
//    std::cout << "line 154" << std::endl;
    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SD.push(_dna_d[i-1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_q[j-1]);
            SD.push('-');
            --j;
        } else {
            if( 0 == trackMatrix ) { // come from V
                {
                    if ( TM[matrixsize * 0 + i * (length2 + 1) + j] == -1 ) {
                        SQ.push(_dna_q[j-1]);
                        SD.push(_dna_d[i-1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i; --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 3 ) {
                        SD.push(_dna_d[i-1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 4 ) {
                        SD.push(_dna_d[i-1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 1 ) {
                        SD.push('-');
                        SQ.push(_dna_q[j-1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 2 ){
                        SD.push('-');
                        SQ.push(_dna_q[j-1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    }else{
                        std::cout << "line 203" << std::endl;
                    }
                }
            } else if ( trackMatrix == -1 ) {
                SQ.push(_dna_q[j-1]);
                SD.push(_dna_d[i-1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i; --j;
            } else if( 1 == trackMatrix ){
                SD.push('-');
                SQ.push(_dna_q[j-1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            }else if( 2 == trackMatrix ){
                SD.push('-');
                SQ.push(_dna_q[j-1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            }else if( 3 == trackMatrix ){
                SD.push(_dna_d[i-1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            }else if( 4 == trackMatrix ){
                SD.push(_dna_d[i-1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            }else{
                std::cout << "line 232" << std::endl;
            }
        }
    }
//    std::cout << "line 226" << std::endl;
    delete[] MM;
  //  std::cout << "line 228" << std::endl;
    delete[] TM;
    //std::cout << "line 229" << std::endl;
    return maxScore;
}





int32_t needleAlignment(const std::string& _dna_ref2, const std::string & _dna_query, const std::string& _dna_ref1, std::stack<char> & SQ, std::stack<char> & SR1, std::stack<char> & SR2,
                        const int32_t & mismatchingPenalty, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2){
    int32_t matchingScore = 0;
    int32_t length1 = _dna_ref1.length();
    int32_t length2 = _dna_ref2.length();
//    std::cout << "line 51" << std::endl;
    int64_t matrixsize = (length1 +1) * (length2+1);
    //  std::cout << "line 53" << std::endl;
    int32_t * MM = new int32_t [matrixsize * 6];
    int8_t * TM = new int8_t [matrixsize * 6];
    //std::cout << "line 56" << std::endl;
    std::fill_n(MM, matrixsize * 6, 0);
    //std::cout << "line 58" << std::endl;
    std::fill_n(TM+matrixsize*0, matrixsize*2, 0); //V and M
    std::fill_n(TM+matrixsize*2, matrixsize, 1); //E1
    std::fill_n(TM+matrixsize*3, matrixsize, 3); //F1
    std::fill_n(TM+matrixsize*4, matrixsize, 2); //E2
    std::fill_n(TM+matrixsize*5, matrixsize, 4); //F2

//    std::cout << "line 65" << std::endl;
    int32_t i=0, j;
    for (j = 0; j < (length2 + 1); ++j  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }
//    std::cout << "line 81" << std::endl;
    j=0;
    for (i = 0; i < (length1 + 1); ++i  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }
//    std::cout << "line 97" << std::endl;
    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;
//    std::cout << "line 106" << std::endl;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            mScore = 0;
            mScore += _dna_ref2[j-1] == '-' ? 0 : (_dna_ref1[i-1] == _dna_ref2[j-1] ? matchingScore : mismatchingPenalty);
            mScore += _dna_query[j-1] == '-' ? 0 : (_dna_ref1[i-1] == _dna_query[j-1] ? matchingScore : mismatchingPenalty);
            mScore = mScore /2;

            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i-1) * (length2 + 1) + j-1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if( _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            }else{
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if( _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            }else{
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if(_extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            }else{
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if(_extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            }else{
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j], TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }
//    std::cout << "line 144" << std::endl;

    int32_t endPosition1 = length1;
    int32_t endPosition2 = length2;
    int32_t maxScore = MM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];

    i = endPosition1;
    j = endPosition2;

    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];
//    std::cout << "line 154" << std::endl;
    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SR2.push('-');
            SR1.push(_dna_ref1[i-1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_query[j-1]);
            SR2.push(_dna_ref2[j-1]);
            SR1.push('-');
            --j;
        } else {
            if( 0 == trackMatrix ) { // come from V
                {
                    if ( TM[matrixsize * 0 + i * (length2 + 1) + j] == -1 ) {
                        SQ.push(_dna_query[j-1]);
                        SR2.push(_dna_ref2[j-1]);
                        SR1.push(_dna_ref1[i-1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i; --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 3 ) {
                        SR1.push(_dna_ref1[i-1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 4 ) {
                        SR1.push(_dna_ref1[i-1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 1 ) {
                        SR1.push('-');
                        SQ.push(_dna_query[j-1]);
                        SR2.push(_dna_ref2[j-1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 2 ){
                        SR1.push('-');
                        SQ.push(_dna_query[j-1]);
                        SR2.push(_dna_ref2[j-1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    }else{
                        std::cout << "line 203" << std::endl;
                    }
                }
            } else if ( trackMatrix == -1 ) {
                SQ.push(_dna_query[j-1]);
                SR2.push(_dna_ref2[j-1]);
                SR1.push(_dna_ref1[i-1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i; --j;
            } else if( 1 == trackMatrix ){
                SR1.push('-');
                SQ.push(_dna_query[j-1]);
                SR2.push(_dna_ref2[j-1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            }else if( 2 == trackMatrix ){
                SR1.push('-');
                SQ.push(_dna_query[j-1]);
                SR2.push(_dna_ref2[j-1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            }else if( 3 == trackMatrix ){
                SR1.push(_dna_ref1[i-1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            }else if( 4 == trackMatrix ){
                SR1.push(_dna_ref1[i-1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            }else{
                std::cout << "line 232" << std::endl;
            }
        }
    }
//    std::cout << "line 226" << std::endl;
    delete[] MM;
    //  std::cout << "line 228" << std::endl;
    delete[] TM;
    //std::cout << "line 229" << std::endl;
    return maxScore;
}


int32_t alignment(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                  const int32_t & matchingScore, int32_t mismatchingPenalty,
                  int32_t _open_gap_penalty1, int32_t _extend_gap_penalty1,
                  int32_t _open_gap_penalty2, int32_t _extend_gap_penalty2){
    int32_t length1 = _dna_d.length();
    int32_t length2 = _dna_q.length();
//    std::cout << "line 51" << std::endl;
    int64_t matrixsize = (length1 +1) * (length2+1);
    //  std::cout << "line 53" << std::endl;
    int32_t * MM = new int32_t [matrixsize * 6];
    int8_t * TM = new int8_t [matrixsize * 6];
    //std::cout << "line 56" << std::endl;
    std::fill_n(MM, matrixsize * 6, 0);
    //std::cout << "line 58" << std::endl;
    std::fill_n(TM+matrixsize*0, matrixsize*2, 0); //V and M
    std::fill_n(TM+matrixsize*2, matrixsize, 1); //E1
    std::fill_n(TM+matrixsize*3, matrixsize, 3); //F1
    std::fill_n(TM+matrixsize*4, matrixsize, 2); //E2
    std::fill_n(TM+matrixsize*5, matrixsize, 4); //F2

//    std::cout << "line 65" << std::endl;
    int32_t i=0, j;
    for (j = 0; j < (length2 + 1); ++j  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }
//    std::cout << "line 81" << std::endl;
    j=0;
    for (i = 0; i < (length1 + 1); ++i  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }
//    std::cout << "line 97" << std::endl;
    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;
//    std::cout << "line 106" << std::endl;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            mScore = _dna_d[i-1] == _dna_q[j-1] ? matchingScore : mismatchingPenalty;
            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i-1) * (length2 + 1) + j-1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if( _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            }else{
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if( _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            }else{
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if(_extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            }else{
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if(_extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            }else{
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j], TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }


//    std::cout << "line 401" << std::endl;
    int32_t endPosition1=0;
    int32_t endPosition2=0;
    int32_t maxScore = -1000000;

    i = length1;
    for( j=1; j<=length2; ++j ){
        if( MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore ){ // please do not change >= to >, since we are doing global alignment
            // >= will omit the first similar fragments
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=j;
        }
    }

    j = length2;
    for ( i=1; i<=length1; ++i ) {
        if (MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore) { // please do not change >= to >, since we are doing global alignment
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2 = j;
        }
    }
    i = endPosition1;
    j = endPosition2;


    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];
//    std::cout << "line 154" << std::endl;
    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SD.push(_dna_d[i-1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_q[j-1]);
            SD.push('-');
            --j;
        } else {
            if( 0 == trackMatrix ) { // come from V
                {
                    if ( TM[matrixsize * 0 + i * (length2 + 1) + j] == -1 ) {
                        SQ.push(_dna_q[j-1]);
                        SD.push(_dna_d[i-1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i; --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 3 ) {
                        SD.push(_dna_d[i-1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 4 ) {
                        SD.push(_dna_d[i-1]);
                        SQ.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 1 ) {
                        SD.push('-');
                        SQ.push(_dna_q[j-1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 2 ){
                        SD.push('-');
                        SQ.push(_dna_q[j-1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    }else{
                        std::cout << "line 203" << std::endl;
                    }
                }
            } else if ( trackMatrix == -1 ) {
                SQ.push(_dna_q[j-1]);
                SD.push(_dna_d[i-1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i; --j;
            } else if( 1 == trackMatrix ){
                SD.push('-');
                SQ.push(_dna_q[j-1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            }else if( 2 == trackMatrix ){
                SD.push('-');
                SQ.push(_dna_q[j-1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            }else if( 3 == trackMatrix ){
                SD.push(_dna_d[i-1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            }else if( 4 == trackMatrix ){
                SD.push(_dna_d[i-1]);
                SQ.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            }else{
                std::cout << "line 232" << std::endl;
            }
        }
    }
//    std::cout << "line 226" << std::endl;
    delete[] MM;
    //  std::cout << "line 228" << std::endl;
    delete[] TM;
    //std::cout << "line 229" << std::endl;
    return maxScore;
}










int32_t alignment(const std::string& _dna_ref2, const std::string & _dna_query, const std::string& _dna_ref1, std::stack<char> & SQ, std::stack<char> & SR1, std::stack<char> & SR2, const int32_t & matchingScore,
                  const int32_t & mismatchingPenalty, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2 ){
    int32_t length1 = _dna_ref1.length();
    int32_t length2 = _dna_ref2.length();
//    std::cout << "line 51" << std::endl;
    int64_t matrixsize = (length1 +1) * (length2+1);
    //  std::cout << "line 53" << std::endl;
    int32_t * MM = new int32_t [matrixsize * 6];
    int8_t * TM = new int8_t [matrixsize * 6];
    //std::cout << "line 56" << std::endl;
    std::fill_n(MM, matrixsize * 6, 0);
    //std::cout << "line 58" << std::endl;
    std::fill_n(TM+matrixsize*0, matrixsize*2, 0); //V and M
    std::fill_n(TM+matrixsize*2, matrixsize, 1); //E1
    std::fill_n(TM+matrixsize*3, matrixsize, 3); //F1
    std::fill_n(TM+matrixsize*4, matrixsize, 2); //E2
    std::fill_n(TM+matrixsize*5, matrixsize, 4); //F2

//    std::cout << "line 65" << std::endl;
    int32_t i=0, j;
    for (j = 0; j < (length2 + 1); ++j  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + j*_extend_gap_penalty1, _open_gap_penalty2 + j*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + j*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + j*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 1;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 1;
    }
//    std::cout << "line 81" << std::endl;
    j=0;
    for (i = 0; i < (length1 + 1); ++i  ){
        MM[matrixsize * 0 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 1 + i * (length2 + 1) + j] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + i*_extend_gap_penalty1;
        MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + i*_extend_gap_penalty2;
        TM[matrixsize * 0 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 1 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 2 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 4 + i * (length2 + 1) + j] = 3;
        TM[matrixsize * 5 + i * (length2 + 1) + j] = 3;
    }
//    std::cout << "line 97" << std::endl;
    MM[0] = 0;
    TM[0] = 0; //V
    TM[matrixsize * 1] = 0;//M
    TM[matrixsize * 2] = 1;//E1
    TM[matrixsize * 3] = 3;//F1
    TM[matrixsize * 4] = 2;//E2
    TM[matrixsize * 5] = 4;//F2
    int32_t mScore;
//    std::cout << "line 106" << std::endl;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            mScore = 0;
            mScore += _dna_ref2[j-1] == '-' ? 0 : (_dna_ref1[i-1] == _dna_ref2[j-1] ? matchingScore : mismatchingPenalty);
            mScore += _dna_query[j-1] == '-' ? 0 : (_dna_ref1[i-1] == _dna_query[j-1] ? matchingScore : mismatchingPenalty);
            mScore = mScore /2;

            MM[matrixsize * 1 + i * (length2 + 1) + j] = mScore + MM[matrixsize * 0 + (i-1) * (length2 + 1) + j-1];
            TM[matrixsize * 1 + i * (length2 + 1) + j] = 0;
            if( _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 2 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = 1;
            }else{
                MM[matrixsize * 2 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 2 + i * (length2 + 1) + j] = -1;
            }
            if( _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1] ){
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 4 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = 2;
            }else{
                MM[matrixsize * 4 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + i * (length2 + 1) + j-1];
                TM[matrixsize * 4 + i * (length2 + 1) + j] = -1;
            }

            if(_extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _extend_gap_penalty1 + MM[matrixsize * 3 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = 3;
            }else{
                MM[matrixsize * 3 + i * (length2 + 1) + j] = _open_gap_penalty1 + _extend_gap_penalty1 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 3 + i * (length2 + 1) + j] = -1;
            }
            if(_extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j] >= _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j] ){
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _extend_gap_penalty2 + MM[matrixsize * 5 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = 4;
            }else{
                MM[matrixsize * 5 + i * (length2 + 1) + j] = _open_gap_penalty2 + _extend_gap_penalty2 + MM[matrixsize * 1 + (i -1)* (length2 + 1) + j];
                TM[matrixsize * 5 + i * (length2 + 1) + j] = -1;
            }
            MM[matrixsize * 0 + i * (length2 + 1) + j] = max(MM[matrixsize * 1 + i * (length2 + 1) + j], MM[matrixsize * 2 + i * (length2 + 1) + j], MM[matrixsize * 3 + i * (length2 + 1) + j], MM[matrixsize * 4 + i * (length2 + 1) + j], MM[matrixsize * 5 + i * (length2 + 1) + j], TM[matrixsize * 0 + i * (length2 + 1) + j]);
        }
    }


//    std::cout << "line 401" << std::endl;
    int32_t endPosition1=0;
    int32_t endPosition2=0;
    int32_t maxScore = -1000000;

    i = length1;
    for( j=1; j<=length2; ++j ){
        if( MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore ){ // please do not change >= to >, since we are doing global alignment
            // >= will omit the first similar fragments
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=j;
        }
    }

    j = length2;
    for ( i=1; i<=length1; ++i ) {
        if (MM[matrixsize * 0 + i * (length2 + 1) + j] >= maxScore) { // please do not change >= to >, since we are doing global alignment
            maxScore = MM[matrixsize * 0 + i * (length2 + 1) + j];
            endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2 = j;
        }
    }
    i = endPosition1;
    j = endPosition2;


    int8_t trackMatrix = TM[matrixsize * 0 + endPosition1 * (length2 + 1) + endPosition2];
//    std::cout << "line 154" << std::endl;
    while (i != 0 || j != 0) {
        if (j == 0) {
            SQ.push('-');
            SR2.push('-');
            SR1.push(_dna_ref1[i-1]);
            --i;
        } else if (i == 0) {
            SQ.push(_dna_query[j-1]);
            SR2.push(_dna_ref2[j-1]);
            SR1.push('-');
            --j;
        } else {
            if( 0 == trackMatrix ) { // come from V
                {
                    if ( TM[matrixsize * 0 + i * (length2 + 1) + j] == -1 ) {
                        SQ.push(_dna_query[j-1]);
                        SR2.push(_dna_ref2[j-1]);
                        SR1.push(_dna_ref1[i-1]);
                        trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                        --i; --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 3 ) {
                        SR1.push(_dna_ref1[i-1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 4 ) {
                        SR1.push(_dna_ref1[i-1]);
                        SQ.push('-');
                        SR2.push('-');
                        trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                        --i;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 1 ) {
                        SR1.push('-');
                        SQ.push(_dna_query[j-1]);
                        SR2.push(_dna_ref2[j-1]);
                        trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                        --j;
                    } else if (  TM[matrixsize * 0 + i * (length2 + 1) + j] == 2 ){
                        SR1.push('-');
                        SQ.push(_dna_query[j-1]);
                        SR2.push(_dna_ref2[j-1]);
                        trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                        --j;
                    }else{
                        std::cout << "line 203" << std::endl;
                    }
                }
            } else if ( trackMatrix == -1 ) {
                SQ.push(_dna_query[j-1]);
                SR2.push(_dna_ref2[j-1]);
                SR1.push(_dna_ref1[i-1]);
                trackMatrix = TM[matrixsize * 1 + i * (length2 + 1) + j];
                --i; --j;
            } else if( 1 == trackMatrix ){
                SR1.push('-');
                SQ.push(_dna_query[j-1]);
                SR2.push(_dna_ref2[j-1]);
                trackMatrix = TM[matrixsize * 2 + i * (length2 + 1) + j];
                --j;
            }else if( 2 == trackMatrix ){
                SR1.push('-');
                SQ.push(_dna_query[j-1]);
                SR2.push(_dna_ref2[j-1]);
                trackMatrix = TM[matrixsize * 4 + i * (length2 + 1) + j];
                --j;
            }else if( 3 == trackMatrix ){
                SR1.push(_dna_ref1[i-1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 3 + i * (length2 + 1) + j];
                --i;
            }else if( 4 == trackMatrix ){
                SR1.push(_dna_ref1[i-1]);
                SQ.push('-');
                SR2.push('-');
                trackMatrix = TM[matrixsize * 5 + i * (length2 + 1) + j];
                --i;
            }else{
                std::cout << "line 232" << std::endl;
            }
        }
    }
//    std::cout << "line 226" << std::endl;
    delete[] MM;
    //  std::cout << "line 228" << std::endl;
    delete[] TM;
    //std::cout << "line 229" << std::endl;
    return maxScore;
}






int32_t alignment_position(const std::string& _dna_q, const std::string& _dna_d, std::stack<char> & SQ, std::stack<char> & SD,
                  const int32_t & matchingScore, const int32_t & mismatchingPenalty,
                  const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                  const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2){ // far from done


    int32_t length1 = _dna_d.length(), length2 = _dna_q.length(), i, j, f1, f2, mScore, endPosition1=0, endPosition2=0, maxScore=INT32_MIN;

    int32_t * MP = new int32_t [length1 + 1];
    int32_t * MC = new int32_t [length1 + 1];
    int32_t * E1 = new int32_t [length1 + 1];
    int32_t * E2 = new int32_t [length1 + 1];

    std::fill_n(MP, length1+1, 0);
    std::fill_n(MC, length1+1, 0);
    std::fill_n(E1, length1+1, 0);
    std::fill_n(E2, length1+1, 0);


    for (i = 0; i < (length1 + 1); ++i  ){
        MP[i] = max(_open_gap_penalty1 + i*_extend_gap_penalty1, _open_gap_penalty2 + i*_extend_gap_penalty2);
        E1[i] = _open_gap_penalty1;
        E2[i] = _open_gap_penalty2;
    }

    MP[0] = 0;
    MC[0] = 0;
    for ( j=1; j<=length2; ++j ){
        f1 = _open_gap_penalty1;
        f2 = _open_gap_penalty2;
        for ( i=1; i<=length1; ++i ){

            f1 = (_open_gap_penalty1 + _extend_gap_penalty1 + MP[j-1]) > (_extend_gap_penalty1 + f1) ? (_open_gap_penalty1 + _extend_gap_penalty1 + MP[j-1]) : (_extend_gap_penalty1 + f1);
            f2 = (_open_gap_penalty2 + _extend_gap_penalty2 + MP[j-1]) > (_extend_gap_penalty2 + f2) ? (_open_gap_penalty2 + _extend_gap_penalty2 + MP[j-1]) : (_extend_gap_penalty2 + f2);

            mScore = _dna_d[i-1] == _dna_q[j-1] ? matchingScore : mismatchingPenalty;
            MC[i] = mScore + MP[i-1];

            if( _extend_gap_penalty1 + E1[i] >= _open_gap_penalty1 + _extend_gap_penalty1 + MP[i] ){
                E1[i] = _extend_gap_penalty1 + E1[i];
            }else{
                E1[i] = _open_gap_penalty1 + _extend_gap_penalty1 + MP[i];
            }
            if( _extend_gap_penalty2 + E2[i] >= _open_gap_penalty2 + _extend_gap_penalty2 + MP[i] ){
                E2[i] = _extend_gap_penalty2 + E2[i];
            }else{
                E2[i] = _open_gap_penalty2 + _extend_gap_penalty2 + MP[i];
            }
            MC[i] = max(MC[i], E1[i], f1, E2[i], f2);
        }
        if( MP[length1] < maxScore ){
            endPosition1 = length1;
            endPosition2 = j;
        }
        int32_t * TEMP = MP;
        MP = MC;
        MC = TEMP;
    }
    for ( j=1; j<=length2; ++j ){
        if( MP[length1] < maxScore ){
            endPosition1 = length1;
        }
    }
    return maxScore;
}






int64_t alignSlidingWindow( const std::string& dna_q, const std::string& dna_d, int64_t _length_of_q, int64_t _length_of_d,
                            std::string & _alignment_q, std::string & _alignment_d, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2){
//    std::cout << dna_d << std::endl;
//    std::cout << dna_q << std::endl;
    //2^15 = 32768
    //of the maximum length of the windowSize of is about 32000/2 = 16000
    int32_t databaseStart=1;
    int32_t databaseEnd = 0;
    int32_t queryStart=1;
    int32_t queryEnd = 0;
    int64_t totalScore = 0;
    if( _length_of_d * _length_of_q <= slidingWindowSize * slidingWindowSize){
        std::stack<char> SQ;
        std::stack<char> SD;

        totalScore += needleAlignment(dna_q, dna_d, SQ, SD, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//        std::cout << "line 450" << std::endl;
        while (!SQ.empty()) {
            _alignment_q += SQ.top();
            _alignment_d += SD.top();
            if( SQ.top() != '-' ){
                ++queryStart;
            }
            if( SD.top() != '-' ){
                ++databaseStart;
            }
            SQ.pop();
            SD.pop();
        }
    }else{
//        std::cout << "line 464" << std::endl;
        while( databaseStart<=_length_of_d && queryStart<=_length_of_q){
            databaseEnd=databaseStart+slidingWindowSize;
            queryEnd=queryStart+slidingWindowSize;
            if( databaseEnd>_length_of_d ){
                databaseEnd=_length_of_d;
            }
            if( queryEnd>_length_of_q ){
                queryEnd=_length_of_q;
            }
            std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd );
            std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd );
            std::stack<char> SQ;
            std::stack<char> SD;
            if( slidingWindowSize>1073741824 ){
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            }else{
//                std::cout << "line 580" << std::endl;
//                std::cout << qSeq << "\t" << qSeq << std::endl << "line 505" << std::endl;
                if( extendGapPenalty2+matchingScore < 0){
                    totalScore += alignment(qSeq, dSeq, SQ, SD, matchingScore, mismatchingPenalty+matchingScore, openGapPenalty1+matchingScore, extendGapPenalty1+matchingScore, openGapPenalty2+matchingScore, extendGapPenalty2+matchingScore);
                }else{
                    totalScore += alignment(qSeq, dSeq, SQ, SD, matchingScore, mismatchingPenalty+matchingScore-1, openGapPenalty1+matchingScore-1, extendGapPenalty1+matchingScore-1, openGapPenalty2+matchingScore-1, -1);
                }
            }
            while (!SQ.empty()) {
                _alignment_q += SQ.top();
                _alignment_d += SD.top();
                if( SQ.top() != '-' ){
                    ++queryStart;
                }
                if( SD.top() != '-' ){
                    ++databaseStart;
                }
                SQ.pop();
                SD.pop();
            }
        }
//        std::cout << "line 601" << std::endl;
    }
    int32_t final_indel_length = 0;
    while( databaseStart<=_length_of_d ){
        _alignment_q += '-';
        _alignment_d += dna_d[databaseStart-1];
        ++databaseStart;
        ++final_indel_length;
    }
    while( queryStart<=_length_of_q ){
        _alignment_q += dna_q[queryStart-1];
        _alignment_d += '-';
        ++queryStart;
        ++final_indel_length;
    }
    if( final_indel_length > 0 ){
        totalScore += max(openGapPenalty1 + extendGapPenalty1*final_indel_length, openGapPenalty2+extendGapPenalty2*final_indel_length);
    }
//    std::cout << "line 618" << std::endl;
    return totalScore;
}





int64_t alignSlidingWindow( std::string& align_ref2, std::string& align_query, const std::string& dna_ref1,
                            std::string & align_ref1, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2){

    int32_t databaseStart=1;
    int32_t databaseEnd = 0;
    int32_t queryStart=1;
    int32_t queryEnd = 0;
    int64_t totalScore = 0;
    int64_t _length_of_q = align_ref2.size();
    int64_t _length_of_d = dna_ref1.size();
    std::string dna_ref2 = align_ref2;
    std::string dna_query = align_query;
    align_ref2 = "";
    align_query = "";

    if( _length_of_d * _length_of_q <= slidingWindowSize * slidingWindowSize){
        std::stack<char> SQ;
        std::stack<char> SR1;
        std::stack<char> SR2;
//        std::cout << "line 1099" << std::endl;
        totalScore += needleAlignment(dna_ref2, dna_query, dna_ref1, SQ, SR1, SR2, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        while (!SQ.empty()) {
            align_query += SQ.top();
            align_ref1 += SR1.top();
            align_ref2 += SR2.top();
            if( SQ.top() != '-' || SR2.top() != '-' ){
                ++queryStart;
            }
            if( SR1.top() != '-'  ){
                ++databaseStart;
            }
            SQ.pop();
            SR2.pop();
            SR1.pop();
        }
    }else{
//        std::cout << "line 1116" << std::endl;
        while( databaseStart<=_length_of_d && queryStart<=_length_of_q){
            databaseEnd=databaseStart+slidingWindowSize;
            queryEnd=queryStart+slidingWindowSize;
            if( databaseEnd>_length_of_d ){
                databaseEnd=_length_of_d;
            }
            if( queryEnd>_length_of_q ){
                queryEnd=_length_of_q;
            }
//            std::cout << dna_query.size() << "\t" << queryStart << "\t" << queryEnd << std::endl;
//            std::cout << dna_ref2.size() << "\t" << queryStart << "\t" << queryEnd << std::endl;
//            std::cout << dna_ref1.size() << "\t" << databaseStart << "\t" << databaseEnd << std::endl;
            std::string qSeq = getSubsequence(dna_query, queryStart, queryEnd );
            std::string d2Seq = getSubsequence(dna_ref2, queryStart, queryEnd );
            std::string d1Seq = getSubsequence(dna_ref1, databaseStart, databaseEnd );
//            std::cout << "line 1129" << std::endl;
            std::stack<char> SQ;
            std::stack<char> SR1;
            std::stack<char> SR2;
            if( slidingWindowSize>1073741824 ){
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            }else{
                if( extendGapPenalty2+matchingScore < 0){
                    totalScore += alignment(d2Seq, qSeq, d1Seq, SQ, SR1, SR2, matchingScore, mismatchingPenalty+matchingScore, openGapPenalty1+matchingScore, extendGapPenalty1+matchingScore, openGapPenalty2+matchingScore, extendGapPenalty2+matchingScore);
                }else{
                    totalScore += alignment(d2Seq, qSeq, d1Seq, SQ, SR1, SR2, matchingScore, mismatchingPenalty+matchingScore-1, openGapPenalty1+matchingScore-1, extendGapPenalty1+matchingScore-1, openGapPenalty2+matchingScore-1, -1);
                }
            }
            while (!SQ.empty()) {
                align_query += SQ.top();
                align_ref1 += SR1.top();
                align_ref2 += SR2.top();
                if( SQ.top() != '-' || SR2.top() != '-' ){
                    ++queryStart;
                }
                if( SR1.top() != '-'  ){
                    ++databaseStart;
                }
                SQ.pop();
                SR2.pop();
                SR1.pop();
            }
        }
    }
    int32_t final_indel_length = 0;
    while( databaseStart<=_length_of_d ){
        align_query += '-';
        align_ref2 += '-';
        align_ref1 += dna_ref1[databaseStart-1];
        ++databaseStart;
        ++final_indel_length;
    }
    while( queryStart<=_length_of_q ){
        align_query += dna_query[queryStart-1];
        align_ref2 += dna_ref2[queryStart-1];
        align_ref1 += '-';
        ++queryStart;
        ++final_indel_length;
    }
    if( final_indel_length > 0 ){
        totalScore += max(openGapPenalty1 + extendGapPenalty1*final_indel_length, openGapPenalty2+extendGapPenalty2*final_indel_length);
    }
    return totalScore;
}

/*

int64_t alignSlidingWindow( std::vector<std::string> & align_refs, std::vector<std::string> & align_queries, const std::vector<std::string> & dna_refs,
                            const std::vector<std::string> & dna_queries, const int64_t & slidingWindowSize,
                            const int32_t & matchingScore, const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2){

    int32_t databaseStart = 1;
    int32_t databaseEnd = 0;
    int32_t queryStart = 1;
    int32_t queryEnd = 0;
    int64_t totalScore = 0;
    int64_t _length_of_q = dna_queries[0].size();
    int64_t _length_of_d = dna_refs[0].size();
    int32_t n = dna_refs.size();
    int32_t m = dna_queries.size();

    align_refs.clear();
    align_queries.clear();
    for( int32_t i = 0; i < n; ++i ){
        align_refs.push_back("");
    }
    for( int32_t i = 0; i < m; ++i ){
        align_queries.push_back("");
    }

    if( _length_of_d * _length_of_q <= slidingWindowSize * slidingWindowSize){
        std::vector<std::stack<char>> SR (n);
        std::vector<std::stack<char>> SQ (m);

        totalScore += needleAlignment(dna_ref2, dna_query, dna_ref1, SQ, SR1, SR2, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
        while (!SQ.empty()) {
            bool sri = false;
            bool sqi = false;
            for( int32_t i = 0; i < n; ++i ){
                align_refs[i] += SR[i].top();
                if( SR[i].top() != '-' ){
                    sri = true;
                }
            }
            for( int32_t i = 0; i < m; ++i ){
                align_queries[i] += SQ[i].top();
                if( SQ[i].top() != '-' ){
                    sqi = true;
                }
            }
            if( sri ){
                ++databaseStart;
            }
            if( sqi ){
                ++queryStart;
            }

            for( int32_t i = 0; i < n; ++i ){
                SR[i].pop();
            }
            for( int32_t i = 0; i < m; ++i ){
                SQ[i].pop();
            }
        }
    }else{
        std::cout << "line 1116" << std::endl;
        while( databaseStart<=_length_of_d && queryStart<=_length_of_q){
            databaseEnd=databaseStart+slidingWindowSize;
            queryEnd=queryStart+slidingWindowSize;
            if( databaseEnd>_length_of_d ){
                databaseEnd=_length_of_d;
            }
            if( queryEnd>_length_of_q ){
                queryEnd=_length_of_q;
            }
            std::string qSeq =

            std::string qSeq = getSubsequence(dna_query, queryStart, queryEnd );
            std::string d2Seq = getSubsequence(dna_ref2, queryStart, queryEnd );
            std::string d1Seq = getSubsequence(dna_ref1, databaseStart, databaseEnd );
            std::cout << "line 1129" << std::endl;
            std::stack<char> SQ;
            std::stack<char> SR1;
            std::stack<char> SR2;
            if( slidingWindowSize>1073741824 ){
                std::cout << "the windows size is too large" << std::endl;
                exit(1);
            }else{
                if( extendGapPenalty2+matchingScore < 0){
                    totalScore += alignment(d2Seq, qSeq, d1Seq, SQ, SR1, SR2, matchingScore, mismatchingPenalty+matchingScore, openGapPenalty1+matchingScore, extendGapPenalty1+matchingScore, openGapPenalty2+matchingScore, extendGapPenalty2+matchingScore);
                }else{
                    totalScore += alignment(d2Seq, qSeq, d1Seq, SQ, SR1, SR2, matchingScore, mismatchingPenalty+matchingScore-1, openGapPenalty1+matchingScore-1, extendGapPenalty1+matchingScore-1, openGapPenalty2+matchingScore-1, -1);
                }
            }
            while (!SQ.empty()) {
                align_query += SQ.top();
                align_ref1 += SR1.top();
                align_ref2 += SR2.top();
                if( SQ.top() != '-' || SR2.top() != '-' ){
                    ++queryStart;
                }
                if( SR1.top() != '-'  ){
                    ++databaseStart;
                }
                SQ.pop();
                SR2.pop();
                SR1.pop();
            }
        }
    }
    int32_t final_indel_length = 0;
    while( databaseStart<=_length_of_d ){
        align_query += '-';
        align_ref2 += '-';
        align_ref1 += dna_ref1[databaseStart-1];
        ++databaseStart;
        ++final_indel_length;
    }
    while( queryStart<=_length_of_q ){
        align_query += dna_query[queryStart-1];
        align_ref2 += dna_ref2[queryStart-1];
        align_ref1 += '-';
        ++queryStart;
        ++final_indel_length;
    }
    if( final_indel_length > 0 ){
        totalScore += max(openGapPenalty1 + extendGapPenalty1*final_indel_length, openGapPenalty2+extendGapPenalty2*final_indel_length);
    }
    return totalScore;
}
*/



int64_t alignSlidingWindow_minimap2( const std::string& dna_q, const std::string& dna_d, int64_t _length_of_q, int64_t _length_of_d,
                            std::string & _alignment_q, std::string & _alignment_d, const int64_t & slidingWindowSize,
                            const int32_t & mismatchingPenalty, const int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2){
    int8_t a = 0, b = mismatchingPenalty < 0? mismatchingPenalty : -mismatchingPenalty; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = dna_d.length(), ql = dna_q.length();
    uint8_t *ts, *qs, c[256];

    const char *tseq = dna_d.c_str();
    const char *qseq = dna_q.c_str();

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
//    std::cout << "line 369" << std::endl;
//    void *km = 0;
#ifdef __AVX512BW__
//    std::cout << "using AVX512" << std::endl;
    ksw_extd2_avx512(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, & ez);
#elif __AVX2__
//    std::cout << "using AVX2" << std::endl;
    ksw_extd2_avx2(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, & ez);
#else
//    std::cout << "using SSE" << std::endl;
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, -openGapPenalty1, -extendGapPenalty1, -openGapPenalty2, -extendGapPenalty2, slidingWindowSize, -1, 0, 0, & ez);
#endif
  //  std::cout << "line 371" << std::endl;
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
                _alignment_q += dna_q[text_pos];
                _alignment_d += dna_d[pattern_pos];
                pattern_pos++;
                text_pos++;
            }
        } else if (cLetter == 'X') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += dna_q[text_pos];
                _alignment_d += dna_d[pattern_pos];
                pattern_pos++;
                text_pos++;
            }
        } else if (cLetter == 'I') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += dna_q[text_pos];
                _alignment_d += '-';
                text_pos++;
            }
        } else if (cLetter == 'D') {
            for (int j = 1; j <= cLen; ++j) {
                _alignment_q += '-';
                _alignment_d += dna_d[pattern_pos];
                pattern_pos++;
            }
        }
    }
    int32_t totalScore = ez.score;
    free(ez.cigar); free(ts); free(qs);
//    km_destroy(km);
    return totalScore;
}

//
//const std::string WHITESPACE = " \n\r\t\f\v";
//std::string rtrim(const std::string& s){
//    size_t end = s.find_last_not_of(WHITESPACE);
//    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
//}

std::string GetStdoutFromCommand(std::string & cmd) {

    std::string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1");

    stream = popen(cmd.c_str(), "r");

    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    rtrim(data);
    return data;
}



int64_t alignSlidingWindow_minimap2_or_NW(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                                       affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                                       const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                                       const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m ){

//    std::cout << "line 523" << std::endl;
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int64_t _length_of_q=dna_q.size();
    int64_t _length_of_d=dna_d.size();

    //check all Ns end
//    std::cout << dna_d << std::endl << dna_q << std::endl << "line 659" << std::endl;
    int32_t longerSeqLength = max(_length_of_d, _length_of_q);
    if( _length_of_d*_length_of_q <= (slidingWindowSize*slidingWindowSize*30) ){ // this calculated via RAM cost
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = -1;
//        std::cout << "line 789" << std::endl;
        totalScore =  alignSlidingWindow_minimap2( dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                   mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//        std::cout << "minimap2:" << std::to_string(_length_of_d) << "\t" << std::to_string(_length_of_q) << std::endl;
    }else if( longerSeqLength*slidingWindowSize*2 <= (slidingWindowSize*slidingWindowSize*30) ){ // this calculated with RAM cost
        /*the above parameter settings were based on RAM cost*/
        int32_t adjustBandWidth = (slidingWindowSize*slidingWindowSize*30)/2/longerSeqLength;
    //    std::cout << "line 878 adjustBandWidth:" << std::to_string(adjustBandWidth) << std::endl;

        totalScore =  alignSlidingWindow_minimap2( dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, adjustBandWidth,
                                                   mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//        std::cout << "minimap2_band:" << std::to_string(_length_of_d) << "\t" << std::to_string(_length_of_q) << std::endl;
    }else{
      //  std::cout << "line 882" << std::endl;
        totalScore =  alignSlidingWindow( dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                          mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//        std::cout << "slidingwindow:" << std::to_string(_length_of_d) << "\t" << std::to_string(_length_of_q) << std::endl;
    }
    return totalScore;
}



int64_t alignSlidingWindow_local_wfa_v2(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                                       affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                                       const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                                       const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m ){
//    std::cout << "line 952 alignSlidingWindow_local_wfa_v2" << std::endl;
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
//    std::cout << dna_d << std::endl << dna_q << std::endl;
    int64_t _length_of_q=dna_q.size();
    int64_t _length_of_d=dna_d.size();

    //check all Ns begin
    std::string temp2 = dna_q;
    temp2.erase(std::remove(temp2.begin(), temp2.end(), 'N'), temp2.end());

    std::string temp1 = dna_d;
    temp1.erase(std::remove(temp1.begin(), temp1.end(), 'N'), temp1.end());

    if( temp1.size() ==0 || temp2.size()==0 ){
        _alignment_q =dna_q;
        _alignment_d =dna_d;
        while( _alignment_q.size() < _alignment_d.size() ){
            _alignment_q =  _alignment_q + "-";
        }
        while( _alignment_d.size() < _alignment_q.size() ){
            _alignment_d =  _alignment_d + "-";
        }
        return totalScore;
    }

    //check all Ns end
//    std::cout << dna_d << std::endl << dna_q << std::endl << "line 659" << std::endl;
    int32_t longerSeqLength = max(_length_of_d, _length_of_q);
  //  std::cout << "line 830" << std::endl;
    std::string qSeq = dna_q;
    std::string dSeq = dna_d;
    const char *pattern = dna_d.c_str();
    const char *text = dna_q.c_str();
    //std::cout << "line 987" << std::endl;
    // Allocate
    affine2p_wavefronts_t* const affine2p_wavefronts = affine2p_wavefronts_new_complete( strlen(pattern), strlen(text), penalties, mm_allocator);
    //std::cout << "line 990" << std::endl;

    // for sliding windows NW, every cell costs 30 bytes memory, so for a pair sequences alignment of 30kb, it would cost 27 Gb memory
    // for WFA, each anchorwave cost 27 byte, so we could have 30000*30000*30/21 = 1285714286 anchorwaves at maximum

    int64_t anchorwavenumberthreshold = (slidingWindowSize*slidingWindowSize*30)/21;
//    std::cout << "anchorwavenumberthreshold:" << std::to_string(anchorwavenumberthreshold) << std::endl;
// Align
    int64_t wfa_code = affine2p_wavefronts_align_anchorwavenumberthreshold_version(affine2p_wavefronts, pattern, strlen(pattern), text, strlen(text), anchorwavenumberthreshold);
  //  std::cout << "line 996" << std::endl;
//    int64_t wfa_code=0;
//    affine2p_wavefronts_align(affine2p_wavefronts, pattern,strlen(pattern),text,strlen(text));
    if(0==wfa_code){
//        std::cout << "WFA:" << std::to_string(_length_of_d) << "\t" << std::to_string(_length_of_q) << std::endl;
        cigar_t *const edit_cigar = &(affine2p_wavefronts->cigar);
        char *const operations = edit_cigar->operations;

        int i, pattern_pos = 0, text_pos = 0;
        char lastCigar = '-';
        for (i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i) {
            if (operations[i] == 'M') {
                _alignment_q += qSeq[text_pos];
                _alignment_d += dSeq[pattern_pos];
                pattern_pos++;
                text_pos++;
            }else if (operations[i] == 'X') {
                _alignment_q += qSeq[text_pos];
                _alignment_d += dSeq[pattern_pos];
                pattern_pos++;
                text_pos++;
            }else if (operations[i] == 'I') {
                _alignment_q += qSeq[text_pos];
                _alignment_d += '-';
                text_pos++;
                if( operations[i] != lastCigar ){
                    totalScore += openGapPenalty1;
                }
            }  else if (operations[i] == 'D') {
                _alignment_q += '-';
                _alignment_d += dSeq[pattern_pos];
                pattern_pos++;
                if( operations[i] != lastCigar ){
                    totalScore += openGapPenalty1;
                }
            }
            lastCigar = operations[i];
        }
        totalScore = cigar_score_gap_affine2p(edit_cigar, penalties);
    }else{
//        std::cout << "WFA failed, anchorwave count:" << std::to_string(wfa_code) << std::endl;
//        std::cout << "WFA_failed:" << std::to_string(_length_of_d) << "\t" << std::to_string(_length_of_q) << std::endl;
        totalScore = alignSlidingWindow_minimap2_or_NW(  dna_q,  dna_d, _alignment_q, _alignment_d,
                                                         penalties,  mm_allocator, slidingWindowSize, wfaSize, matchingScore,
                                                         mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2,
                                                         extendGapPenalty2, min_wavefront_length, max_distance_threshold, m );
    }
    affine2p_wavefronts_delete(affine2p_wavefronts);
    return totalScore;
}

int64_t alignSlidingWindow(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m, std::map<std::string, std::string>& parameters ){

//    return alignSlidingWindow_call_sys( dna_q, dna_d, _alignment_q, _alignment_d, penalties, mm_allocator, slidingWindowSize, wfaSize, matchingScore,
//                                   mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
//                                   min_wavefront_length, max_distance_threshold, m, parameters);
//    return alignSlidingWindow_local_wfa( dna_q, dna_d, _alignment_q, _alignment_d, penalties, mm_allocator, slidingWindowSize, wfaSize, matchingScore,
//                                         mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
//                                         min_wavefront_length, max_distance_threshold, m);
    return alignSlidingWindow_local_wfa_v2( dna_q, dna_d, _alignment_q, _alignment_d, penalties, mm_allocator, slidingWindowSize, wfaSize, matchingScore,
                                         mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                                         min_wavefront_length, max_distance_threshold, m);
}


int64_t alignSlidingWindowNW(  std::string& dna_q,  std::string& dna_d, std::string & _alignment_q, std::string & _alignment_d,
                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1,  const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m ){
//    std::cout << dna_d << std::endl << dna_q << std::endl << "line 722" << std::endl;
    int64_t totalScore = 0;
    _alignment_q = "";
    _alignment_d = "";
    int64_t _length_of_q=dna_q.size();
    int64_t _length_of_d=dna_d.size();

    //check all Ns begin
    std::string temp2 = dna_q;
    temp2.erase(std::remove(temp2.begin(), temp2.end(), 'N'), temp2.end());

    std::string temp1 = dna_d;
    temp1.erase(std::remove(temp1.begin(), temp1.end(), 'N'), temp1.end());

    if( temp1.size() ==0 || temp2.size()==0 ){
        _alignment_q =dna_q;
        _alignment_d =dna_d;
        int32_t _alignment_q_length = _alignment_q.length();
        int32_t _alignment_d_length = _alignment_d.length();
        if( _alignment_q_length < _alignment_d_length ){
            int32_t dis = _alignment_d_length - _alignment_q_length;
            _alignment_q =  _alignment_q + std::string(dis, '-');
        }else if( _alignment_d_length < _alignment_q_length ){
            int32_t dis = _alignment_q_length - _alignment_d_length;
            _alignment_d =  _alignment_d + std::string(dis, '-');
        }
        return totalScore;
    }
//    std::cout << "line 722" << std::endl;
    totalScore =  alignSlidingWindow( dna_q, dna_d, _length_of_q, _length_of_d, _alignment_q, _alignment_d, slidingWindowSize, matchingScore,
                                          mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//    std::cout << "realigned_sliding_window:" << std::to_string(_length_of_d) << "\t" << std::to_string(_length_of_q) << std::endl;
    return totalScore;
}



int64_t alignSlidingWindow(  std::string& dna_ref1,  std::string& dna_ref2, std::string& dna_query, std::string & _alignment_ref1, std::string & _alignment_ref2, std::string & _alignment_query,
                             affine2p_penalties_t* const penalties, mm_allocator_t* const mm_allocator, const int64_t & slidingWindowSize, const int32_t & wfaSize, const int32_t & matchingScore,
                             const int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                             const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, const Scorei & m, std::map<std::string, std::string>& parameters ){

    return alignSlidingWindow_local_wfa_v2( dna_query, dna_ref2, _alignment_query, _alignment_ref2, penalties, mm_allocator, slidingWindowSize, wfaSize, matchingScore,
                                            mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                                            min_wavefront_length, max_distance_threshold, m);
}
