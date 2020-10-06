/*
 * =====================================================================================
 *
 *       Filename:  sequenceAlignment.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/25/2019 01:13:39
 *       Revision:  08/12/2019
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 *         Lynn Johnson performed review reading on 08/12/2019
 *
 * =====================================================================================
 */


/*************************************************************************
 *
 *
 * Implementing all the sequence alignment dynamic programming algorithms here
 *
 *
 *
 *
 ************************************************************************/

#include "sequenceAlignment.h"


static const int32_t SCORE_OUT_BANDED_ALIGNMENT_REGION = -536870912;
static const int32_t MAX_ALIGNMENT_SCORE_USING_INT8 = 125;
static const int32_t MAX_ALIGNMENT_SCORE_USING_INT16 = 32767;

/**
 * this is a standarded smith waterman implementation
 * should be used as benchmark to test other implementations
 * reverseAlignment is used to align the reverse sequence for CNS detection
 * if it is true, the algorithm try hardly to start the alignment from the first base pair for both query and target sequence
 * */

std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                   bool reverseAlignment, bool returnCigar){

    int32_t ** M = new int32_t *[length1 + 1];
    int32_t ** E = new int32_t *[length1 + 1];
    int32_t ** F = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
        M[i] = new int32_t [length2 + 1];
        E[i] = new int32_t [length2 + 1];
        F[i] = new int32_t [length2 + 1];
        std::fill_n(M[i], length2+1, 0);
        std::fill_n(E[i], length2+1, 0);
        std::fill_n(F[i], length2+1, 0);
    }

    endPosition1=0;
    endPosition2=0;
    maxScore = 0;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            E[i][j] = (_open_gap_penalty + M[i][j-1]) > (_extend_gap_penalty + E[i][j-1]) ? (_open_gap_penalty + M[i][j-1]) : (_extend_gap_penalty + E[i][j-1]);
//            E[i][j] = E[i][j] > (_open_gap_penalty + F[i][j-1]) ? E[i][j] : (_open_gap_penalty + F[i][j-1]);

            F[i][j] = (_open_gap_penalty + M[i-1][j]) > (_extend_gap_penalty + F[i-1][j]) ?( _open_gap_penalty + M[i-1][j]) : (_extend_gap_penalty + F[i-1][j]);
//            F[i][j] = F[i][j] > (_open_gap_penalty + E[i-1][j]) ? F[i][j] : (_open_gap_penalty + E[i-1][j]);

            M[i][j] = (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) > E[i][j] ? (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) : E[i][j];
            M[i][j] = M[i][j] > F[i][j] ? M[i][j] : F[i][j];
            M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
            F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
            E[i][j] = E[i][j] > 0 ? E[i][j] : 0;

            if( M[i][j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                // >= will omit the first similar fragments
                maxScore = M[i][j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
    }
    //std::cout << "maxScore " << maxScore << std::endl;
    std::vector<uint32_t> cigar;
    if(returnCigar ){ // cigar is only needed for reverse alignment
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int32_t ii = endPosition1;
        int32_t jj = endPosition2;
        while (ii>0 && jj>0) {
            if (M[ii][jj] == M[ii - 1][jj - 1] + m.getScore(seq1[ii - 1], seq2[jj - 1]) ) {
                --ii;
                --jj;
                op = 0;
            } else if (M[ii][jj] == E[ii][jj]) {
                --jj;
                op = 1;
            } else if (M[ii][jj] == F[ii][jj]) {
                --ii;
                op = 2;
            } else if(M[ii][jj]==0){
                break;
            } else {
                std::cout << "M[ii][jj] " << M[ii][jj] << " E[ii][jj] " << E[ii][jj] << " F[ii][jj] " << F[ii][jj] << std::endl;
                std::cout << "there is something wrong with smith-waterman algorithm in line 114" << std::endl;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        if (!reverseAlignment){
            std::reverse(cigar.begin(),cigar.end());
        }
        // trace back end
    }

    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E[i];
        delete[] F[i];
    }
    delete[] M;
    delete[] E;
    delete[] F;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}


// a smithwaterman method which could return the maximum score and the position where the maximum score located
void SmithWaterman_non_avx(int8_t *ref, int8_t *read, const int32_t &refLen,
                   const int32_t &readLen, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & scorei, bool & positions){
    int32_t length1 = refLen;
    int32_t length2 = readLen;
    int8_t * seq1 = ref;
    int8_t * seq2 = read;

    int32_t * M1 = new int32_t [length2 + 1];
    int32_t * M2 = new int32_t [length2 + 1];
    int32_t * F = new int32_t [length2 + 1];
    std::fill_n(M1, length2+1, 0);
    std::fill_n(F, length2+1, 0);
    std::fill_n(M2, length2+1, 0);

    int32_t e;
    int32_t * t;

    int32_t i, j;
    maxScore = 0;
    int32_t mscore;
    if (positions){
        endPosition1=0;
        endPosition2=0;
        for ( i=1; i<=length1; ++i ){
            e=0;
            for ( j=1; j<=length2; ++j ){
                e = (_open_gap_penalty + M2[j-1]) > (_extend_gap_penalty + e) ? (_open_gap_penalty + M2[j-1]) : (_extend_gap_penalty + e);

                F[j] = (_open_gap_penalty + M1[j]) >  (_extend_gap_penalty + F[j]) ?(_open_gap_penalty + M1[j]) :  (_extend_gap_penalty + F[j]);

                mscore = (scorei.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) > e ? (scorei.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) : e;
                mscore = mscore > F[j] ? mscore : F[j];

                M2[j] = mscore > 0 ? mscore : 0;
                F[j] = F[j] > 0 ? F[j] : 0;
                e = e > 0 ? e : 0;
                if( mscore > maxScore ){ // please do not change > to >=, since we are doing local alignment, mscore is no need to change to M2[J]
                    // >= will omit the first similar fragments
                    maxScore = mscore;
                    endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                    endPosition2=j;
                }
            }
            t = M1;
            M1 = M2;
            M2 = t;
        }
        --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
        --endPosition2;
    }else{
        for ( i=1; i<=length1; ++i ){
            e=0;
            for ( j=1; j<=length2; ++j ){
                e = _open_gap_penalty + M2[j-1] > _extend_gap_penalty + e ? _open_gap_penalty + M2[j-1] : _extend_gap_penalty + e;
                F[j] = (_open_gap_penalty + M1[j]) >  (_extend_gap_penalty + F[j]) ?(_open_gap_penalty + M1[j]) :  (_extend_gap_penalty + F[j]);

                mscore = (scorei.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) > e ? (scorei.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) : e;
                mscore = mscore > F[j] ? mscore : F[j];

                M2[j] = mscore > 0 ? mscore : 0; // mscore is no need to change to M2[j]
                F[j] = F[j] > 0 ? F[j] : 0;
                e = e > 0 ? e : 0;
                maxScore = maxScore > mscore ? maxScore : mscore;
            }
            t = M1;
            M1 = M2;
            M2 = t;
        }
    }
    delete[] M1;
    delete[] F;
    delete[] M2;
}

// a smithwaterman method which could return the maximum score and the position where the maximum score located
// this is a SIMD speedup version
void SmithWaterman(int8_t *ref, int8_t *read, const int32_t &refLen,
                   const int32_t &readLen, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & scorei, bool & positions) {
#ifdef __AVX2__
    if ( ( (-refLen*_extend_gap_penalty - _open_gap_penalty)< MAX_ALIGNMENT_SCORE_USING_INT8 || (-readLen*_extend_gap_penalty - _open_gap_penalty) < MAX_ALIGNMENT_SCORE_USING_INT8 )&&
            ( refLen*scorei.getScore(0,0) < MAX_ALIGNMENT_SCORE_USING_INT8 || readLen*scorei.getScore(0,0) < MAX_ALIGNMENT_SCORE_USING_INT8 ))  {
        positions = true;
        maxScore = 0;

        int8_t bias = 0;
        size_t _segLen = (readLen + 31) / 32, _n = 5;
        __m256i *profile = (__m256i *) malloc(_n * _segLen * sizeof(__m256i));
        int32_t ai, i, j, k, n = 5;
        for (ai = 0; ai < n; ++ai) {
            int8_t *s = (int8_t * )(profile + ai * _segLen);
            for (i = 0; i < _segLen; ++i) {
                j = i;
                for (k = 0; k < 32; ++k) {
                    *s++ = j >= readLen ? bias : scorei.getScore(ai, read[j]);
                    j += _segLen;
                }
            }
        }

        int8_t *_similarity_matrix = new int8_t[readLen + 1];
        _similarity_matrix[0] = 0;

        const __m256i vZero = _mm256_set1_epi8(0);

        __m256i *pvHStore = (__m256i *) calloc(_segLen, sizeof(__m256i));
        __m256i *pvHLoad = (__m256i *) calloc(_segLen, sizeof(__m256i));
        __m256i *pvE = (__m256i *) calloc(_segLen, sizeof(__m256i));

        __m256i *pv;
        __m256i *vP;
        __m256i vGapO, vGapE;
        // INDEL begin vector
        vGapO = _mm256_set1_epi8((int8_t)_open_gap_penalty);
        // INDEL extension vector
        vGapE = _mm256_set1_epi8((int8_t)_extend_gap_penalty);

        int16_t z;

        //for outer loop
        __m256i vH, e;
        int8_t deletionScore;
        int8_t *a;
        for (i = 0; i < refLen; ++i) {//important

            vP = profile + ref[i] * _segLen;

            vH = _mm256_loadu_si256(pvHStore + (_segLen - 1));// the last one
            vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);

            pv = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
            // inner loop to process the query sequence
            for (j = 0; j < _segLen; ++j) {
                vH = _mm256_add_epi8(vH, _mm256_loadu_si256(vP + j)); // add the scoring profile to vH
                vH = _mm256_max_epi8(vH, vZero); /* vH will be always > 0 */

                /* Get max from vH, vE and vF. */
                e = _mm256_loadu_si256(pvE + j);
                vH = _mm256_max_epi8(vH, e);

                /* Save vH values. */
                _mm256_storeu_si256(pvHStore + j, vH);

                /* Update vE value. */
                vH = _mm256_add_epi8(vH, vGapO); /* saturation arithmetic, result >= 0 */
                e = _mm256_add_epi8(e, vGapE);
                e = _mm256_max_epi8(e, vH);
                e = _mm256_max_epi8(e, vZero);
                _mm256_storeu_si256(pvE + j, e);

                /* Load the next vH. */
                vH = _mm256_loadu_si256(pvHLoad + j);

            }
            // put the score into the score matrix begin
            for (j = 0; j < _segLen; ++j) {
                a = (int8_t * )(pvHStore + j);
                for (z = 0; z < 32; ++z) {
                    if ((j + z * _segLen) < readLen) {
                        _similarity_matrix[j + z * _segLen + 1] = a[z];
                    }
                }
            }
            deletionScore = 0;
            for (j = 0; j < readLen; ++j) {
                deletionScore = _similarity_matrix[j] + (int8_t)(_open_gap_penalty) > deletionScore + _extend_gap_penalty ?
                                _similarity_matrix[j] + (int8_t)(_open_gap_penalty) : deletionScore + _extend_gap_penalty;
                if (deletionScore > _similarity_matrix[j + 1]) {
                    _similarity_matrix[j + 1] = deletionScore;
                }
                deletionScore = deletionScore > 0 ? deletionScore : 0;
                if (_similarity_matrix[j+1] > maxScore) {
                    endPosition1 = i;
                    endPosition2 = j;
                    maxScore = _similarity_matrix[j+1];
                }
            }

            //updata vH begin
            for (j = 0; j < _segLen; j++) {
                a = (int8_t * )(pvHStore + j);
                for (z = 0; z < 32; ++z) {
                    if ((j + z * _segLen) < readLen) {
                        a[z] = _similarity_matrix[j + z * _segLen + 1];
                    }
                }
            }//updata vH end
        }
        free(pvHStore);
        free(pvHLoad);
        free(profile);
        free(pvE);
    } else if ( ( (-refLen*_extend_gap_penalty - _open_gap_penalty)<MAX_ALIGNMENT_SCORE_USING_INT16 || (-readLen*_extend_gap_penalty - _open_gap_penalty)<MAX_ALIGNMENT_SCORE_USING_INT16 )&&
          ( refLen*scorei.getScore(0,0) < MAX_ALIGNMENT_SCORE_USING_INT16 || readLen*scorei.getScore(0,0) < MAX_ALIGNMENT_SCORE_USING_INT16 ))  {
        positions = true;
        maxScore = 0;
//        std::cout << "simd line 242" << std::endl;
        int8_t bias = 0;
        size_t _segLen = ( readLen+ 15)/16, _n = 5;
        __m256i* profile = (__m256i*)malloc(_n * _segLen * sizeof(__m256i));
        int32_t ai, i, j, k, n=5;
        for( ai =0; ai<n; ++ai ){
            int16_t* s = (int16_t*)(profile+ai*_segLen);
            for( i=0; i<_segLen; ++i ){
                j = i;
                for( k=0; k<16; ++k ){
                    *s++ = j>= readLen ? bias : scorei.getScore(ai, read[j]);
                    j += _segLen;
                }
            }
        }

        int16_t * _similarity_matrix = new int16_t[readLen + 1];
        _similarity_matrix[0] = 0;
        const __m256i vZero = _mm256_set1_epi16(0);

        __m256i* pvHStore = (__m256i*) calloc(_segLen, sizeof(__m256i));
        __m256i* pvHLoad = (__m256i*) calloc(_segLen, sizeof(__m256i));
        __m256i* pvE = (__m256i*) calloc(_segLen, sizeof(__m256i));

        __m256i* pv;
        __m256i* vP;
        __m256i vGapO, vGapE;
        // INDEL begin vector
        vGapO = _mm256_set1_epi16((int16_t)_open_gap_penalty); // here I use epi16, because for the score matrix I used int16_t
        // INDEL extension vector
        vGapE = _mm256_set1_epi16((int16_t)_extend_gap_penalty);

        int16_t z;

        //for outer loop
        __m256i vH, e;
        int16_t deletionScore;
        int16_t* a;
        for (i = 0; i<refLen; ++i) {//important

//        __m256i vF = vZero, vMaxColumn = vZero;

            vP = profile + ref[i] * _segLen;

            vH = _mm256_loadu_si256(pvHStore+(_segLen-1));// the last one
            vH = _mm256_alignr_epi8(vH, _mm256_permute2x128_si256(vH, vH, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);

            pv = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
            // inner loop to process the query sequence
            for (j = 0; j < _segLen; ++j) {
                vH = _mm256_add_epi16(vH, _mm256_loadu_si256(vP+j)); // add the scoring profile to vH
                vH = _mm256_max_epi16(vH, vZero); /* vH will be always > 0 */


                /* Get max from vH, vE and vF. */
                e = _mm256_loadu_si256(pvE + j);
                vH = _mm256_max_epi16(vH, e);
//            vH = _mm256_max_epi16(vH, vF);
//            vMaxColumn = _mm256_add_epi16(vMaxColumn, vH);

                /* Save vH values. */
                _mm256_storeu_si256(pvHStore + j, vH);

                /* Update vE value. */
                vH = _mm256_add_epi16(vH, vGapO); /* saturation arithmetic, result >= 0 */
                e = _mm256_add_epi16(e, vGapE);
                e = _mm256_max_epi16(e, vH);
                e = _mm256_max_epi16(e, vGapO);
                _mm256_storeu_si256(pvE + j, e);

                /* Update vF value. */
//            vF = _mm256_add_epi16(vF, vGapE);
//            vF = _mm256_max_epi16(vF, vH);

                /* Load the next vH. */
                vH = _mm256_loadu_si256(pvHLoad + j);

            }
//        std::cout << "SIMD line 294" << std::endl;
            // put the score into the score matrix begin
            for( j=0; j<_segLen; ++j ){
                a = ( int16_t*)(pvHStore+j);
                for( z=0; z<16; ++z ){
                    if( (j + z*_segLen) < readLen ){
                        _similarity_matrix[j+z*_segLen+1] = a[z];
                    }
                }
            }
//        std::cout << "SIMD line 304" << std::endl;
            deletionScore = 0;
            for( j=0; j<readLen; ++j ){
                deletionScore = _similarity_matrix[j] + (int16_t)(_open_gap_penalty) > deletionScore + _extend_gap_penalty ? _similarity_matrix[j] + (int16_t)(_open_gap_penalty) : deletionScore + _extend_gap_penalty ;
                if( deletionScore > _similarity_matrix[j+1] ){
                    _similarity_matrix[j+1] = deletionScore;
                }
                deletionScore = deletionScore > 0 ? deletionScore : 0;
                if( _similarity_matrix[j+1] > maxScore ){
                    endPosition1 = i;
                    endPosition2 = j;
                    maxScore = _similarity_matrix[j+1];
                }
            }
//        std::cout << "SIMD line 315" << std::endl;
            //updata vH begin
            for( j=0; j<_segLen; j++ ){
                a = (int16_t*)(pvHStore+j);
                for( z=0; z<16; ++z ){
                    if( (j + z*_segLen) < readLen ){
                        a[z] = _similarity_matrix[j+z*_segLen+1];
                    }
                }
            }//updata vH end
        }
        free(pvHStore);
        free(pvHLoad);
        free(profile);
        free(pvE);
    }else{
        SmithWaterman_non_avx(ref, read, refLen, readLen, _open_gap_penalty, _extend_gap_penalty,
                                   maxScore, endPosition1, endPosition2,  scorei,  positions);
    }
#else
    SmithWaterman_non_avx(ref, read, refLen, readLen, _open_gap_penalty, _extend_gap_penalty,
                                   maxScore, endPosition1, endPosition2,  scorei,  positions);
#endif
}






// this is a 2-piece affine gap cost smithwaterman method which return cigar
std::vector<uint32_t> SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                    const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                    const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                    int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                                    bool reverseAlignment, bool returnCigar){

    int32_t ** M = new int32_t *[length1 + 1];
    int32_t ** E1 = new int32_t *[length1 + 1];
    int32_t ** E2 = new int32_t *[length1 + 1];
    int32_t ** F1 = new int32_t *[length1 + 1];
    int32_t ** F2 = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
        M[i] = new int32_t [length2 + 1];
        E1[i] = new int32_t [length2 + 1];
        E2[i] = new int32_t [length2 + 1];
        F1[i] = new int32_t [length2 + 1];
        F2[i] = new int32_t [length2 + 1];
        std::fill_n(M[i], length2+1, 0);
        std::fill_n(E1[i], length2+1, 0);
        std::fill_n(E2[i], length2+1, 0);
        std::fill_n(F1[i], length2+1, 0);
        std::fill_n(F2[i], length2+1, 0);
    }

    endPosition1=0;
    endPosition2=0;
    maxScore = 0;
    for ( i=1; i<=length1; ++i ){
        for ( j=1; j<=length2; ++j ){
            E1[i][j] = (_open_gap_penalty1 + M[i][j-1]) > (_extend_gap_penalty1 + E1[i][j-1]) ? (_open_gap_penalty1 + M[i][j-1]) : (_extend_gap_penalty1 + E1[i][j-1]);
            E2[i][j] = (_open_gap_penalty2 + M[i][j-1]) > (_extend_gap_penalty2 + E2[i][j-1]) ? (_open_gap_penalty2 + M[i][j-1]) : (_extend_gap_penalty2 + E2[i][j-1]);

            F1[i][j] = (_open_gap_penalty1 + M[i-1][j]) > (_extend_gap_penalty1 + F1[i-1][j]) ? (_open_gap_penalty1 + M[i-1][j]) : (_extend_gap_penalty1 + F1[i-1][j]);
            F2[i][j] = (_open_gap_penalty2 + M[i-1][j]) >  (_extend_gap_penalty2 + F2[i-1][j]) ?(_open_gap_penalty2 + M[i-1][j]) : (_extend_gap_penalty2 + F2[i-1][j]);

            M[i][j] = (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) > E1[i][j] ? (m.getScore(seq1[i-1], seq2[j-1]) + M[i-1][j-1]) : E1[i][j];
            M[i][j] = M[i][j] > F1[i][j] ? M[i][j] : F1[i][j];
            M[i][j] = M[i][j] > F2[i][j] ? M[i][j] : F2[i][j];
            M[i][j] = M[i][j] > E2[i][j] ? M[i][j] : E2[i][j];

            M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
            F1[i][j] = F1[i][j] > 0 ? F1[i][j] : 0;
            E1[i][j] = E1[i][j] > 0 ? E1[i][j] : 0;
            F2[i][j] = F2[i][j] > 0 ? F2[i][j] : 0;
            E2[i][j] = E2[i][j] > 0 ? E2[i][j] : 0;

            if( M[i][j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                maxScore = M[i][j];
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
    }
    //std::cout << "maxScore " << maxScore << std::endl;
    std::vector<uint32_t> cigar;

    if( returnCigar ){ // cigar is only needed for reverse alignment
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int32_t ii = endPosition1;
        int32_t jj = endPosition2;
        while (ii>0 && jj>0) {
            if (M[ii][jj] == M[ii - 1][jj - 1] + m.getScore(seq1[ii - 1], seq2[jj - 1]) ) {
                --ii;
                --jj;
                op = 0;
            } else if ( M[ii][jj] == E1[ii][jj]) {
                --jj;
                op = 1;
            } else if ( M[ii][jj] == F1[ii][jj]) {
                --ii;
                op = 2;
            } else if ( M[ii][jj] == F2[ii][jj]) {
                --ii;
                op = 2;
            } else if ( M[ii][jj] == E2[ii][jj]) {
                --jj;
                op = 1;
            } else if(M[ii][jj]==0){
                //std::cerr << "there is something wrong with smith-waterman algorithm in line 110" << std::endl;
                // here should never run, there is some problem with the code
                break;
            } else {
//                std::cout << "M[ii][jj] " << M[ii][jj] << " E[ii][jj] " << E[ii][jj] << " F[ii][jj] " << F[ii][jj] << std::endl;
                std::cout << "there is something wrong with smith-waterman algorithm in line 114" << std::endl;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        if (!reverseAlignment){
            std::reverse(cigar.begin(),cigar.end());
        }
    }

    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E1[i];
        delete[] F1[i];
        delete[] E2[i];
        delete[] F2[i];
    }
    delete[] M;
    delete[] E1;
    delete[] F1;
    delete[] E2;
    delete[] F2;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}


// this is a 2-piece affine gap cost smith-waterman method which give maximum value and positions
void SmithWaterman(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                   const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                   const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                   int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m){

    int32_t * M1 = new int32_t [length2 + 1];
    int32_t * F11 = new int32_t [length2 + 1];
    int32_t * F12 = new int32_t [length2 + 1];

    int32_t * M2 = new int32_t [length2 + 1];
    int32_t * F21 = new int32_t [length2 + 1];
    int32_t * F22 = new int32_t [length2 + 1];

    int32_t * t;

    std::fill_n(M1, length2+1, 0);
    std::fill_n(F11, length2+1, 0);
    std::fill_n(F12, length2+1, 0);

    std::fill_n(M2, length2+1, 0);
    std::fill_n(F21, length2+1, 0);
    std::fill_n(F22, length2+1, 0);

    endPosition1=0;
    endPosition2=0;
    int32_t i, j, e1, e2;
    maxScore = 0;
    int32_t mscore;

    for ( i=1; i<=length1; ++i ){
        e1=0, e2=0;
        for ( j=1; j<=length2; ++j ){
            e1 = (_open_gap_penalty1 + M2[j-1]) > (_extend_gap_penalty1 + e1) ? (_open_gap_penalty1 + M2[j-1]) : (_extend_gap_penalty1 + e1);
            F21[j] = (_open_gap_penalty1 + M1[j]) >  (_extend_gap_penalty1 + F11[j]) ?(_open_gap_penalty1 + M1[j]) :  (_extend_gap_penalty1 + F11[j]);

            e2 = (_open_gap_penalty2 + M2[j-1]) > (_extend_gap_penalty2 + e2) ? (_open_gap_penalty2 + M2[j-1]) : (_extend_gap_penalty2 + e2);
            F22[j] = (_open_gap_penalty2 + M1[j]) >  (_extend_gap_penalty2 + F12[j]) ?(_open_gap_penalty2 + M1[j]) :  (_extend_gap_penalty2 + F12[j]);

            mscore = (m.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) > e1 ? (m.getScore(seq1[i-1], seq2[j-1]) + M1[j-1]) : e1;
            mscore = mscore > F21[j] ? mscore : F21[j];
            mscore = mscore > F22[j] ? mscore : F22[j];
            mscore = mscore > e2 ? mscore :e2;

            M2[j] = mscore > 0 ? mscore : 0;
            F21[j] = F21[j] > 0 ? F21[j] : 0;
            e1 = e1 > 0 ? e1 : 0;
            F22[j] = F22[j] > 0 ? F22[j] : 0;
            e2 = e2 > 0 ? e2 : 0;

            if( mscore > maxScore ){ // please do not change > to >=, since we are doing local alignment
                maxScore = mscore;
                endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2=j;
            }
        }
        t = M1;
        M1 = M2;
        M2 = t;

        t = F11;
        F11 = F21;
        F21 = t;

        t = F12;
        F12 = F22;
        F22 = t;
    }

    delete[] M1;
    delete[] F11;
    delete[] F12;

    delete[] M2;
    delete[] F21;
    delete[] F22;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
}



//standard implementing of SemiGlobal sequence alignment approach
std::vector<uint32_t> SemiGlobal_no_avx(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                 const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                                 bool returnCigar, const int32_t & zdrop, const int32_t & w, Matrix & T){
    int32_t i=0, j=0, mscore=0;

    uint8_t d=0;
    if( returnCigar ) {
        T.reset(length1+1, length2+1);
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F1 = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2
    int32_t * F2 = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(M2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = -1;
    endPosition2 = -1;

    maxScore = 0;
    int32_t e1, e2, start, end, thisMax, thisMaxj, l;
    for ( i=1; i<=length1; ++i ){
        thisMax = 0;
        thisMaxj = 0;

        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;

        e1 = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        e2 = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F1[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F2[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;

        if( returnCigar ) {
            for (j = start; j < end; ++j) {
                mscore = m.getScore(seq1[i - 1], seq2[j - 1]) + M1[j - 1];

                d = mscore > F1[j] ? 0 : 1;
                mscore = mscore > F1[j] ? mscore : F1[j];

                d = mscore > e1 ? d : 2;
                mscore = mscore > e1 ? mscore : e1;

                d = mscore > F2[j] ? d : 3;
                mscore = mscore > F2[j] ? mscore : F2[j];

                d = mscore > e2 ? d : 4;
                mscore = mscore > e2 ? mscore : e2;

                M2[j] = mscore;
                if ( mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }

                int32_t h = mscore + _open_gap_penalty1;
                int32_t f = F1[j] + _extend_gap_penalty1;
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F1[j] = f;

                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;

                int32_t h2 = mscore + _open_gap_penalty2;
                int32_t f2 = F2[j] + _extend_gap_penalty2;
                d |= f2 >= h2 ? 1 << 5 : 0;
                f2 = f2 >= h2 ? f2 : h2;
                F2[j] = f2;

                e2 += _extend_gap_penalty2;
                d |= e2 >= h2 ? 1 << 6 : 0;
                e2 = e2 >= h2 ? e2 : h2;

                T.set(i, j, d);
            }
            j=end;
            {
                mscore = m.getScore(seq1[i - 1], seq2[j - 1]) + M1[j - 1];

                d = mscore > F1[j] ? 0 : 1;
                mscore = mscore > F1[j] ? mscore : F1[j];

                d = mscore > e1 ? d : 2;
                mscore = mscore > e1 ? mscore : e1;

                d = mscore > F2[j] ? d : 3;
                mscore = mscore > F2[j] ? mscore : F2[j];

                d = mscore > e2 ? d : 4;
                mscore = mscore > e2 ? mscore : e2;

                M2[j] = mscore;
                if ( mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }
            }
        }else{
            for ( j=start; j<end; ++j ){
                mscore = m.getScore(seq1[i - 1], seq2[j - 1]) + M1[j - 1];

                mscore = mscore > F1[j] ? mscore : F1[j];

                mscore = mscore > e1 ? mscore : e1;

                mscore = mscore > F2[j] ? mscore : F2[j];

                mscore = mscore > e2 ? mscore : e2;

                M2[j] = mscore;
                if ( mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }
                int32_t h = mscore + _open_gap_penalty1;
                int32_t f = F1[j] + _extend_gap_penalty1;
                f = f >= h ? f : h;
                F1[j] = f;

                e1 += _extend_gap_penalty1;
                e1 = e1 >= h ? e1 : h;

                int32_t h2 = mscore + _open_gap_penalty2;
                int32_t f2 = F2[j] + _extend_gap_penalty2;
                f2 = f2 >= h2 ? f2 : h2;
                F2[j] = f2;

                e2 += _extend_gap_penalty2;
                e2 = e2 >= h2 ? e2 : h2;
            }
            j = end;
            {
                mscore = m.getScore(seq1[i - 1], seq2[j - 1]) + M1[j - 1];

                mscore = mscore > F1[j] ? mscore : F1[j];

                mscore = mscore > e1 ? mscore : e1;

                mscore = mscore > F2[j] ? mscore : F2[j];

                mscore = mscore > e2 ? mscore : e2;

                M2[j] = mscore;
                if ( mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }
            }
        }

        if( thisMax > maxScore ){ // please do not change > to >=, since we are doing local alignment
            // >= will omit the first similar fragments
            maxScore = thisMax;
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=thisMaxj;
        }//else if( thisMaxj > endPosition2 && i > endPosition1 ){
            l = (i - endPosition1) - (thisMaxj-endPosition2);
            l = l>0 ? l : -l;
            if( maxScore-thisMax > zdrop + l * _extend_gap_penalty2) {
                break;
            }
        //}
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
        M1[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        M2[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
    }

    std::vector<uint32_t> cigar;
    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int32_t ii = endPosition1;
        int32_t jj = endPosition2;
        int tmp, state=0;
        while (ii>0 && jj>0) {
            tmp = T.get(ii, jj);
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        while( ii>0 ){
            op =2;
            --ii;
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        while( jj>0 ){
            op =1;
            --jj;
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
    }

//    std::cout << " line 511 maxScore " << maxScore << std::endl;
    delete[] M1;
    delete[] M2;
    delete[] F2;
    delete[] F1;

    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
//    std::cout << "line 889 endPosition1 " << endPosition1 << " endPosition2: " << endPosition2 << " maxScore: " << maxScore << std::endl;

    return cigar;
}



std::vector<uint32_t> SemiGlobal(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                        const int32_t &length2, const int &_open_gap_penalty1, const int &_extend_gap_penalty1,
                                        const int &_open_gap_penalty2, const int &_extend_gap_penalty2,
                                        int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                                        bool returnCigar, const int32_t & zdrop, const int32_t & w, Matrix & T){
   return SemiGlobal_no_avx(seq1, seq2, length1, length2, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2,
           _extend_gap_penalty2, maxScore, endPosition1, endPosition2, m, returnCigar, zdrop, w, T);


}




// this is a gapless alignment extension approach, and do not return cigar
void gapless_extend(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                        const int32_t &length2, int32_t &maxScore,
                        int32_t &endPosition1, int32_t &endPosition2, const Scorei & m){// not sure this band approach and avx which is faster

    int32_t length = length2 < length1 ? length2 : length1;
    int32_t i, j, currentScore=maxScore, mscore;
    int32_t thisMaxi=0;
    for ( i=1; i<=length; ++i ){
        currentScore += m.getScore(seq1[i-1], seq2[i-1]);
        if( currentScore > maxScore ){
            maxScore = currentScore;
            thisMaxi = i;
        }else if( currentScore < maxScore + 3* m.getScore(1, 2)){ //m.getScore(1, 2) is a mismatch
            break;
        }
    }
    endPosition1 += thisMaxi;
    endPosition2 += thisMaxi;
}
/*
void gapless_extend(int8_t *seq1, int8_t *seq2, int32_t &maxScore, const int8_t & yDrop,
                    int32_t & start1, int32_t & end1, int32_t & start2, int32_t & end2, const Scorei & m){// not sure this band approach and avx which is faster

    int32_t length = length2 < length1 ? length2 : length1;
    int32_t i, j, currentScore=maxScore, mscore;
    int32_t thisMaxi=0;
    for ( i=1; i<=length; ++i ){
        currentScore += m.getScore(seq1[i-1], seq2[i-1]);
        if( currentScore > maxScore ){
            maxScore = currentScore;
            thisMaxi = i;
        }else if( currentScore < maxScore + 3* m.getScore(1, 2)){ //m.getScore(1, 2) is a mismatch
            break;
        }
    }
    endPosition1 += thisMaxi;
    endPosition2 += thisMaxi;
}
*/


// this is a x-drop sequence alignment extension approach, and do not return cigar
void SemiGlobal_xextend(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                                 const int32_t & xdrop, const int32_t & w){// not sure this band approach and avx which is faster
    int32_t i, j, mscore;
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(M2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e;
    int32_t thisMaxj;
    int32_t thisMax;
    int32_t start;
    int32_t end;
    for ( i=1; i<=length1; ++i ){
        thisMax = 0;
        thisMaxj=-1;

        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;

        e = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        for ( j=start; j<=end; ++j ){
            mscore = m.getScore(seq1[i-1], seq2[j-1]) + M1[j-1];
            mscore = mscore > F[j] ? mscore : F[j];
            mscore = mscore > e ? mscore : e;
            M2[j] = mscore;
            if (mscore > thisMax){
                thisMax = mscore;
                thisMaxj = j;
            }
            int32_t h = mscore + _open_gap_penalty;
            int32_t f = F[j] + _extend_gap_penalty;
            f  = f >= h? f    : h;
            F[j] = f;
            e += _extend_gap_penalty;
            e  = e >= h? e  : h;
        }
        if( thisMaxj != length2 && thisMaxj != 1 && (thisMaxj == start ||  thisMaxj == end) ){
            int32_t thisw = 2*w;
            SemiGlobal_xextend(seq1, seq2, length1, length2, _open_gap_penalty, _extend_gap_penalty, maxScore,
                    endPosition1, endPosition2, m, xdrop, thisw);
            return;
        }
        if( thisMax > maxScore ){ // please do not change > to >=, since we are doing local alignment
            // >= will omit the first similar fragments
            maxScore = thisMax;
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=thisMaxj;
        }//else if( thisMaxj > endPosition2 && i > endPosition1 ){
            if( maxScore-thisMax > xdrop ) {
                break;
            }
        //}
        if ( thisMax<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
        M1[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        M2[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
    }
    delete[] M1;
    delete[] M2;
    delete[] F;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
}


//same with above one, but return cigar
std::vector<uint32_t> SemiGlobal_xextend(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                        const int32_t &length2, const int &_open_gap_penalty, const int &_extend_gap_penalty,
                        int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2, const Scorei & m,
                        const int32_t & xdrop, const int32_t & w, Matrix & T, int32_t & iii, int32_t  &jjj){

    int32_t i=0, j=0, mscore=0, temp;

    uint8_t d=0;

    T.reset(length1+1, length2+1);

    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(M2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = 0;
    endPosition2 = 0;

    maxScore = 0;
    int32_t e;
    int32_t thisMaxj, thisMax, start, end, h, f;
    for ( i=1; i<=length1; ++i ){
        thisMax = 0;
        thisMaxj=-1;

        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;

        e = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        for (j = start; j <= end; ++j) {
            mscore = m.getScore(seq1[i-1], seq2[j-1]) + M1[j - 1];

            d = mscore > F[j] ? 0 : 1;
            mscore = mscore > F[j] ? mscore : F[j];

            d = mscore > e ? d : 2;
            mscore = mscore > e ? mscore : e;
            M2[j] = mscore;
            if (mscore > thisMax){
                thisMax = mscore;
                thisMaxj = j;
            }

            h = mscore + _open_gap_penalty;
            f = F[j] + _extend_gap_penalty;
            d |= f >= h ? 1 << 3 : 0;
            f = f >= h ? f : h;
            F[j] = f;

            e += _extend_gap_penalty;
            d |= e >= h ? 1 << 4 : 0;
            e = e >= h ? e : h;

            T.set(i, j, d);
        }
        if( thisMaxj != length2 && thisMaxj != 1 && (thisMaxj == start ||  thisMaxj == end) ){
            int32_t thisw = 2*w;
            return SemiGlobal_xextend(seq1, seq2, length1, length2, _open_gap_penalty, _extend_gap_penalty, maxScore,
                    endPosition1, endPosition2, m, xdrop, thisw, T, iii, jjj);
        }

        if( thisMax > maxScore ){ // please do not change > to >=, since we are doing local alignment
            // >= will omit the first similar fragments
            maxScore = thisMax;
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=thisMaxj;
        }//else if( thisMaxj > endPosition2 && i > endPosition1 ){   // the idea here is, even the maximum is too small, all the values in the desired area are too small
            if( maxScore-thisMax > xdrop) {
                break;
            }
        //}
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
        M1[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        M2[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
    }
//    std::cout << "line 861 endPosition1 " << endPosition1 << " endPosition2: " << endPosition2 << " maxScore: " << maxScore << std::endl;
    std::vector<uint32_t> cigar;

    uint32_t op = 0;
    uint32_t length = 1;
    // trace back begin
    // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
    int32_t ii = endPosition1;
    int32_t jj = endPosition2;
    int tmp, state=0;
    while (ii>0 && jj>0) {
        tmp = T.get(ii, jj);
        if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
        else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
        if (state == 0) state = tmp & 7;
        if (state == 0){
            op=0;
            --ii;
            --jj;
        }else if (state == 1 || state == 3){
            op =2;
            --ii;
        }else{
            op =1;
            --jj;
        }

        if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
            cigar.push_back(length << 4 | op);
        }else{
            cigar[cigar.size()-1] += length<<4;
        }
    }

//    std::reverse(cigar.begin(),cigar.end());

    delete[] M1;
    delete[] M2;
    delete[] F;

    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    iii=ii;
    jjj=jj;
    return cigar;
}




/*************************************************************************

a weighted dynamic programming sequence alignment method ZDP (Zebric dynamic programming)

The standard needleman-wunsch algorithm set up a score matrix and generate alignment by tracing back
Since we have different score weight for each base-pair, the tracing back step would query the weighted score again, which is time consuming.
Here when set up the score matrix, we setup a trace matrix at the same time, and we could trace back using the trace matrix
        (this is a kind of approximation, might could not get the identical alignment result comparing the the standard traceback approach,
while the standard traceback might not work for the weighted sequence alignment approach, since different weight could generate indistinct trace path)

************************************************************************/



// seqAchar is the reference sequence
// seqBChar is the query sequence

// lengthA the the length of the reference sequence (if the sequence is a java/cpp string, this variable is not necessary)
// lengthA the the length of the query sequence

// weight is a vector of weight of each reference sequence bp
// Score is a custom class, together with the 'weight' variable,  we could get difference match/mis-match/insertion/deletion scores for each base-pair

// A is the alignment result using CIGAR letters [MDI]


// this is the semiglobal weighted sequence alignment approach
std::vector<uint32_t> SemiGlobal(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int16_t * weights, Score & score,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                                 bool returnCigar, const int32_t & z, const int32_t & w, Matrix & T){
    int32_t i, j, mscore;
    uint8_t d;
    if( returnCigar ) {
        T.reset(length1+1, length2+1);
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F1 = new int32_t [length2 + 1]; //F1 and F2 is for gapscore1 and gapscore2
    int32_t * F2 = new int32_t [length2 + 1];

    std::fill_n(M1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(M2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = -1;
    endPosition2 = -1;

    maxScore = 0;
    int32_t e1, e2, zdrop, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2;
    int32_t ** m;
    int32_t thisMax;
    int32_t thisMaxj;
    int32_t start;
    int32_t end;
    int32_t l;
    for ( i=1; i<=length1; ++i ){
        thisMax = 0;
        thisMaxj = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;
        e1 = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        e2 = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F1[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        F2[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        if( weights[i-1] == weights[i-2] && i>1 ){

        }else{
            _open_gap_penalty1 = score.getOpenGapPenalty1(weights[i-1]);
            _extend_gap_penalty1 = score.getExtendGapPenalty1(weights[i-1]);
            _open_gap_penalty2 = score.getOpenGapPenalty2(weights[i-1]);
            _extend_gap_penalty2 = score.getExtendGapPenalty2(weights[i-1]);
            m = score.getM(weights[i-1]);
            //zdrop = z * (weights[i-1]/2); // 1/2===1
            zdrop = score.getZdrop(weights[i-1]);
        }

        if( returnCigar ) {
            for (j = start; j <= end; ++j) {
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                d = mscore > F1[j] ? 0 : 1;
                mscore = mscore > F1[j] ? mscore : F1[j];

                d = mscore > e1 ? d : 2;
                mscore = mscore > e1 ? mscore : e1;

                d = mscore > F2[j] ? d : 3;
                mscore = mscore > F2[j] ? mscore : F2[j];

                d = mscore > e2 ? d : 4;
                mscore = mscore > e2 ? mscore : e2;
                M2[j]=mscore;
                if (mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }

                int32_t h = mscore + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F1[j] + score.getExtendGapPenalty1(weights[i]);
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F1[j] = f;

                h = mscore + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;

                int32_t h2 = mscore + score.getOpenGapPenalty2(weights[i]);
                int32_t f2 = F2[j] + score.getExtendGapPenalty2(weights[i]);
                d |= f2 >= h2 ? 1 << 5 : 0;
                f2 = f2 >= h2 ? f2 : h2;
                F2[j] = f2;

                h2 = mscore + _open_gap_penalty2;
                e2 += _extend_gap_penalty2;
                d |= e2 >= h2 ? 1 << 6 : 0;
                e2 = e2 >= h2 ? e2 : h2;
                T.set(i, j, d);
            }
        }else{
            for ( j=start; j<=end; ++j ){
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                mscore = mscore > F1[j] ? mscore : F1[j];

                mscore = mscore > e1 ? mscore : e1;

                mscore = mscore > F2[j] ? mscore : F2[j];

                mscore = mscore > e2 ? mscore : e2;
                M2[j]=mscore;

                if ( mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }

                int32_t h = mscore + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F1[j] + score.getExtendGapPenalty1(weights[i]);
                f  = f >= h? f    : h;
                F1[j] = f;

                h = mscore + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                e1  = e1 >= h? e1  : h;

                int32_t h2 = mscore + score.getOpenGapPenalty2(weights[i]);
                int32_t f2 = F2[j] + score.getExtendGapPenalty2(weights[i]);
                f2 = f2 >= h2? f2 : h2;
                F2[j] = f2;

                h2 = mscore + _open_gap_penalty2;
                e2 += _extend_gap_penalty2;
                e2 = e2 >= h2? e2 : h2;
            }
        }

        if( thisMax > maxScore ){ // please do not change > to >=, since we are doing local alignment
            // >= will omit the first similar fragments
            maxScore = thisMax;
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=thisMaxj;
        }//else if ( thisMaxj > endPosition2 && i > endPosition1 ){
            l = (i - endPosition1) - (thisMaxj-endPosition2);
            l = l>0 ? l : -l;
            if( maxScore-thisMax > zdrop + l * _extend_gap_penalty2) {
                break;
            }
       // }
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
        M1[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        M2[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
    }

    std::vector<uint32_t> cigar;

    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int32_t ii = endPosition1;
        int32_t jj = endPosition2;
        int tmp, state=0;

        while (ii>0 && jj>0) {
            tmp = T.get(ii, jj);
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        while( ii>0 ){
            op =2;
            --ii;
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        while( jj>0 ){
            op =1;
            --jj;
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
    }

    delete[] M1;
    delete[] M2;
    delete[] F2;
    delete[] F1;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}



// this is the semiglobal weighted sequence alignment approach
// similar with the above function by use single gap penalty
std::vector<uint32_t> SemiGlobal_single_gap_penalty(int8_t *seq1, int8_t *seq2, const int32_t &length1,
                                 const int32_t &length2, const int16_t * weights, Score & score,
                                 int32_t &maxScore, int32_t &endPosition1, int32_t &endPosition2,
                                 bool returnCigar, const int32_t & z, const int32_t & w, Matrix & T){
    int32_t i, j, mscore;
    uint8_t d;
    if( returnCigar ) {
        T.reset(length1 + 1, length2 + 1);
    }
    int32_t * t;
    int32_t * M1 = new int32_t [length2 + 1]; //M1 and M2 is for the previous column and the current column
    int32_t * M2 = new int32_t [length2 + 1];

    int32_t * F = new int32_t [length2 + 1]; //F

    std::fill_n(M1, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(M2, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    std::fill_n(F, length2+1, SCORE_OUT_BANDED_ALIGNMENT_REGION);
    M1[0]=0;
    M2[0]=0;
    endPosition1 = -1;
    endPosition2 = -1;

    maxScore = 0;
    int32_t e1, zdrop, _open_gap_penalty1, _extend_gap_penalty1;
    int32_t ** m;

    int32_t start;
    int32_t end;

    for ( i=1; i<=length1; ++i ){
        int32_t thisMax = 0;
        int32_t thisMaxj=-1;

        start = i - w > 1 ? i - w : 1;
        end = i + w < length2 ?  i + w : length2;

        e1 = SCORE_OUT_BANDED_ALIGNMENT_REGION;

        F[end] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        if( weights[i-1] == weights[i-2] && i>1 ){  // if this one is same with last one, we do not need to change parameters

        }else{
            _open_gap_penalty1 = score.getOpenGapPenalty1(weights[i-1]);
            _extend_gap_penalty1 = score.getExtendGapPenalty1(weights[i-1]);
            m = score.getM(weights[i-1]);
            zdrop = score.getZdrop(weights[i-1]);
        }


        if( returnCigar ) {
            for (j = start; j <= end; ++j) {
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];

                d = mscore > F[j] ? 0 : 1;
                mscore = mscore > F[j] ? mscore : F[j];

                d = mscore > e1 ? d : 2;
                mscore = mscore > e1 ?  mscore : e1;
                M2[j]=mscore;
                if (mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }

                int32_t h = mscore + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F[j] + score.getExtendGapPenalty1(weights[i]);
                d |= f >= h ? 1 << 3 : 0;
                f = f >= h ? f : h;
                F[j] = f;

                h = mscore + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                d |= e1 >= h ? 1 << 4 : 0;
                e1 = e1 >= h ? e1 : h;
                T.set(i, j ,d);
            }
        }else{
            for ( j=start; j<=end; ++j ){
                mscore = m[seq1[i - 1]][seq2[j - 1]] + M1[j - 1];
                mscore = mscore > F[j] ? mscore : F[j];
                mscore = mscore > e1 ? mscore : e1;
                M2[j]=mscore;
                if ( mscore > thisMax){
                    thisMax = mscore;
                    thisMaxj = j;
                }

                int32_t h = mscore + score.getOpenGapPenalty1(weights[i]);
                int32_t f = F[j] + score.getExtendGapPenalty1(weights[i]);
                f  = f >= h? f    : h;
                F[j] = f;

                h = mscore + _open_gap_penalty1;
                e1 += _extend_gap_penalty1;
                e1  = e1 >= h? e1  : h;
            }
        }
        if( thisMax > maxScore ){ // please do not change > to >=, since we are doing local alignment
            // >= will omit the first similar fragments
            maxScore = thisMax;
            endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
            endPosition2=thisMaxj;
        }//else if( thisMaxj > endPosition2 && i > endPosition1 ){
            if( maxScore-thisMax > zdrop ) {
                break;
            }
        //}
        if ( thisMaxj<=0 ){
            break;
        }
        t = M1;
        M1 = M2;
        M2 = t;
        M1[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
        M2[0] = SCORE_OUT_BANDED_ALIGNMENT_REGION;
    }

    std::vector<uint32_t> cigar;

    if( returnCigar ){
        uint32_t op = 0;
        uint32_t length = 1;
        // trace back begin
        // For speed up and RAM saving purpose, this is just am approximation tracing back implementation
        int32_t ii = endPosition1;
        int32_t jj = endPosition2;
        int tmp, state=0;

        while (ii>0 && jj>0) {
            tmp = T.get(ii, jj);
            if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
            else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
            if (state == 0) state = tmp & 7;
            if (state == 0){
                op=0;
                --ii;
                --jj;
            }else if (state == 1 || state == 3){
                op =2;
                --ii;
            }else {
                op =1;
                --jj;
            }

            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        while( ii>0 ){
            op =2;
            --ii;
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        while( jj>0 ){
            op =1;
            --jj;
            if( cigar.size() == 0 || op != (cigar[cigar.size() - 1]&0xf) ){
                cigar.push_back(length << 4 | op);
            }else{
                cigar[cigar.size()-1] += length<<4;
            }
        }
        // trace back end
        std::reverse(cigar.begin(),cigar.end());
    }

    delete[] M1;
    delete[] M2;
    delete[] F;
    --endPosition1;// change the endPosition1 and endPosition2 to 0 based coordinate
    --endPosition2;
    return cigar;
}
