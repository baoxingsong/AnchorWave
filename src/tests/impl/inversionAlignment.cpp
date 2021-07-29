//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include <string>
#include <iostream>
#include <sstream>
#include "../../myImportandFunction/alignSlidingWindow.h"
#include "../../minimap2/minimap.h"
#include "../../minimap2/kseq.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

TEST(inversionAlignment, c1){ // just to make sure that every line has been analysed


    std::string dna_d = "GAAGTACTAATCGATGTT";
    std::string dna_q = "GTAAGTACAACATCGATTA";

    std::string _alignment_q;
    std::string _alignment_d;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -2;
    int32_t openGapPenalty1 = -3;
    int32_t extendGapPenalty1 = -1;
    int32_t invertyPenalty = 0;

    int32_t maxDnaDIndex = 0;

    inversionAlignment( dna_q, dna_d, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, invertyPenalty, maxDnaDIndex );
    std::cout << "maxDnaDIndex:\t" << maxDnaDIndex << std::endl;
    int32_t length1 = dna_d.length();
    std::cout << dna_d.substr (0, maxDnaDIndex) << "\t" << getReverseComplementary(dna_d.substr (maxDnaDIndex, length1-maxDnaDIndex)) << std::endl;
    std::string _dna_d2 = dna_d.substr (0, maxDnaDIndex) + getReverseComplementary(dna_d.substr (maxDnaDIndex, length1-maxDnaDIndex));

    double inversion_PENALTY = -1;
    double MIN_ALIGNMENT_SCORE = 3;
    bool considerInversion = false;
    int32_t wfaSize = 40000;
    int32_t wfaSize2 = 100000;
    int32_t slidingWindowSize = 40000;
    int expectedCopies = 1;
    double maximumSimilarity = 0.6;
    int32_t min_wavefront_length = 20;
    int32_t max_distance_threshold = 100;

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 30;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    Scorei m(matchingScore, mismatchingPenalty);

    affine_penalties_t affine_penalties = {
            .match = 0,
            .mismatch = -mismatchingPenalty+matchingScore,
            .gap_opening = -openGapPenalty1+matchingScore,
            .gap_extension = -extendGapPenalty1+matchingScore,
    };
    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_512M);
    std::cout << _dna_d2 << std::endl << dna_q << std::endl << std::endl;

    int32_t score = alignSlidingWindow(  dna_q,  _dna_d2, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    std::cout << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;

    mm_allocator_delete(mm_allocator);
    ASSERT_EQ(0, 0);
}


//maxDnaDIndex:	7
//maxDnaQIndex:	0
//GAAGTAC	AACATCGATTA
//GAAGTACAACATCGATTA
//GTAAGTACAACATCGATTA
//
//G-AAGTACAACATCGATTA
//GTAAGTACAACATCGATTA
//32

//std::string dna_d = "GAAGTACTAATCGATGTT";
//std::string dna_q = "GTAAGTACAACATCGATTA";


TEST(inversionAlignment2, c1){ // just to make sure that every line has been analysed

//    std::string dna_d = "GAAGTACTAATCGATGTTGTCACT";
//    std::string dna_q = "GTAAGTACAACATCGATTAGTCACT";

    std::string dna_d = "GAAGTACTAATCGATGTT";
    std::string dna_q = "GTAAGTACAACATCGATTA";

    std::string _alignment_q;
    std::string _alignment_d;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -2;
    int32_t openGapPenalty1 = -3;
    int32_t extendGapPenalty1 = -1;
    int32_t invertyPenalty = -2;

    int32_t maxDnaDIndex = 0;
    int32_t maxDnaDIndexe = 0;
    inversionAlignment( dna_q, dna_d, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, invertyPenalty, maxDnaDIndex, maxDnaDIndexe);
    std::cout << "maxDnaDIndex:\t" << maxDnaDIndex << std::endl;
    std::cout << "maxDnaQIndexe:\t" << maxDnaDIndexe << std::endl;
    int32_t length1 = dna_d.length();

    std::cout << dna_d.substr (0, maxDnaDIndex) << "\t" << getReverseComplementary(dna_d.substr (maxDnaDIndex, maxDnaDIndexe-maxDnaDIndex)) << "\t" <<  dna_d.substr(maxDnaDIndexe, length1-maxDnaDIndexe) << std::endl;

    std::string _dna_d2 = dna_d.substr (0, maxDnaDIndex) + getReverseComplementary(dna_d.substr (maxDnaDIndex, maxDnaDIndexe-maxDnaDIndex)) +  dna_d.substr(maxDnaDIndexe, length1-maxDnaDIndexe);
    if( maxDnaDIndexe == length1 ){
        _dna_d2 = dna_d.substr (0, maxDnaDIndex) + getReverseComplementary(dna_d.substr (maxDnaDIndex, maxDnaDIndexe-maxDnaDIndex));
    }
    std::cout << dna_d << std::endl;
    double inversion_PENALTY = -1;
    double MIN_ALIGNMENT_SCORE = 3;
    bool considerInversion = false;
    int32_t wfaSize = 40000;
    int32_t wfaSize2 = 100000;
    int32_t slidingWindowSize = 40000;
    int expectedCopies = 1;
    double maximumSimilarity = 0.6;
    int32_t min_wavefront_length = 20;
    int32_t max_distance_threshold = 100;

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 30;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    Scorei m(matchingScore, mismatchingPenalty);

    affine_penalties_t affine_penalties = {
            .match = 0,
            .mismatch = -mismatchingPenalty+matchingScore,
            .gap_opening = -openGapPenalty1+matchingScore,
            .gap_extension = -extendGapPenalty1+matchingScore,
    };
    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_512M);
    std::cout << _dna_d2 << std::endl << dna_q << std::endl << std::endl;

    int32_t score = alignSlidingWindow(  dna_q,  _dna_d2, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    std::cout << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;

    mm_allocator_delete(mm_allocator);
    ASSERT_EQ(0, 0);
}



TEST(inversionAlignment2, c2){ // just to make sure that every line has been analysed

//    std::string dna_d = "GAAGTACTAATCGATGTTGTCACT";
//    std::string dna_q = "GTAAGTACAACATCGATTAGTCACT";

    std::string dna_d = "";
    std::string dna_q = "";

//std::string dna_d = "TGGTGGGTGGCAGCGAGCGACCGTTTCTTGACGTCCTTGAGCTCTCTGCCCCTACTTGTCCTTCTTCCTCCTCCTCTGCTGCTTCTTGTCCTGTTCGCTATATGCATGCCTGCATTTGGTAGCCGCCAGTGATGGCTTGCTTTCGGAGGCAACCGCTCACCTTCCTCCTCTCCGCTCCTCCTCTTGGAGGCCGGTGTGCTGCTGCTGTTGCTTGTAGGAGGAGGGAGGGGAGGGGAGAGGGGTTGCAAGTT";
//    std::string dna_q = "TGGTGGGTGGCAGCGAGCGACCGTTTCTTGACGTCCTTGAGCTCTCTGCCCCTACTTgtccttcttcctcctcctctgCTGCTTCTTGTCCTGTTCGCTATATGCATGCCTGCATTTGGTAGCCGCCAGTGATGGCTTGCTTTcggaggcaaccgctcaccttccTCCTCTCCGCTCCTCCTCTTGGAGGCCGGTGTGCTGCTGCTGTTGCTtgtaggaggagggaggggaggggagaggggtTGCA";


    std::ofstream ofile;
    ofile.open("/media/bs674/ppi8t/testWAF/maizesk/inversionALignment");

    for (auto & c: dna_d) c = toupper(c);
    for (auto & c: dna_q) c = toupper(c);

    std::string _alignment_q;
    std::string _alignment_d;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -2;
    int32_t openGapPenalty1 = -3;
    int32_t extendGapPenalty1 = -1;
    int32_t invertyPenalty = -2;


    int32_t wfaSize = 10000;
    int32_t wfaSize2 = 50000;
    int32_t slidingWindowSize = 30000;
    int expectedCopies = 1;
    double maximumSimilarity = 0.6;
    int32_t min_wavefront_length = 20;
    int32_t max_distance_threshold = 100;

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 30;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    Scorei m(matchingScore, mismatchingPenalty);

    affine_penalties_t affine_penalties = {
            .match = 0,
            .mismatch = -mismatchingPenalty+matchingScore,
            .gap_opening = -openGapPenalty1+matchingScore,
            .gap_extension = -extendGapPenalty1+matchingScore,
    };

    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_512M);
    int32_t score = alignSlidingWindow(  dna_q,  dna_d, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );


    int32_t dindex = 0;
    int32_t qindex = 0;
    int32_t dindexCurrent = 0;
    int32_t qindexCurrent = 0;
    int32_t currentScore = 0;
    int32_t maximumScore = 0;
    int32_t maxi = 0;
    bool ifLastDeletion = false;
    bool ifLastInsertion = false;
    for ( int32_t i=0; i<_alignment_d.size(); ++i ){
        if ( _alignment_d[i] == _alignment_q[i] && _alignment_q[i] != '-' ){
            currentScore = currentScore + matchingScore;
            ifLastInsertion = false;
            ifLastDeletion = false;
        } else if (  _alignment_d[i] == '-' ) {
            if( ! ifLastInsertion ){
                currentScore = currentScore + openGapPenalty1;
            }
            currentScore = currentScore + extendGapPenalty1;
            ifLastInsertion = true;
            ifLastDeletion = false;
        } else if (  _alignment_q[i] == '-' ) {
            if( ! ifLastDeletion ){
                currentScore = currentScore + openGapPenalty1;
            }
            currentScore = currentScore + extendGapPenalty1;
            ifLastInsertion = false;
            ifLastDeletion = true;
        } else {
            ifLastInsertion = false;
            ifLastDeletion = false;
            currentScore = currentScore + mismatchingPenalty;
        }
        if (  _alignment_q[i] != '-' ) {
            qindexCurrent =  qindexCurrent + 1;
        }
        if (  _alignment_d[i] != '-' ) {
            dindexCurrent =  dindexCurrent + 1;
        }
        if( maximumScore < currentScore ){
            maximumScore = currentScore;
            dindex = dindexCurrent;
            qindex = qindexCurrent;
            maxi = i;
        }
        if( currentScore <= 0 ){
            i = _alignment_d.size();
        }
    }
    ++maxi;

    ofile << dna_d.substr(0, dindex) << std::endl << dna_q.substr(0, qindex) << std::endl;
    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;
    ofile << dindex << "\t" << qindex << std::endl;
    ofile << _alignment_d.substr(0, maxi) << std::endl << _alignment_q.substr(0, maxi) << std::endl;
    ofile << "above is the first alignment" << std::endl;


    dna_q = dna_q.substr(qindex, dna_q.size() - qindex);
    dna_d = dna_d.substr(dindex, dna_d.size() - dindex);

    dna_d = "";
    dna_q = "";

    std::string dna_d_v = getReverseComplementary(dna_d);
    std::string dna_q_v = getReverseComplementary(dna_q);
    score = alignSlidingWindow(  dna_q_v,dna_d_v, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );


    dindex = 0;
    qindex = 0;
    dindexCurrent = 0;
    qindexCurrent = 0;
    currentScore = 0;
    maximumScore = 0;
    ifLastDeletion = false;
    ifLastInsertion = false;
    for ( int32_t i=0; i<_alignment_d.size(); ++i ){
        if ( _alignment_d[i] == _alignment_q[i] && _alignment_q[i] != '-' ){
            currentScore = currentScore + matchingScore;
            ifLastInsertion = false;
            ifLastDeletion = false;
        } else if (  _alignment_d[i] == '-' ) {
            if( ! ifLastInsertion ){
                currentScore = currentScore + openGapPenalty1;
            }
            currentScore = currentScore + extendGapPenalty1;
            ifLastInsertion = true;
            ifLastDeletion = false;
        } else if (  _alignment_q[i] == '-' ) {
            if( ! ifLastDeletion ){
                currentScore = currentScore + openGapPenalty1;
            }
            currentScore = currentScore + extendGapPenalty1;
            ifLastInsertion = false;
            ifLastDeletion = true;
        } else {
            ifLastInsertion = false;
            ifLastDeletion = false;
            currentScore = currentScore + mismatchingPenalty;
        }
        if (  _alignment_q[i] != '-' ) {
            qindexCurrent =  qindexCurrent + 1;
        }
        if (  _alignment_d[i] != '-' ) {
            dindexCurrent =  dindexCurrent + 1;
        }
        if( maximumScore < currentScore ){
            maximumScore = currentScore;
            dindex = dindexCurrent;
            qindex = qindexCurrent;
            maxi = i;
        }
        if( currentScore <= 0 ){
            i = _alignment_d.size();
        }
    }
    ++maxi;

    ofile << dna_d_v.substr(0, dindex) << std::endl << dna_q_v.substr(0, qindex) << std::endl;
    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;
    ofile << dindex << "\t" << qindex << std::endl;
    ofile << _alignment_d.substr(0, maxi) << std::endl << _alignment_q.substr(0, maxi) << std::endl;

    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;
    ofile << "above is the second alignment" << std::endl;


    dna_q_v = dna_q_v.substr(qindex, dna_q_v.size() - qindex);
    dna_d_v = dna_d_v.substr(dindex, dna_d_v.size() - dindex);
    dna_q_v = "";
    dna_d_v = "";

    dna_d = getReverseComplementary(dna_d_v);
    dna_q = getReverseComplementary(dna_q_v);

    dna_d_v = getReverseComplementary(dna_d);


    score = alignSlidingWindow(  dna_q,dna_d_v, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    ofile << dna_d_v << std::endl << dna_q << std::endl;
    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;


    dna_d_v = "";
    dna_q = "";
    dna_d = getReverseComplementary(dna_d_v);

    score = alignSlidingWindow(  dna_q,dna_d_v, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    ofile << dna_d_v << std::endl << dna_q << std::endl;
    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;


    score = alignSlidingWindow(  dna_q,dna_d, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    ofile << dna_d << std::endl << dna_q << std::endl;
    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;


    int32_t maxDnaDIndex = 0;
    int32_t maxDnaDIndexe = 0;
    inversionAlignment( dna_q, dna_d, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, invertyPenalty, maxDnaDIndex, maxDnaDIndexe);
    ofile << "maxDnaDIndex:\t" << maxDnaDIndex << std::endl;
    ofile << "maxDnaQIndexe:\t" << maxDnaDIndexe << std::endl;
    int32_t length1 = dna_d.length();

    ofile << dna_d.substr (0, maxDnaDIndex) << "\t" << getReverseComplementary(dna_d.substr (maxDnaDIndex, maxDnaDIndexe-maxDnaDIndex)) << "\t" <<  dna_d.substr(maxDnaDIndexe, length1-maxDnaDIndexe) << std::endl;

    std::string _dna_d2 = dna_d.substr (0, maxDnaDIndex) + getReverseComplementary(dna_d.substr (maxDnaDIndex, maxDnaDIndexe-maxDnaDIndex)) +  dna_d.substr(maxDnaDIndexe, length1-maxDnaDIndexe);
    if( maxDnaDIndexe == length1 ){
        _dna_d2 = dna_d.substr (0, maxDnaDIndex) + getReverseComplementary(dna_d.substr (maxDnaDIndex, maxDnaDIndexe-maxDnaDIndex));
    }
    ofile << dna_d << std::endl;
    double inversion_PENALTY = -1;
    double MIN_ALIGNMENT_SCORE = 3;
    bool considerInversion = false;
    ofile << _dna_d2 << std::endl << dna_q << std::endl << std::endl;

    score = alignSlidingWindow(  dna_q,  _dna_d2, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    ofile << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;

    mm_allocator_delete(mm_allocator);
    ofile.close();
    ASSERT_EQ(0, 0);
}

