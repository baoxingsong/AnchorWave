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

TEST(alignSlidingWindow, c1){ // just to make sure that every line has been analysed

    std::string dna_d = "";
    std::string dna_q = "";

    std::string _alignment_q; std::string _alignment_d;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -2;
    int32_t openGapPenalty1 = -3;
    int32_t extendGapPenalty1 = -1;

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
    int32_t score = alignSlidingWindow(  dna_q,  dna_d, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, slidingWindowSize, wfaSize,  matchingScore,  mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m );
    std::cout << _alignment_d << std::endl << _alignment_q << std::endl << score << std::endl;
    mm_allocator_delete(mm_allocator);
    ASSERT_EQ(0, 0);
}

TEST(minimap2, c1){ // just to make sure that every line has been analysed
    std::string querySeq = "";
    std::string refSeq = "";

    int w = 1;
    int k = 11;
    int is_hpc = 0; // no, do not use  homopolymer-compressed (HPC) minimizers.
    int bucket_bits = 2;
    int n = 1;

    int32_t matchingScore = 1;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -5;
    int32_t extendGapPenalty1 = -1;

    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
//    mopt.flag |= MM_F_CIGAR; // perform alignment
//    mopt.flag &= ~ MM_F_CIGAR; // DO NOT perform alignment
//    mopt.best_n = 1;
    mopt.flag |= MM_F_CIGAR; // DO NOT perform alignment
    mopt.flag |= MM_F_NO_PRINT_2ND;

//    mopt.max_gap_ref = 5000;
//    mopt.max_gap = 5000;
    mopt.bw = 400;
    mopt.flag |= MM_F_NO_LJOIN;
//    mopt.chain_gap_scale=150;
    mopt.a = matchingScore;
    mopt.b = -mismatchingPenalty;
    mopt.q = -openGapPenalty1;
    mopt.e = -extendGapPenalty1;
    mopt.mid_occ = 2; // ignore seeds with occurrences above this threshold
    mopt.min_cnt = 2;// min number of minimizers on each chain
    int32_t referenceSeqLength = refSeq.length();
    int32_t querySeqLength = querySeq.length();
    std::cout << referenceSeqLength << "\t" << querySeqLength << std::endl;
    char * reference_seq_array = new char [referenceSeqLength+1];

    //char query_seq_array[querySeqLength + 1];
    char * query_seq_array = new char[querySeqLength + 1];
    strcpy(reference_seq_array, refSeq.c_str());
    const char *refseq[1] = {reference_seq_array};
    strcpy(query_seq_array, querySeq.c_str());

    mm_idx_t *mi = mm_idx_str(w, k, is_hpc, bucket_bits, n, refseq, 0);

    mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    mm_reg1_t *reg;
    int j, n_reg;
    reg = mm_map(mi, querySeqLength, query_seq_array, &n_reg, tbuf, &mopt, 0); // get all hits for the query
    int i;
    if (n_reg > 0  && (&reg[0])->rev == 0) {
        mm_reg1_t *r = &reg[0];
        printf("\t%d\t%d\t",  r->rs, r->re);
        printf("\t%d\t%d\t%c\t",  r->qs, r->qe, "+-"[r->rev]);
    }
    std::cout  << std::endl << "line 153" << std::endl << std::endl;
    for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
        mm_reg1_t *r = &reg[j];
        printf("\t%d\t%d\t",  r->rs, r->re);
        printf("\t%d\t%d\t%c\t",  r->qs, r->qe, "+-"[r->rev]);

        for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
            printf("%d%c", r->p->cigar[i]>>4, "MIDNSH"[r->p->cigar[i]&0xf]);
//        free(r->p);
        std::cout << std::endl << std::endl;
    }
    free(reg);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
    std::cout << "line 744" << std::endl;
    delete(reference_seq_array);
    delete(query_seq_array);
    ASSERT_EQ(0, 0);
}