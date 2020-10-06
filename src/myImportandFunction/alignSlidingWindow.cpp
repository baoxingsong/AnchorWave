//
// Created by song on 8/5/18.
//

#include "alignSlidingWindow.h"


#include "../../WFA/gap_affine/affine_wavefront_align.h"
#include <stdlib.h>


int64_t alignSlidingWindow( const std::string& dna_q, const std::string& dna_d,
                         std::string & _alignment_q, std::string & _alignment_d, const int & slidingWindowSize, const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
                            const int32_t & extendGapPenalty1 ){
    int64_t totalScore = 0;

    _alignment_q = "";
    _alignment_d = "";

    size_t _length_of_q=dna_q.size();
    size_t _length_of_d=dna_d.size();

//    std::cout << dna_d << std::endl;
//    std::cout << dna_q << std::endl;
    //2^15 = 32768
    //of the maximum length of the windowSize of is about 32000/2 = 16000
    size_t databaseStart=1;
    size_t databaseEnd = 0;
    size_t queryStart=1;
    size_t queryEnd = 0;

    //mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_32G);
    // Set penalties
    affine_penalties_t affine_penalties = {
            .match = 0,
            .mismatch = -mismatchingPenalty-matchingScore,
            .gap_opening = -openGapPenalty1-matchingScore,
            .gap_extension = -extendGapPenalty1-matchingScore,
    };
    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_32M);
    try {
        while (databaseStart < _length_of_d && queryStart < _length_of_q) {
            databaseEnd = databaseStart + slidingWindowSize;
            queryEnd = queryStart + slidingWindowSize;

            if ( (databaseEnd >= _length_of_d || queryEnd >= _length_of_q) && (_length_of_d-databaseEnd)<=slidingWindowSize && (_length_of_q-queryEnd)<=slidingWindowSize ) {
                databaseEnd = _length_of_d;
                queryEnd = _length_of_q;
            }else{
                if (databaseEnd > _length_of_d) {
                    databaseEnd = _length_of_d;
                }
                if (queryEnd > _length_of_q) {
                    queryEnd = _length_of_q;
                }
            }

            std::string qSeq = getSubsequence(dna_q, queryStart, queryEnd);
            std::string dSeq = getSubsequence(dna_d, databaseStart, databaseEnd);

            const char *pattern = dSeq.c_str();
            const char *text = qSeq.c_str();
            // Init Affine-WFA
            affine_wavefronts_t *affine_wavefronts = affine_wavefronts_new_complete(
                    strlen(pattern), strlen(text), &affine_penalties, NULL, mm_allocator);
            //std::cout << "line 766" << std::endl;
            // Align
            affine_wavefronts_align(
                    affine_wavefronts, pattern, strlen(pattern), text, strlen(text));
            //std::cout << "line 770" << std::endl;
            const int score = edit_cigar_score_gap_affine(
                    &affine_wavefronts->edit_cigar,&affine_penalties);
            totalScore += score;
            edit_cigar_t *const edit_cigar = &(affine_wavefronts->edit_cigar);
            //            std::cout << "line 762" << std::endl;
            char *const operations = edit_cigar->operations;

            int i, pattern_pos = 0, text_pos = 0;
            for (i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i) {
                switch (operations[i]) {
                    case 'M':
                        _alignment_q += qSeq[text_pos];
                        _alignment_d += dSeq[pattern_pos];
                        ++queryStart;
                        ++databaseStart;
                        pattern_pos++;
                        text_pos++;
                        break;
                    case 'X':
                        _alignment_q += qSeq[text_pos];
                        _alignment_d += dSeq[pattern_pos];
                        ++queryStart;
                        ++databaseStart;
                        pattern_pos++;
                        text_pos++;
                        break;
                    case 'I':
                        _alignment_q += qSeq[text_pos];
                        _alignment_d += '-';
                        text_pos++;
                        ++queryStart;
                        break;
                    case 'D':
                        _alignment_q += '-';
                        _alignment_d += dSeq[pattern_pos];
                        pattern_pos++;
                        ++databaseStart;
                        break;
                    default:
                        break;
                }
            }
            affine_wavefronts_delete(affine_wavefronts);

        }
    } catch(std::bad_alloc& ex) {
        mm_allocator_delete(mm_allocator);
        alignSlidingWindow( dna_q, dna_d, _alignment_q,  _alignment_d, slidingWindowSize/2 );
        return 0;
    }catch (std::exception const& e){
        mm_allocator_delete(mm_allocator);
        alignSlidingWindow( dna_q, dna_d, _alignment_q,  _alignment_d, slidingWindowSize/2 );
        return 0;
    }
    catch (...) {
        std::cout << "We caught an exception of an undetermined type\n";
        return 0;
    }
    while( databaseStart<=_length_of_d ){
        _alignment_q += '-';
        _alignment_d += dna_d[databaseStart-1];
        ++databaseStart;
    }
    while( queryStart<=_length_of_q ){
        _alignment_q += dna_q[queryStart-1];
        _alignment_d += '-';
        ++queryStart;
    }
    //std::cout << _alignment_d << std::endl;
    //std::cout << _alignment_q << std::endl << std::endl << std::endl;
    mm_allocator_delete(mm_allocator);
    return totalScore;
}




