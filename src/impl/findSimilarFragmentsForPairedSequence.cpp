
/*
 * =====================================================================================
 *
 *       Filename:  findSimilarFragmentsForPairedSequence.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/25/2018 01:13:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

/**
 * all the pair sequence alignment operations
 *
 * */
#include "findSimilarFragmentsForPairedSequence.h"


// the position of an alignment match, for unique aim using set
class StartAlignment{
    private:
        int32_t refPosition;
        int32_t queryPosition;
    public:
        StartAlignment(int32_t & _refPosition, int32_t & _queryPosition){
            refPosition = _refPosition;
            queryPosition = _queryPosition;
        }
        int32_t getRefPosition() const {
            return refPosition;
        }
        void setRefPosition(int32_t & _refPosition) {
            refPosition = _refPosition;
        }
        int32_t getQueryPosition() const {
            return queryPosition;
        }
        void setQueryPosition(int32_t & _queryPosition) {
            queryPosition = _queryPosition;
        }
        bool operator< (const StartAlignment & c) const {
            if ( this->refPosition < c.getRefPosition() ){
                return true;
            }else if( this->refPosition == c.getRefPosition() ){
                return this->queryPosition < c.getQueryPosition();
            }
            return false;
        }
};



//set up seeds using smither-waterman method
std::set<Seed> getSeeds( int8_t * seq1, int8_t * seq2, int8_t * seq1_rev_com, int8_t * seq2_rev_com,
                           const int32_t & length1, const int32_t & length2, const int32_t & miniReferenceSize, const int32_t & _open_gap_penalty1,
                           const int32_t & _extend_gap_penalty1,
                           const int32_t & startPosition2, const int32_t & mini_cns_score,
                         const int32_t & windowsSize, const Scorei & m, const int32_t & step_size ){

    std::set<Seed> seeds;
    int32_t maxScore;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;


    bool returePosition = false;
    SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                  maxScore, endPosition1, endPosition2, m, returePosition);
    if (mini_cns_score <= maxScore) {
        if( !returePosition ){
            returePosition=true;
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                          maxScore, endPosition1, endPosition2, m, returePosition);
        }
//        std::cout << "line 182 end1:" << endPosition1 << " start2:" << startPosition2 << " end2:" << endPosition2 << std::endl;
        int32_t e1 = 1 + endPosition1< this_windowsSize+windowsSize? 1 + endPosition1 : this_windowsSize+windowsSize;
        int32_t e2 = 1 + endPosition2< this_windowsSize+windowsSize? 1 + endPosition2 : this_windowsSize+windowsSize;
        //int32_t oldMaxScore = maxScore;
        SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1),
                      seq2_rev_com + (length2 - (startPosition2 + endPosition2) - 1),
                      e1, e2, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                      endPosition1_rc, endPosition2_rc, m, returePosition);// this is the reverse alignment



        int32_t end1 = endPosition1; // 0based
        int32_t end2 = startPosition2+endPosition2; //0-based
        int32_t start1 = endPosition1-endPosition1_rc; //0-based coordinate
        int32_t start2 = startPosition2+endPosition2-endPosition2_rc; //0-based coordinate
        Seed seed0(start1, end1, start2, end2);
        seeds.insert(seed0);
        if( start1 >= miniReferenceSize ) {
            std::set<Seed> seeds1 = getSeeds(seq1, seq2, seq1_rev_com + (length1 - start1), seq2_rev_com, start1,
                                             length2, miniReferenceSize, _open_gap_penalty1, _extend_gap_penalty1,
                                             startPosition2, mini_cns_score, windowsSize, m, step_size);
            if (!seeds1.empty()) {
                for (Seed seed1 : seeds1) {
                    seeds.insert(seed1);
                }
            }
        }
        if( length1-end1-1 >= miniReferenceSize ) {
            std::set<Seed> seeds2 = getSeeds( seq1+end1+1, seq2, seq1_rev_com, seq2_rev_com, length1-end1-1,length2, miniReferenceSize,
                    _open_gap_penalty1, _extend_gap_penalty1, startPosition2, mini_cns_score, windowsSize, m, step_size );
            if ( !seeds2.empty() ){
                for( Seed seed : seeds2 ){
                    int32_t start11 = seed.getStart1()+end1+1;
                    int32_t start21 = seed.getStart2();
                    int32_t end11 = seed.getEnd1()+end1+1;
                    int32_t end21 = seed.getEnd2();
                    Seed seed2(start11, end11, start21, end21);
                    seeds.insert(seed2);
                }
            }
        }
    }
    return seeds;
}

//same with the above functions just give different parameters, and call the above function
void getSeeds( int8_t * seq1, int8_t * seq2, int8_t * seq1_rev_com, int8_t * seq2_rev_com, const int32_t & length1,
        const int32_t & length2, const int32_t & miniReferenceSize, const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
        const int32_t & matrix_boundary_distance, const int32_t & startPosition2, const int32_t & mini_cns_score,
        const int32_t & windowsSize, int32_t & this_windowsSize_return, const Scorei & m,
        const int32_t & step_size, std::set<Seed> & seeds){

    int32_t maxScore;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t this_windowsSize = (startPosition2 + windowsSize) > length2 ? (length2 - startPosition2) : windowsSize;
    bool returePosition = false;
    SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                  maxScore, endPosition1, endPosition2, m, returePosition);
    if (mini_cns_score <= maxScore) {
        if( !returePosition ){
            returePosition=true;
            SmithWaterman(seq1, seq2 + startPosition2, length1, this_windowsSize, _open_gap_penalty1, _extend_gap_penalty1,
                          maxScore, endPosition1, endPosition2, m, returePosition);
        }
        int32_t e1 = 1 + endPosition1< this_windowsSize+windowsSize? 1 + endPosition1 : this_windowsSize+windowsSize;
        int32_t e2 = 1 + endPosition2< this_windowsSize+windowsSize? 1 + endPosition2 : this_windowsSize+windowsSize;
        //int32_t oldMaxScore = maxScore;
        // it should be possible to avoid run this round of smith waterman algorithm, but maybe would not save a lot of times
        SmithWaterman(seq1_rev_com + (length1 - 1 - endPosition1),
                      seq2_rev_com + (length2 - (startPosition2 + endPosition2) - 1),
                      e1, e2, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                      endPosition1_rc, endPosition2_rc, m, returePosition);// this is the reverse alignment

        //assert(oldMaxScore == maxScore);
        //if ((endPosition1_rc ) < mini_cns_size && (endPosition2_rc ) < mini_cns_size) {
        int32_t end1 = endPosition1; // 0based
        int32_t end2 = startPosition2+endPosition2; //0-based
        int32_t start1 = endPosition1-endPosition1_rc; //0-based coordinate
        int32_t start2 = startPosition2+endPosition2-endPosition2_rc; //0-based coordinate
//        std::cout << "line 202 start1:" << start1 << " end1:" << end1 << " start2:" << start2 << " end2:" << end2 << " maxScore:" << maxScore << std::endl;
        Seed seed0(start1, end1, start2 , end2);

        if( seeds.find(seed0) == seeds.end()){
            seeds.insert(seed0);
            if( start1>= miniReferenceSize){
                std::set<Seed> seeds1 = getSeeds( seq1, seq2, seq1_rev_com+(length1-start1), seq2_rev_com, start1,
                        length2, miniReferenceSize, _open_gap_penalty1, _extend_gap_penalty1,
                        startPosition2, mini_cns_score, windowsSize, m, step_size );
                if ( !seeds1.empty() ){
                    for( Seed seed1 : seeds1 ){
                        seeds.insert(seed1);
                    }
                }
            }
            if( length1-end1-1 >= miniReferenceSize ) {
                std::set<Seed> seeds2 = getSeeds(seq1 + end1 + 1, seq2, seq1_rev_com, seq2_rev_com, length1 - end1 - 1,
                                                 length2, miniReferenceSize, _open_gap_penalty1, _extend_gap_penalty1,
                                                 startPosition2, mini_cns_score, windowsSize, m, step_size);
                if (!seeds2.empty()) {
                    for (Seed seed : seeds2) {
                        int32_t start11 = seed.getStart1() + end1 + 1;
                        int32_t start21 = seed.getStart2();
                        int32_t end11 = seed.getEnd1() + end1 + 1;
                        int32_t end21 = seed.getEnd2();
                        Seed seed2(start11, end11, start21, end21);
                        seeds.insert(seed2);
                    }
                }
            }
        }
        //}
    }
    this_windowsSize_return = this_windowsSize;
}


// using x-drop algorithm to extend the seeds alignment
void x_extend_seed(int32_t & start1, int32_t & end1, int32_t & start2, int32_t & end2,
                   int32_t & maxScore, int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com,
                   const double & length1d, const double & length2d,const int32_t & length1,
                   const int32_t & length2,const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                   const int32_t & xdrop, const int32_t & w,
                   const double & kValue, const double & lambda, double & eValue, double & pvalue, const Scorei & m ){

    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t s1 = start1;
    int32_t s2 = start2;
    // TODO this SemiGlobal_xextend implementation could be faster
    SemiGlobal_xextend(seq1+start1, seq2+start2, length1-start1, length2-start2, _open_gap_penalty1, _extend_gap_penalty1,
    maxScore, endPosition1, endPosition2, m, xdrop, w);
    if( endPosition1<0 || endPosition2 <0 ){
        pvalue = 100;
        return;
    }

    SemiGlobal_xextend(seq1_rev_com + (length1 - 1 - (start1+endPosition1)),
                  seq2_rev_com + (length2 - 1 - (start2 + endPosition2)),
                  start1+1 + endPosition1, start2+1+endPosition2, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                  endPosition1_rc, endPosition2_rc, m,  xdrop, w);// this is the reverse alignment
    if( endPosition1_rc<0 || endPosition2_rc <0 ){
        pvalue = 100;
        return;
    }
    end1 = s1 + endPosition1; // 0based
    end2 = s2 + endPosition2; //0-based
    start1 = s1+endPosition1-endPosition1_rc; //0-based coordinate
    start2 = s2+endPosition2-endPosition2_rc; //0-based coordinate
    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
    pvalue = 1.0 - exp(-1.0*eValue);
}


// same with the above one, but return cigar
std::vector<uint32_t > x_extend_seed(int32_t & start1, int32_t & end1, int32_t & start2, int32_t & end2,
                   int32_t & maxScore, int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com,
                   const double & length1d, const double & length2d, const int32_t & length1,
                   const int32_t & length2,const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                   const int32_t & xdrop, const int32_t & w,
                   const double & kValue, const double & lambda, double & eValue, double & pvalue, const Scorei & m, class Matrix & T){

    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t s1 = start1;
    int32_t s2 = start2;
    SemiGlobal_xextend(seq1+start1, seq2+start2, length1-start1, length2-start2, _open_gap_penalty1, _extend_gap_penalty1,
                       maxScore, endPosition1, endPosition2, m, xdrop, w);

    if( endPosition1<0 || endPosition2 <0 ){
        pvalue = 100;
        std::vector<uint32_t > cigar;
        return cigar;
    }
    int32_t iii, jjj;
    std::vector<uint32_t > cigar = SemiGlobal_xextend(seq1_rev_com + (length1 - 1 - (start1+endPosition1)),
                       seq2_rev_com + (length2 - 1 - (start2 + endPosition2)),
                       start1+1 + endPosition1, start2+1+endPosition2, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                       endPosition1_rc, endPosition2_rc, m,  xdrop, w, T, iii, jjj);// this is the reverse alignment
    if( endPosition1_rc<0 || endPosition2_rc <0 ){
        pvalue = 100;
        return cigar;
    }
    end1 = s1 + endPosition1; // 0based
    end2 = s2 + endPosition2; //0-based
    start1 = s1+endPosition1-endPosition1_rc; //0-based coordinate
    start2 = s2+endPosition2-endPosition2_rc; //0-based coordinate
    end1 -= iii;
    end2 -= jjj;
    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
    pvalue = 1.0 - exp(-1.0*eValue);
    return cigar;
}

// merge seeds into larger ones
void link_seeds( std::set<Seed> & seeds ){
    std::vector<Seed> old_seeds;
    for( Seed s : seeds ){
        old_seeds.push_back(s);
    }
    std::sort(old_seeds.begin(), old_seeds.end(), [](Seed a, Seed b) {
        return a.getStart1() < b.getStart1();
    });

    bool everChange = true;
    int i, j;
    while(everChange){
        std::vector<Seed> new_seeds;
        everChange = false;
        std::set<int> usedList;
        for( i=0; i<old_seeds.size(); ++i ){
            if( usedList.find(i) == usedList.end() ){
                Seed s1 = old_seeds[i];
                bool ever_s1_used = false;
                for( j=i+1; j<old_seeds.size(); ++j ){
                    if( usedList.find(j) == usedList.end() ) {
                        Seed s2 = old_seeds[j];
                        if( s1.getStart2() <= s2.getStart2() && s1.getEnd1()>=s2.getEnd1() && s1.getEnd2()>=s2.getEnd2()){ // s2 is a part of s1
                            usedList.insert(j);
                        }else if ( s1.getStart2() <= s2.getStart2() && s1.getEnd1()<=s2.getEnd1() && s1.getEnd2()<=s2.getEnd2()  ) {
                            int a = s2.getStart1() - s1.getStart1();
                            if ( (a==s2.getStart2()-s1.getStart2()) && ( a==s2.getEnd1()-s1.getEnd1() ) && (a==s2.getEnd2()-s1.getEnd2()) ){ // merger s1 and s2
                                int32_t start1 = s1.getStart1();
                                int32_t start2 = s1.getStart2();
                                int32_t end1 = s2.getEnd1();
                                int32_t end2 = s2.getEnd2();
                                Seed s(start1, end1, start2, end2);
                                new_seeds.push_back(s);
                                usedList.insert(i);
                                usedList.insert(j);
                                ever_s1_used = true;
                                j = old_seeds.size(); //stop here
                            }
                        }else if ( s2.getStart1() > s1.getEnd1() || s2.getStart2() > s1.getEnd2() ) {
                            break;
                        }
                    }
                }
                if(  ever_s1_used ){
                    everChange = true;
                }else{
                    new_seeds.push_back(s1);
                }
            }
        }
        old_seeds = new_seeds;
    }
    seeds.clear();
    for( Seed s : old_seeds ){
        seeds.insert(s);
    }
}


// using smith-waterman approach to find seeds and using x-drop to extend the seeds
void getAllExtendSeed( int8_t * seq1, int8_t * seq1_rev_com,
                       int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                       const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                       const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & matchingScore,
                       const int32_t & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                       std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                       const double & lambda, const double & kValue, const int32_t & zDrop,
                       const int32_t & w, const int32_t & xDrop,  std::vector<Seed> & x_extend_seeds ){


    int32_t miniReferenceSize = mini_cns_score/matchingScore;
    int32_t startPosition2 = 0;
    int32_t this_windowsSize_return;
    int32_t maxScore;


    int32_t start1;
    int32_t end1;
    int32_t start2 ;
    int32_t end2;

    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    std::set<Seed> seeds;

    bool notEnd = true;
    while( notEnd ) { // could not put == here, since the first one is startPosition2 + 0
        getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2, miniReferenceSize,
                 _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                 startPosition2, mini_cns_score, windowsSize,
                 this_windowsSize_return, m, step_size, seeds);
//        startPosition2 += (this_windowsSize_return - windowsSize  + step_size);
        startPosition2 += step_size;
        if( (startPosition2 + windowsSize) > length2 ){
            getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2, miniReferenceSize,
                     _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                     startPosition2, mini_cns_score, windowsSize,
                     this_windowsSize_return, m, step_size, seeds);
            notEnd = false;
            break;
        }
    }
    link_seeds( seeds );  // this function is very fast, by using less seeds, the program is significantly faster

    for (Seed seed : seeds) {
        start1 = seed.getStart1()+1;
        end1 = seed.getEnd1()+1;
        start2 = seed.getStart2()+1;
        end2 = seed.getEnd2()+1;

        x_extend_seed(start1, end1, start2, end2, maxScore, seq1, seq1_rev_com, seq2, seq2_rev_com,
                      length1d, length2d, length1, length2, _open_gap_penalty1, _extend_gap_penalty1, xDrop, w,
                      kValue, lambda, eValue, pvalue, m );

//        std::cout << "line 439:" <<  start1 << " " << end1 << " " << start2 << " " << end2 << std::endl;
        if (pvalue < pvalues ) {
            bool never_used = true;
            for ( int32_t i=0; i<x_extend_seeds.size(); ++i ){
                if( x_extend_seeds[i].getStart1() == start1 && x_extend_seeds[i].getStart2() == start2 ){
                    never_used = false;
                    if( end1 > x_extend_seeds[i].getEnd1() && end2 > x_extend_seeds[i].getEnd2() ){
                        x_extend_seeds[i].setEnd1(end1);
                        x_extend_seeds[i].setEnd2(end2);
                    }
                }else if( x_extend_seeds[i].getEnd1() == end1 && x_extend_seeds[i].getEnd2() == end2 ){
                    if( start1 < x_extend_seeds[i].getStart1() && start2 < x_extend_seeds[i].getStart2() ){
                        x_extend_seeds[i].setStart1(start1);
                        x_extend_seeds[i].setStart2(start2);
                    }
                    never_used = false;
                }
            }
            if( never_used ){
                Seed x_seed(start1, end1, start2, end2);
                x_extend_seeds.push_back(x_seed);
            }
        }
    }
}

// the first alignment algorithm
// using smith-waterman approach to find seeds and using x-drop to extend the seeds, return result for output
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                   const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                   const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,const int32_t & matchingScore,
                   const int32_t & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                   std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                   const double & lambda, const double & kValue, const int32_t & w, const int32_t & xDrop){
    int32_t miniReferenceSize = mini_cns_score/matchingScore;

    std::vector<PairedSimilarFragment> pairedSimilarFragments;

    std::set<StartAlignment> startAlignments;

    int32_t startPosition2 = 0;
    int32_t this_windowsSize_return;
    int32_t maxScore;

    int32_t start1;
    int32_t end1;
    int32_t start2 ;
    int32_t end2;

    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    std::set<Seed> seeds;
    bool notEnd = true;
    while( notEnd ) {
        getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2, miniReferenceSize,
                 _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                 startPosition2, mini_cns_score, windowsSize,
                 this_windowsSize_return, m, step_size, seeds);
        startPosition2 += step_size;
        if( startPosition2 > length2 ){
            notEnd = false;
            break;
        }
        if( (startPosition2 + windowsSize) > length2 ){
            getSeeds(seq1, seq2, seq1_rev_com, seq2_rev_com, length1, length2, miniReferenceSize,
                     _open_gap_penalty1, _extend_gap_penalty1, matrix_boundary_distance,
                     startPosition2, mini_cns_score, windowsSize,
                     this_windowsSize_return, m, step_size, seeds);
            notEnd = false;
            break;
        }
    }
    std::cout << "seeds size" << seeds.size() << std::endl;
//    std::cout << "line 461 begin to link seeds" << std::endl;
    link_seeds( seeds );  // this function is very fast, by using less seeds, the program is significantly faster
//    std::cout << "line 463 generating seeds done" << std::endl;

//    for (Seed seed : seeds) {
//
//    }

    class Matrix T (length1 + 1, length2 + 1);
    for (Seed seed : seeds) {
        start1 = seed.getStart1()+1;
        end1 = seed.getEnd1()+1;
        start2 = seed.getStart2()+1;
        end2 = seed.getEnd2()+1;
        std::vector<uint32_t > cigar = x_extend_seed(start1, end1, start2, end2, maxScore, seq1, seq1_rev_com, seq2, seq2_rev_com,
                                        length1d, length2d, length1, length2, _open_gap_penalty1, _extend_gap_penalty1,
                                                xDrop, w, kValue, lambda,  eValue,  pvalue,  m,  T);

        if (pvalue < pvalues ) {
            bool never_used = true;
            PairedSimilarFragment pairedSimilarFragment(start1+1, end1+1, start2+1, end2+1, maxScore, cigar, pvalue, eValue);
            for ( int32_t i=0; i<pairedSimilarFragments.size(); ++i ){
                if( pairedSimilarFragments[i].getStart1() == (start1+1) && pairedSimilarFragments[i].getStart2() == (start2+1) ){
                    never_used = false;
                    if( end1 >= pairedSimilarFragments[i].getEnd1() && end2 >= pairedSimilarFragments[i].getEnd2() ){
                        pairedSimilarFragments[i] = pairedSimilarFragment;
                    }
                    break;
                }else if( pairedSimilarFragments[i].getEnd1() == (end1+1) && pairedSimilarFragments[i].getEnd2() == (end2+1) ){
                    never_used = false;
                    if( start1 <= pairedSimilarFragments[i].getStart1() && start2 <= pairedSimilarFragments[i].getStart2() ){
                        pairedSimilarFragments[i] = pairedSimilarFragment;
                    }
                    break;
                }
            }
            if( never_used ){
                pairedSimilarFragments.push_back(pairedSimilarFragment);
            }
        }
    }
    return pairedSimilarFragments;
}




// the second algorithm, continue from the output of first result
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a 2-piece gap cost approach
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                                                                           const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2, const int32_t & matchingScore,
                                                                           const int32_t & mismatchingPenalty, const Scorei & m,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                                                                           std::vector<Seed> & x_extend_seeds, class Matrix & T){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t end1;
    int32_t start2;
    int32_t end2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();

        SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1),
                   seq2_rev_com + (length2 - endPosition2 - 1),
                   1 + endPosition1, endPosition2 + 1, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2,
                   _extend_gap_penalty2, maxScore, endPosition1_rc, endPosition2_rc, m, false, zDrop, bandwidth, T);// this is the reverse alignment
//        std::cout << "line 779 " << maxScore << std::endl;
        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based

        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {
            startAlignments.insert(startAlignment);
//            std::cout <<"line 786" << std::endl;
            std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1, seq2 + start2, length1 - start1,
                                                     length2 - start2, _open_gap_penalty1, _extend_gap_penalty1,
                                                     _open_gap_penalty2, _extend_gap_penalty2,
                                                     maxScore, end1, end2, m, true, zDrop, bandwidth, T);
            if( cigar.size() > 0  && end1>0 && end2>0) {
                end1 += start1; // 0based
                end2 += start2; // 0based

                int32_t toRemove = -1;
                bool found = false;
                for (int32_t i = 0; i < pairedSimilarFragments.size(); ++i) {
                    if (pairedSimilarFragments[i].getEnd1() == end1 + 1 &&
                        pairedSimilarFragments[i].getEnd2() == end2 + 1) {
                        found = true;
                        if (start1 + 1 < pairedSimilarFragments[i].getStart1() &&
                            start2 + 1 < pairedSimilarFragments[i].getStart2()) {
                            toRemove = i;
                        }
                    }
                }
                if (!found or toRemove >= 0) {

                    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                    pvalue = 1.0 - exp(-1.0 * eValue);

                    PairedSimilarFragment pairedSimilarFragment(start1 + 1, end1 + 1, start2 + 1, end2 + 1, maxScore,
                                                                cigar,
                                                                pvalue,
                                                                eValue); // start1， start2, end1 and end2 are 1 based
                    // here the ordinate put into pairedSimilarFragment is 1 based
                    if (toRemove >= 0) {
                        pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                    } else {
                        pairedSimilarFragments.push_back(pairedSimilarFragment);
                    }
                }
            }
        }
    }
    return pairedSimilarFragments;
}


// the second algorithm, continue from the output of first result
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a 2-piece gap cost approach
// maximumAlignLength is the maximum length of sequence being used for dynamic programming algorithm
// if the input sequence is very long, we will have out of RAM problem, so here I set maximumAlignLength parameter to only align a fragment of the long sequence
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                                                                           const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2, const int32_t & matchingScore,
                                                                           const int32_t & mismatchingPenalty, const Scorei & m,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                                                                           std::vector<Seed> & x_extend_seeds, class Matrix & T, int32_t & maximumAlignLength){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t end1;
    int32_t start2;
    int32_t end2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    for( Seed x_seed :  x_extend_seeds){
        bool useThisSeed = true;
        for ( PairedSimilarFragment pairedSimilarFragment : pairedSimilarFragments ){
            if( pairedSimilarFragment.getStart1() <= x_seed.getStart1() && pairedSimilarFragment.getEnd1() >= x_seed.getEnd1() &&
                    pairedSimilarFragment.getStart2() <= x_seed.getStart2() && pairedSimilarFragment.getEnd2() >= x_seed.getEnd2()){
                useThisSeed = false;
                break;
            }
        }
        if ( useThisSeed ){
            endPosition1 = x_seed.getEnd1();
            endPosition2 = x_seed.getEnd2();
            int32_t length11 = 1 + endPosition1 < maximumAlignLength ? 1 + endPosition1 : maximumAlignLength;
            int32_t length12 = 1 + endPosition2 < maximumAlignLength ? 1 + endPosition2 : maximumAlignLength;
            if( static_cast<int16_t>(*(seq1_rev_com + (length1 - 1 - endPosition1))) != static_cast<int16_t>(*(seq2_rev_com+(length2 - endPosition2 - 1))) ){
                //todo, this is a patch, it should never run, but it ever ran. This suggests the seed x-drop extand code has small bug
                --endPosition1;
                --endPosition2;
            }
            SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1),
                       seq2_rev_com + (length2 - endPosition2 - 1),
                       length11, length12, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2,
                       _extend_gap_penalty2, maxScore, endPosition1_rc, endPosition2_rc, m, false, zDrop, bandwidth, T);// this is the reverse alignment

//            for ( int32_t  seq_index = (length1 - 1 - endPosition1); seq_index <= (length1 - 1 - endPosition1 + 100); ++seq_index){
//                std::cout << static_cast<int16_t>(*(seq1_rev_com+seq_index)) ;
//            }
//            std::cout << std::endl;
//
//            for ( int32_t  seq_index = (length2 - 1 - endPosition2); seq_index <= (length2 - 1 - endPosition2 + 100); ++seq_index){
//                std::cout << static_cast<int16_t>(*(seq2_rev_com+seq_index)) ;
//            }
//            std::cout << std::endl;

//            std::cout << "line 655 " << maxScore << " endPosition1_rc:" << endPosition1_rc << " endPosition2_rc:" << endPosition2_rc << std::endl;




            start1 = endPosition1 - endPosition1_rc; //0 based
            start2 = endPosition2 - endPosition2_rc ; // 0 based

            StartAlignment startAlignment(start1, start2);
            if (startAlignments.find(startAlignment) == startAlignments.end()) {
                startAlignments.insert(startAlignment);
    //            std::cout <<"line 638" << std::endl;
                int32_t length21 = length1 - start1 < maximumAlignLength ? length1 - start1 : maximumAlignLength;
                int32_t length22 = length2 - start2 < maximumAlignLength ? length2 - start2 : maximumAlignLength;

                std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1, seq2 + start2, length21,
                                                         length22, _open_gap_penalty1, _extend_gap_penalty1,
                                                         _open_gap_penalty2, _extend_gap_penalty2,
                                                         maxScore, end1, end2, m, true, zDrop, bandwidth, T);
    //            std::cout << "line 670 " << maxScore << std::endl;
                if( cigar.size() > 0  && end1>0 && end2>0){
                    end1 += start1; // 0based
                    end2 += start2; // 0based

                    int32_t toRemove = -1;
                    bool found = false;
                    for ( int32_t i=0; i< pairedSimilarFragments.size(); ++i) {
                        if ( pairedSimilarFragments[i].getEnd1() == end1+1 && pairedSimilarFragments[i].getEnd2() == end2+1 ){
                            found = true;
                            if( start1+1 < pairedSimilarFragments[i].getStart1() && start2+1 < pairedSimilarFragments[i].getStart2() ){
                                toRemove = i;
                            }
                        }
                    }
                    if( !found or toRemove>=0 ){

                        eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                        pvalue = 1.0 - exp(-1.0*eValue);

                        PairedSimilarFragment pairedSimilarFragment(start1+1, end1 + 1, start2+1, end2 + 1, maxScore, cigar,
                                                                    pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                        // here the ordinate put into pairedSimilarFragment is 1 based
                        if( toRemove>=0  ){
                            pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                        }else{
                            pairedSimilarFragments.push_back(pairedSimilarFragment);
                        }
                    }
                }else{
                    int32_t s1 = length1 - 1 - endPosition1;
                    int32_t s2 = length2 - endPosition2 - 1;

                    std::cerr << "error at findSimilarFragmentsForPairedSequence line 700 s1:" << s1 << " s2:" << s2 << " start1:" << start1 << " start2:" << start2 << " ";
                    std::cerr << "x_seed.start1:" << x_seed.getStart1() << "x_seed.end1:" << x_seed.getEnd1() << "x_seed.start2:" << x_seed.getStart2() << "x_seed.end2:" << x_seed.getEnd1() << std::endl;
                }
            }
        }
    }

    return pairedSimilarFragments;
}



// the second algorithm
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a 2-piece gap cost approach
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2, const int32_t & matchingScore,
                const int32_t & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                const int32_t & w, const int32_t & xDrop){

    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed(seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2, windowsSize,
                     mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                     matchingScore, mismatchingPenalty, m, step_size,
                     seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, w,
                     xDrop, x_extend_seeds );
    class Matrix T(length1+1, length2 + 1);
    return findSimilarFragmentsForPairedSequence ( seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2, windowsSize,
            _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2, _extend_gap_penalty2, matchingScore,
            mismatchingPenalty, m, seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, bandwidth, x_extend_seeds, T);

}





// the third algorithm
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a weighted sequence alignment approach

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted_1gap ( int8_t * seq1, int8_t * seq1_rev_com,
                       int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                       const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                       const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,  const int32_t & matchingScore, const int32_t & mismatchingPenalty,
                       const Scorei & m, const int32_t & step_size, std::string & seq1_string, std::string & seq2_string,
                       const double & pvalues, const double & lambda, const double & kValue, const int32_t & zDrop,
                       const int32_t & bandwidth, int32_t & w, const int32_t & xDrop, Score & score,
                       int16_t * weight, int16_t * weight_rev){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t start2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;
    Matrix T(length1+1, length2 + 1);

    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed( seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2, windowsSize,
                      mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                      matchingScore, mismatchingPenalty, m, step_size,
                      seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, w,
                      xDrop, x_extend_seeds );

//    std::cout << "line 970:" << x_extend_seeds.size() << std::endl;
    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();
//        std::cout << "line 975: endPosition1:" << endPosition1 << " endPosition2:" << endPosition2  << std::endl;
        SemiGlobal_single_gap_penalty(seq1_rev_com + (length1 - 1 - endPosition1),
                   seq2_rev_com + (length2 - endPosition2 - 1),
                   1 + endPosition1, endPosition2 + 1,
                   weight_rev + (length1 - 1 - endPosition1), score, maxScore,
                   endPosition1_rc, endPosition2_rc, false, zDrop, bandwidth, T);
//        std::cout << "line 981: maxScore:" << maxScore << " endPosition1_rc:" << endPosition1_rc << " endPosition2_rc:" << endPosition2_rc << std::endl;

        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based
//        std::cout << "line 985 start1:" << start1 << " start2:" << start2 << std::endl;
        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {

            startAlignments.insert(startAlignment);

            std::vector<uint32_t> cigar = SemiGlobal_single_gap_penalty(seq1 + start1, seq2 + start2, length1 - start1, length2 - start2,
                                                     weight + start1, score, maxScore, endPosition1, endPosition2, true, zDrop, bandwidth, T);
//            std::cout << "line 993 maxScore:" << maxScore << " endPosition1:" << endPosition1 << " endPosition2:" << endPosition2 << std::endl;

            endPosition1 += start1; // 0based
            endPosition2 += start2; // 0based

            int32_t toRemove = -1;
            bool found = false;
            for ( int32_t i=0; i< pairedSimilarFragments.size(); ++i) {
                if ( pairedSimilarFragments[i].getEnd1() == endPosition1+1 && pairedSimilarFragments[i].getEnd2() == endPosition2+1 ){
                    found = true;
                    if( start1+1 < pairedSimilarFragments[i].getStart1() && start2+1 < pairedSimilarFragments[i].getStart2() ){
                        toRemove = i;
                    }
                }
            }
//            std::cout << "line 808" << std::endl;
            if( !found or toRemove>=0 ){
                eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                pvalue = 1.0 - exp(-1.0*eValue);

                PairedSimilarFragment pairedSimilarFragment(start1+1, endPosition1 + 1, start2+1, endPosition2 + 1, maxScore, cigar,
                                                            pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                // here the ordinate put into pairedSimilarFragment is 1 based
                if( toRemove>=0  ){
                    pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                }else{
                    pairedSimilarFragments.push_back(pairedSimilarFragment);
                }
            }
//            std::cout << "line 822" << std::endl;
        }
    }

    return pairedSimilarFragments;
}



// the forth algorithm
// using smith-waterman approach to find seeds and using x-drop to extend the seeds
// extend the x-drop result using a weighted 2-piece gap cost sequence alignment approach

std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence_wighted ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                                   int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                                   const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                                                                                   const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1, const int32_t & matchingScore, const int32_t & mismatchingPenalty,
                                                                                   const Scorei & m, const int32_t & step_size, std::string & seq1_string, std::string & seq2_string,
                                                                                   const double & pvalues, const double & lambda, const double & kValue, const int32_t & zDrop,
                                                                                   const int32_t & bandwidth, int32_t & w,  const int32_t & xDrop, Score & score,
                                                                                   int16_t * weight, int16_t * weight_rev){

    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t start2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;
    Matrix T(length1+1, length2 + 1);


    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed(seq1, seq1_rev_com, seq2, seq2_rev_com, reinterpret_cast<int32_t &>(length1), length2, windowsSize,
                     mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                     matchingScore, mismatchingPenalty, m, step_size,
                     seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, w,
                     xDrop, x_extend_seeds );
    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();
        SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1),
                   seq2_rev_com + (length2 - endPosition2 - 1),
                   1 + endPosition1, endPosition2 + 1,
                   weight_rev + (length1 - 1 - endPosition1), score, maxScore,
                   endPosition1_rc, endPosition2_rc, false, zDrop, bandwidth, T);
        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based

        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {

            startAlignments.insert(startAlignment);

            std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1, seq2 + start2, length1 - start1, length2 - start2,
                                                     weight + start1, score, maxScore, endPosition1, endPosition2, true, zDrop, bandwidth, T);
            if( cigar.size() > 0 && endPosition1>0 && endPosition2>0 ){
                endPosition1 += start1; // 0based
                endPosition2 += start2; // 0based
                int32_t toRemove = -1;
                bool found = false;
                for ( int32_t i=0; i< pairedSimilarFragments.size(); ++i) {
                    if ( pairedSimilarFragments[i].getEnd1() == endPosition1+1 && pairedSimilarFragments[i].getEnd2() == endPosition2+1 ){
                        found = true;
                        if( start1+1 < pairedSimilarFragments[i].getStart1() && start2+1 < pairedSimilarFragments[i].getStart2() ){
                            toRemove = i;
                        }
                    }
                }

                if( !found or toRemove>=0 ){
                    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                    pvalue = 1.0 - exp(-1.0*eValue);

                    PairedSimilarFragment pairedSimilarFragment(start1+1, endPosition1 + 1, start2+1, endPosition2 + 1, maxScore, cigar,
                                                                pvalue, eValue); // start1， start2, end1 and end2 are 1 based
                    // here the ordinate put into pairedSimilarFragment is 1 based
                    if( toRemove>=0  ){
                        pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                    }else{
                        pairedSimilarFragments.push_back(pairedSimilarFragment);
                    }
                }
            }
        }
    }
    return pairedSimilarFragments;
}


/*



// TODO
std::vector<PairedSimilarFragment> mapCNSToGenome ( int8_t * seq1, int8_t * seq1_rev_com,
                                                                           int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                                                                           const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                                                                           const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                                                                           const int32_t & _open_gap_penalty2, const int32_t & _extend_gap_penalty2, const int32_t & matchingScore,
                                                                           const int32_t & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                                                                           std::string & seq1_string, std::string & seq2_string, const double & pvalues,
                                                                           const double & lambda, const double & kValue, const int32_t & zDrop, const int32_t & bandwidth,
                                                                           const int32_t & w, const int32_t & xDrop){

    std::vector<Seed> x_extend_seeds;

    getAllExtendSeed(seq1, seq1_rev_com, seq2, seq2_rev_com, length1, length2, windowsSize,
                     mini_cns_score, matrix_boundary_distance, _open_gap_penalty1, _extend_gap_penalty1,
                     matchingScore, mismatchingPenalty, m, step_size,
                     seq1_string, seq2_string, pvalues, lambda, kValue, zDrop, bandwidth, w,
                     xDrop, x_extend_seeds );





    std::vector<PairedSimilarFragment> pairedSimilarFragments;
    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;

    int32_t maxScore;

    int32_t start1;
    int32_t end1;
    int32_t start2;
    int32_t end2;
    std::set<StartAlignment> startAlignments;
    double length1d = double(length1);
    double length2d = double(length2);
    double eValue;
    double pvalue;

    for( Seed x_seed :  x_extend_seeds){

        endPosition1 = x_seed.getEnd1();
        endPosition2 = x_seed.getEnd2();

        SemiGlobal(seq1_rev_com + (length1 - 1 - endPosition1),
                   seq2_rev_com + (length2 - endPosition2 - 1),
                   1 + endPosition1, endPosition2 + 1, _open_gap_penalty1, _extend_gap_penalty1, _open_gap_penalty2,
                   _extend_gap_penalty2, maxScore, endPosition1_rc, endPosition2_rc, m, false, zDrop, bandwidth, T);// this is the reverse alignment
//        std::cout << "line 779 " << maxScore << std::endl;
        start1 = endPosition1 - endPosition1_rc; //0 based
        start2 = endPosition2 - endPosition2_rc ; // 0 based

        StartAlignment startAlignment(start1, start2);
        if (startAlignments.find(startAlignment) == startAlignments.end()) {
            startAlignments.insert(startAlignment);
//            std::cout <<"line 786" << std::endl;
            std::vector<uint32_t> cigar = SemiGlobal(seq1 + start1, seq2 + start2, length1 - start1,
                                                     length2 - start2, _open_gap_penalty1, _extend_gap_penalty1,
                                                     _open_gap_penalty2, _extend_gap_penalty2,
                                                     maxScore, end1, end2, m, true, zDrop, bandwidth, T);
            if( cigar.size() > 0  && end1>0 && end2>0) {
                end1 += start1; // 0based
                end2 += start2; // 0based

                int32_t toRemove = -1;
                bool found = false;
                for (int32_t i = 0; i < pairedSimilarFragments.size(); ++i) {
                    if (pairedSimilarFragments[i].getEnd1() == end1 + 1 &&
                        pairedSimilarFragments[i].getEnd2() == end2 + 1) {
                        found = true;
                        if (start1 + 1 < pairedSimilarFragments[i].getStart1() &&
                            start2 + 1 < pairedSimilarFragments[i].getStart2()) {
                            toRemove = i;
                        }
                    }
                }
                if (!found or toRemove >= 0) {

                    eValue = kValue * length1d * length2d * exp(-1.0 * lambda * double(maxScore));
                    pvalue = 1.0 - exp(-1.0 * eValue);

                    PairedSimilarFragment pairedSimilarFragment(start1 + 1, end1 + 1, start2 + 1, end2 + 1, maxScore,
                                                                cigar,
                                                                pvalue,
                                                                eValue); // start1， start2, end1 and end2 are 1 based
                    // here the ordinate put into pairedSimilarFragment is 1 based
                    if (toRemove >= 0) {
                        pairedSimilarFragments[toRemove] = pairedSimilarFragment;
                    } else {
                        pairedSimilarFragments.push_back(pairedSimilarFragment);
                    }
                }
            }
        }
    }





    int32_t maxScore;
    int32_t endPosition1;
    std::cout << "mapCNSToGenome line 19, length1:" << length1 << " length2:"<< length2 << std::endl;
    mapCnsToGenome(seq1, seq2, length1, length2, _open_gap_penalty1, _extend_gap_penalty1,
                   _open_gap_penalty2, _extend_gap_penalty2,
                   maxScore, endPosition1, m);

    int32_t length11 = 1 + endPosition1;
    int32_t  endPosition1_rc;
    std::cout << "mapCNSToGenome line 23, length11:" << length11 << " length2:"<< length2 << std::endl;
    mapCnsToGenome(seq1_rev_com+(length1 - 1 - endPosition1), seq2_rev_com, length11, length2, _open_gap_penalty1, _extend_gap_penalty1,
                   _open_gap_penalty2, _extend_gap_penalty2,
                   maxScore, endPosition1_rc, m);
    std::cout << "mapCNSToGenome line 30" << std::endl;
    int32_t start1 = endPosition1 - endPosition1_rc; //0 based
    int32_t length21 = endPosition1_rc + 1;

    std::vector<uint32_t> cigar;

    std::cout << "length21:" << length21 << " length2:" << length2 << std::endl;
//    Matrix T (length21 + 1, length2 + 1);
//    std::vector<uint32_t> cigar = mapCnsToGenome(seq1+start1, seq2, length21, length2, _open_gap_penalty1, _extend_gap_penalty1,
//                                                 _open_gap_penalty2, _extend_gap_penalty2,
//                                                 maxScore, endPosition1, m, T);


    double p=0;
    PairedSimilarFragment pairedSimilarFragment(start1, endPosition1 + 1, 1, length2, maxScore, cigar, p, p);
    return pairedSimilarFragment;
}
*/