
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


// same with the above one, but return cigar
std::vector<uint32_t > x_extend_seed(int32_t & start1, int32_t & end1, int32_t & start2, int32_t & end2,
                   int32_t & maxScore, int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com,
                   const double & length1d, const double & length2d, const int32_t & length1,
                   const int32_t & length2,const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,
                   const int32_t & xdrop, const int32_t & w, const Scorei & m, class Matrix & T){

    int32_t endPosition1;
    int32_t endPosition2;
    int32_t endPosition1_rc;
    int32_t endPosition2_rc;
    int32_t s1 = start1;
    int32_t s2 = start2;
    SemiGlobal_xextend(seq1+start1, seq2+start2, length1-start1, length2-start2, _open_gap_penalty1, _extend_gap_penalty1,
                       maxScore, endPosition1, endPosition2, m, xdrop, w);

    if( endPosition1<0 || endPosition2 <0 ){
        maxScore = 0;
        std::vector<uint32_t > cigar;
        return cigar;
    }
    int32_t iii, jjj;
    std::vector<uint32_t > cigar = SemiGlobal_xextend(seq1_rev_com + (length1 - 1 - (start1+endPosition1)),
                       seq2_rev_com + (length2 - 1 - (start2 + endPosition2)),
                       start1+1 + endPosition1, start2+1+endPosition2, _open_gap_penalty1, _extend_gap_penalty1, maxScore,
                       endPosition1_rc, endPosition2_rc, m,  xdrop, w, T, iii, jjj);// this is the reverse alignment
    if( endPosition1_rc<0 || endPosition2_rc <0 ){
        maxScore = 0;
        return cigar;
    }
    end1 = s1 + endPosition1; // 0based
    end2 = s2 + endPosition2; //0-based
    start1 = s1+endPosition1-endPosition1_rc; //0-based coordinate
    start2 = s2+endPosition2-endPosition2_rc; //0-based coordinate
    end1 -= iii;
    end2 -= jjj;
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


// the first alignment algorithm
// using smith-waterman approach to find seeds and using x-drop to extend the seeds, return result for output
std::vector<PairedSimilarFragment> findSimilarFragmentsForPairedSequence ( int8_t * seq1, int8_t * seq1_rev_com,
                   int8_t * seq2, int8_t * seq2_rev_com, int32_t & length1, int32_t & length2, int32_t & windowsSize,
                   const int32_t & mini_cns_score, const int32_t & matrix_boundary_distance,
                   const int32_t & _open_gap_penalty1, const int32_t & _extend_gap_penalty1,const int32_t & matchingScore,
                   const int32_t & mismatchingPenalty, const Scorei & m, const int32_t & step_size,
                   std::string & seq1_string, std::string & seq2_string, const int32_t & scoreThreshold,
                    const int32_t & w, const int32_t & xDrop){
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
//    double eValue;
//    double pvalue;

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
//    std::cout << "seeds size" << seeds.size() << std::endl;
//    std::cout << "line 461 begin to link seeds" << std::endl;
    link_seeds( seeds );  // this function is very fast, by using less seeds, the program is significantly faster


    class Matrix T (length1 + 1, length2 + 1);
    for (Seed seed : seeds) {
        start1 = seed.getStart1()+1;
        end1 = seed.getEnd1()+1;
        start2 = seed.getStart2()+1;
        end2 = seed.getEnd2()+1;
        std::vector<uint32_t > cigar = x_extend_seed(start1, end1, start2, end2, maxScore, seq1, seq1_rev_com, seq2, seq2_rev_com,
                                        length1d, length2d, length1, length2, _open_gap_penalty1, _extend_gap_penalty1,
                                                xDrop, w, m,  T);

        if (maxScore >= scoreThreshold ) {
            bool never_used = true;
            PairedSimilarFragment pairedSimilarFragment(start1+1, end1+1, start2+1, end2+1, maxScore, cigar);
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



