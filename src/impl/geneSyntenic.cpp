//
// Created by Baoxing Song on 2019-03-13.
//

#include "geneSyntenic.h"

/**
 * this function try to keep those genes in the syntenic region using a longest path algorithm
 * which is a kind of global alignment method
 *
 * longest path algorithm from here:
 * https://www.geeksforgeeks.org/find-longest-path-directed-acyclic-graph/
 *
 */


/**
 * my version of Quota-alignment using the longest alignment method
 * This longest alignment is designed as a local alignment model
 * And a lot of ideas were borrowed from the DAGchainder method
 * */
void myOrthologPairsSortQuota(std::vector<AlignmentMatch> &pairedSimilarFragments) {
    std::sort(pairedSimilarFragments.begin(), pairedSimilarFragments.end(), [](AlignmentMatch a, AlignmentMatch b) {
        return a < b;
    });
}

void myOrthologPairsSortQueryQuota(std::vector<AlignmentMatch> &pairedSimilarFragments) {
    std::sort(pairedSimilarFragments.begin(), pairedSimilarFragments.end(), [](AlignmentMatch a, AlignmentMatch b) {
        return a.getQueryStartPos() < b.getQueryStartPos();
    });
}

struct Path {
    double score;
    int index;
};
// since we will change pairedSimilarFragments, so do not use reference C++ data type here

/**
 * this function try to keep those genes in the syntenic region using a longest path algorithm
 * which is a kind of global alignment method
 *
 * longest path algorithm from here:
 * https://www.geeksforgeeks.org/find-longest-path-directed-acyclic-graph/
 *
 */
// this function put the forward entries in the increasing order and put the reversion entries in a decrease order
void myAlignmentMatchSort(std::vector<AlignmentMatch> &pairedSimilarFragments, const double &penalty, const double &scoreThreshold, const bool &keepTandemDuplication, const bool &considerInversion) {
    std::sort(pairedSimilarFragments.begin(), pairedSimilarFragments.end(), [](AlignmentMatch a, AlignmentMatch b) {
        return a.getRefStartPos() < b.getRefStartPos();
    });
    if (!considerInversion) {
        return;
    }

    // the following part is for reversion
    int startIndex = 0;
    int endIndex = 0;
    double maxScore = 0;
    double currentScore = 0;
    for (int idx = 0; idx < pairedSimilarFragments.size(); ++idx) {
        if (NEGATIVE == pairedSimilarFragments[idx].getStrand()) {  // reverse strand  we need to do something special for inversion
            // look for all following pairs that are reverse strand,
            // ie get all reverse strand entries in this group of reverse strands
            for (int jdx = idx; jdx < pairedSimilarFragments.size(); ++jdx) {
                if (NEGATIVE == pairedSimilarFragments[jdx].getStrand()) {
                    if (idx == jdx) { // the first one.
                        currentScore += pairedSimilarFragments[jdx].getScore();
                    } else { // for the reverse alignments, check the strand of the previous alignment
                        // If both current and previous are reverse strand, and
                        // if previous assembly start is greater than current assembly start, increase score.
                        //  else, apply penalty (they are out of order)
                        if (pairedSimilarFragments[jdx - 1].getQueryStartPos() > pairedSimilarFragments[jdx].getQueryStartPos() && pairedSimilarFragments[idx].getQueryStartPos() > pairedSimilarFragments[jdx].getQueryStartPos()) { // GOOD INVERSION
                            currentScore += pairedSimilarFragments[jdx].getScore();
                        } else if (keepTandemDuplication && pairedSimilarFragments[jdx - 1].getReferenceGeneName() == pairedSimilarFragments[jdx].getReferenceGeneName()) { // tandem duplication
                            currentScore += pairedSimilarFragments[jdx].getScore();
                        } else {
                            currentScore += penalty * pairedSimilarFragments[jdx].getScore(); // GIVE PENALTY
                        }
                    }
                } else {
                    currentScore += penalty * pairedSimilarFragments[jdx].getScore(); // penalty because are now forward strand
                }
                if (maxScore < currentScore) {
                    maxScore = currentScore;
                    endIndex = jdx; // keeps track of where to stop the reverse strand grouping
                }
                // If score is negative, stop the loop.  This will happen as we find more
                // forward vs reverse strands.  If there was just 1 reverse strand followed
                // by a forward strand, currentScore goes from 3 to -1 (3 plus -4)
                if (currentScore < 0) {
                    break;
                }
            }
            // if maxScore is larger than scoreThreshold, it means there were multiple reverse
            // strand entries in this group, or a single reverse alignment of significant length.
            // If we have several reverse alignments, we think it is real.  If just 1, may be false
            // alignment.  If the maxScore is greater than the scoreThreshold, we treat it as real.
            // Flip all elements in the range so assembly coordinates are in increasing order.
            if (maxScore > scoreThreshold) {
                maxScore = 0.0;
                currentScore = 0.0;
                // loop to find start index of elements we want to flip
                for (int jdx = endIndex; jdx >= idx; --jdx) {
                    if (NEGATIVE == pairedSimilarFragments[jdx].getStrand()) {
                        if (jdx > idx) {
                            // Verify the reverse alignments are in order to each other.  If not,
                            // apply penalty.  This doesn't prevent overlaps, which will be dealt with later.
                            if (pairedSimilarFragments[jdx - 1].getQueryStartPos() > pairedSimilarFragments[jdx].getQueryStartPos()) {
                                currentScore += pairedSimilarFragments[jdx].getScore();
                            } else if (keepTandemDuplication && pairedSimilarFragments[jdx - 1].getReferenceGeneName() == pairedSimilarFragments[jdx].getReferenceGeneName()) { // tandem duplication
                                currentScore += pairedSimilarFragments[jdx].getScore();
                            } else {
                                currentScore += penalty * pairedSimilarFragments[jdx].getScore(); // GIVE PENALTY
                            }
                        } else {
                            currentScore += pairedSimilarFragments[jdx].getScore();
                        }
                    } else {
                        currentScore += penalty * pairedSimilarFragments[jdx].getScore(); // GIVE PENALTY
                    }
                    if (maxScore < currentScore) {
                        maxScore = currentScore;
                        startIndex = jdx;
                    }
                    if (currentScore < 0) {
                        break;
                    }
                }

                int length = (endIndex - startIndex + 1) / 2;
                if (length > 0) {
                    for (int j = 0; j < length; ++j) {
                        AlignmentMatch temp = pairedSimilarFragments[startIndex + j];
                        pairedSimilarFragments[startIndex + j] = pairedSimilarFragments[endIndex - j];
                        pairedSimilarFragments[endIndex - j] = temp;
                    }
                    // Flip the elements in the list so the asm reverse alignments
                    // all have increasing asm values. (ie the last element of the reverse
                    // grouping is now the first, and the first is the last.
                    bool thereAreReverseAlignments = true;
                    while (thereAreReverseAlignments) {
                        thereAreReverseAlignments = false;
                        for (int j = 1; j < length; ++j) {
                            // If a single reference position has multiple assembly alignments mapping to it,
                            // swap the order until the assembly positions are all increasing.
                            if (pairedSimilarFragments[startIndex + j - 1].getReferenceGeneName() == pairedSimilarFragments[startIndex + j].getReferenceGeneName() &&
                                pairedSimilarFragments[startIndex + j - 1].getQueryStartPos() > pairedSimilarFragments[startIndex + j].getQueryStartPos()) {
                                thereAreReverseAlignments = true;
                                AlignmentMatch temp = pairedSimilarFragments[startIndex + j];
                                pairedSimilarFragments[startIndex + j] = pairedSimilarFragments[startIndex + j - 1];
                                pairedSimilarFragments[startIndex + j - 1] = temp;
                            }
                        }
                    }
                }
                // Flip the elements in the list so the asm reverse ailgnments
                // all have increasing asm values. (ie the last element of the reverse
                // grouping is now the first, and the first is the last.

                idx = endIndex;
            }
            maxScore = 0.0;
            currentScore = 0.0;
        }
    }
}

void longestPath(std::vector<AlignmentMatch> &pairedSimilarFragments, std::vector<AlignmentMatch> &sortedOrthologPairs, const bool &keepTandemDuplication, double &scoreThreshold) { //double checked used for self to self alignment
    double maxSore = 0;
    int bestEnd = 0;
    double scoreArray[pairedSimilarFragments.size()]; // arrays of scores
    int prev[pairedSimilarFragments.size()];  // index of previous node in longest path
    scoreArray[0] = pairedSimilarFragments[0].getScore();
    prev[0] = -1;

    if (scoreArray[0] > maxSore) {
        bestEnd = 0;
        maxSore = scoreArray[0];
    }

    for (int idx = 1; idx < pairedSimilarFragments.size(); ++idx) {
        scoreArray[idx] = pairedSimilarFragments[idx].getScore();
        prev[idx] = -1;
        for (int jdx = idx - 1; jdx >= 0; --jdx) {// checking all previous nodes
            // Because we swapped asm/query start position so that inversions were all increasing,
            // we should always be on the diagonal.  If not, then we filter it.
            // This gets rid of the noise, while preserving the inversions on
            // the diagonal
            // Are only looking at positions previous to our current "idx" position
            if ((scoreArray[jdx] + pairedSimilarFragments[idx].getScore()) > scoreArray[idx] &&
                pairedSimilarFragments[jdx].getQueryEndPos() < pairedSimilarFragments[idx].getQueryStartPos()) {
                if (pairedSimilarFragments[jdx].getStrand() == POSITIVE && pairedSimilarFragments[idx].getStrand() == POSITIVE
                    && pairedSimilarFragments[jdx].getRefEndPos() < pairedSimilarFragments[idx].getRefStartPos()
                    && pairedSimilarFragments[jdx].getQueryEndPos() < pairedSimilarFragments[idx].getQueryStartPos()) {
                    scoreArray[idx] = scoreArray[jdx] + pairedSimilarFragments[idx].getScore();
                    prev[idx] = jdx;
                } else if (pairedSimilarFragments[jdx].getStrand() == NEGATIVE && pairedSimilarFragments[idx].getStrand() == NEGATIVE
                           && pairedSimilarFragments[jdx].getRefStartPos() > pairedSimilarFragments[idx].getRefEndPos()
                           && pairedSimilarFragments[jdx].getQueryEndPos() < pairedSimilarFragments[idx].getQueryStartPos()) {
                    scoreArray[idx] = scoreArray[jdx] + pairedSimilarFragments[idx].getScore();
                    prev[idx] = jdx;
                } else if (pairedSimilarFragments[jdx].getStrand() == NEGATIVE && pairedSimilarFragments[idx].getStrand() == NEGATIVE
                           && pairedSimilarFragments[jdx].getRefEndPos() < pairedSimilarFragments[idx].getRefStartPos()
                           && pairedSimilarFragments[jdx].getQueryEndPos() < pairedSimilarFragments[idx].getQueryStartPos() && scoreArray[jdx] > scoreThreshold) { // this block is under test, todo remind PHG guies
                    scoreArray[idx] = scoreArray[jdx] + pairedSimilarFragments[idx].getScore();
                    prev[idx] = jdx;
                } else if (pairedSimilarFragments[jdx].getStrand() != pairedSimilarFragments[idx].getStrand()
                           && pairedSimilarFragments[jdx].getRefEndPos() < pairedSimilarFragments[idx].getRefStartPos()
                           && pairedSimilarFragments[jdx].getQueryEndPos() < pairedSimilarFragments[idx].getQueryStartPos()) {
                    scoreArray[idx] = scoreArray[jdx] + pairedSimilarFragments[idx].getScore()/*-0.505*/ ;
                    prev[idx] = jdx;
                }
            }
        }
        if (scoreArray[idx] > maxSore) {
            bestEnd = idx;
            maxSore = scoreArray[idx];
        }
    }
//    std::cout << "maxSore: " << maxSore << std::endl;
    int idx = bestEnd; // bestEnd is where to stop the longest path
    sortedOrthologPairs.push_back(pairedSimilarFragments[idx]);
    int jdx = prev[idx]; // prev[] is index on the longest path
    while (jdx >= 0) {
        sortedOrthologPairs.push_back(pairedSimilarFragments[jdx]);
        jdx = prev[jdx];
    }
    // Reversing the order
    std::reverse(std::begin(sortedOrthologPairs), std::end(sortedOrthologPairs));
    std::vector<int> toRemoveIndex;
    std::vector<int> thisRoundOfInversions;
    toRemoveIndex.clear();
    thisRoundOfInversions.clear();
    double score = 0;
//    std::cout << "line 238 sortedOrthologPairs size:" << sortedOrthologPairs.size() << std::endl;
    for (int32_t i = 0; i < sortedOrthologPairs.size(); ++i) {

        if (sortedOrthologPairs[i].getStrand() == NEGATIVE &&
            (i == 0 || (sortedOrthologPairs[i - 1].getStrand() == POSITIVE) ||
             (sortedOrthologPairs[i - 1].getStrand() == NEGATIVE && sortedOrthologPairs[i - 1].getRefStartPos() > sortedOrthologPairs[i].getRefEndPos()
              && sortedOrthologPairs[i - 1].getQueryEndPos() < sortedOrthologPairs[i].getQueryStartPos())
            )
                ) {
//            std::cout << "line 260" << std::endl;
            score += sortedOrthologPairs[i].getScore();
            thisRoundOfInversions.push_back(i);
        } else {
            if (score <= scoreThreshold) {
                for (int j: thisRoundOfInversions) {
                    toRemoveIndex.push_back(j);
//                    std::cout << "line 268, remove:" << j << std::endl;
                }
            }
            thisRoundOfInversions.clear();
            score = 0;
        }
    }

    // in case the last one is an inversion begin
    if (score > 0 && score <= scoreThreshold) {
        for (int j: thisRoundOfInversions) {
            toRemoveIndex.push_back(j);
        }
    }
    thisRoundOfInversions.clear();
    score = 0;
    // in case the last one is an inversion end

    for (int j = toRemoveIndex.size() - 1; j >= 0; j--) {
        sortedOrthologPairs.erase(sortedOrthologPairs.begin() + toRemoveIndex[j]);
    }
    //  std::cout << "line 259 sortedOrthologPairs size:" << sortedOrthologPairs.size() << std::endl;
//    for( int i=0; i<sortedOrthologPairs.size(); ++i ){
//
//        std::cout << "line:287\t" << sortedOrthologPairs[i].getRefChr() << "\t"
//                  << sortedOrthologPairs[i].getRefStartPos() << "\t"
//                  << sortedOrthologPairs[i].getRefEndPos() << "\t"
//                  << sortedOrthologPairs[i].getQueryChr() << "\t"
//                  << sortedOrthologPairs[i].getQueryStartPos() << "\t"
//                  << sortedOrthologPairs[i].getQueryEndPos() << "\t"
//                  << sortedOrthologPairs[i].getStrand() << "\t"
//                  << sortedOrthologPairs[i].getReferenceGeneName() << "\t" <<
//                  i <<"\t" << sortedOrthologPairs[i].getScore() << std::endl;
//    }
}


// since we will change pairedSimilarFragments, so do not use reference C++ data type here
void longestPathQuotav2(std::vector<AlignmentMatch> pairedSimilarFragments, std::vector<std::vector<AlignmentMatch>> &sortedOrthologPairChains,
                        std::map<std::string, std::map<int64_t, AlignmentMatch>> &refIndexMap /*chr, index, AlignmentMatch*/, std::map<std::string, std::map<int64_t, AlignmentMatch>> &queryIndexMap,
                        double &INDEL_SCORE, double &GAP_OPEN_PENALTY,
                        double &MIN_ALIGNMENT_SCORE, const int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                        double &calculateIndelDistance, bool withNovelAnchors) {
    sortedOrthologPairChains.clear();
    std::map<std::string, int64_t> refTimes;
    std::map<std::string, std::map<int64_t, int64_t>> queryTimes;  // key is the id set above

    std::vector<Path> high;
    std::vector<int64_t> ans;
    bool done;
    int n, i, j;
    do {
        done = true;
        n = pairedSimilarFragments.size();
//        std::cout << "longest path n:" << n << std::endl;
        double scoreArray[n]; // arrays of scores
        int prev[n];  // index of previous node in longest path

        for (int idx = 0; idx < n; ++idx) {
            scoreArray[idx] = pairedSimilarFragments[idx].getScore();
            prev[idx] = -1;
        }

        for (int idx = 1; idx < n; ++idx) {
            double thisIndexScore = scoreArray[idx];
            for (int jdx = idx - 1; jdx >= 0; --jdx) {// checking all previous nodes
                // Because we swapped asm/query start position so that inversions were all increasing,
                // we should always be on the diagonal.  If not, then we filter it.
                // This gets rid of the noise, while preserving the inversions on
                // the diagonal
                // Are only looking at positions previous to our current "idx" position
                if (pairedSimilarFragments[idx].getQueryChr() == pairedSimilarFragments[jdx].getQueryChr()
                    && pairedSimilarFragments[idx].getRefChr() == pairedSimilarFragments[jdx].getRefChr()) {
                    if (pairedSimilarFragments[idx].getStrand() == pairedSimilarFragments[jdx].getStrand()) {
                        // the node one the chain should be in the same STRAND, if not this is an INDEL
                        if (pairedSimilarFragments[idx].getStrand() == POSITIVE && pairedSimilarFragments[idx].getQueryId() > pairedSimilarFragments[jdx].getQueryId() &&
                            pairedSimilarFragments[jdx].getQueryEndPos() < pairedSimilarFragments[idx].getQueryStartPos() &&
                            pairedSimilarFragments[jdx].getRefEndPos() < pairedSimilarFragments[idx].getRefStartPos()) { //same strand
//                            int ref_del = pairedSimilarFragments[idx].getRefId() - pairedSimilarFragments[jdx].getRefId()-1;
//                            int query_del = pairedSimilarFragments[idx].getQueryId() - pairedSimilarFragments[jdx].getQueryId()-1;

//                            int ref_del = 0;
//                            for( int idi = pairedSimilarFragments[jdx].getRefId()+1; idi<pairedSimilarFragments[idx].getRefId(); ++idi){
//                                ref_del += (refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getRefEndPos() - refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getRefStartPos() + 1 );
//                            }
//
//                            int query_del =  0;
//                            for( int idi = pairedSimilarFragments[jdx].getQueryId()+1; idi<pairedSimilarFragments[idx].getQueryId(); ++idi){
//                                query_del += (queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryEndPos() - queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryStartPos() + 1 );
//                            }

                            double ref_del = 0;
//                            int smallerIdi = pairedSimilarFragments[jdx].getRefId();
//                            int largerIdi = pairedSimilarFragments[idx].getRefId();
//                            if ( smallerIdi > largerIdi ){
//                                int temp = smallerIdi;
//                                smallerIdi = largerIdi;
//                                largerIdi = temp;
//                            }
//                            for( int idi = smallerIdi+1; idi<largerIdi; ++idi){
                            for (int idi = pairedSimilarFragments[jdx].getRefId() + 1; idi < pairedSimilarFragments[idx].getRefId(); ++idi) {
                                if (refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName().find("localAlignment") == std::string::npos) {
                                    ++ref_del;
                                    if (refTimes.find(refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName()) != refTimes.end() && refTimes[refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName()] >= refMaximumTimes) {
                                        ref_del += MAX_DIST_BETWEEN_MATCHES + MAX_DIST_BETWEEN_MATCHES;
                                    }
                                }
                            }

//                            smallerIdi = pairedSimilarFragments[jdx].getQueryId();
//                            largerIdi = pairedSimilarFragments[idx].getQueryId();
//                            if ( smallerIdi > largerIdi ){
//                                int temp = smallerIdi;
//                                smallerIdi = largerIdi;
//                                largerIdi = temp;
//                            }
                            double query_del = 0;
//                            for( int idi = smallerIdi+1; idi<largerIdi; ++idi){
                            for (int idi = pairedSimilarFragments[jdx].getQueryId() + 1; idi < pairedSimilarFragments[idx].getQueryId(); ++idi) {
                                if (queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryGeneName().find("localAlignment") == std::string::npos) {
                                    ++query_del;
//                                    std::cout << "line 410: " << pairedSimilarFragments[idx].getQueryChr() << "\t" << queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getReferenceGeneName() << " pairedSimilarFragments[idx].getQueryId():" << pairedSimilarFragments[idx].getQueryId() << " pairedSimilarFragments[jdx].getQueryId():" << pairedSimilarFragments[jdx].getQueryId() <<" idi:" << idi << " queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryId():" << queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryId() << std::endl;
                                    assert(queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryId() == idi);
                                    if (queryTimes.find(pairedSimilarFragments[idx].getQueryChr()) != queryTimes.end()
                                        && queryTimes[pairedSimilarFragments[idx].getQueryChr()].find(idi) != queryTimes[pairedSimilarFragments[idx].getQueryChr()].end()
                                        && queryTimes[pairedSimilarFragments[idx].getQueryChr()][idi] >= queryMaximumTimes) {
                                        query_del += MAX_DIST_BETWEEN_MATCHES + MAX_DIST_BETWEEN_MATCHES;
                                    }
                                }
                            }

                            assert(ref_del >= 0);
                            assert(query_del >= 0);
                            double distance = (((ref_del + query_del) + std::abs(ref_del - query_del)) / (calculateIndelDistance));
//                            if( withNovelAnchros ){
//                                distance = (ref_del+query_del)/2;
//                            }

                            if (std::abs(ref_del) > MAX_DIST_BETWEEN_MATCHES || std::abs(query_del) > MAX_DIST_BETWEEN_MATCHES) {
                                // if this position is too large, then the last node of j could not be i and the chain restart from j
                                break;
                            }
                            if (std::abs(std::abs(ref_del) - std::abs(query_del)) > MAX_DIST_BETWEEN_MATCHES) {
                                break;
                            }

                            double thisScore = thisIndexScore;
                            if (pairedSimilarFragments[idx].getRefId() != pairedSimilarFragments[jdx].getRefId()) {
                                thisScore += scoreArray[jdx];
                            }
                            if (distance > 0) {
                                if (pairedSimilarFragments[idx].getQueryGeneName().find("localAlignment") != std::string::npos || pairedSimilarFragments[jdx].getQueryGeneName().find("localAlignment") != std::string::npos) {
                                    thisScore = thisScore + INDEL_SCORE * distance;
                                } else {
                                    thisScore = thisScore + GAP_OPEN_PENALTY + INDEL_SCORE * distance;
                                }
                            }

//                            std::cout << "line 599 thisScore:" << thisScore << " distance:" << distance << " jdx:" << jdx << " idx:" << idx << " thisIndexScore:" << thisIndexScore << " scoreArray[jdx]:" << scoreArray[jdx] << " pairedSimilarFragments[idx].getRefId()" << pairedSimilarFragments[idx].getRefId() << " pairedSimilarFragments[jdx].getRefId():" << pairedSimilarFragments[jdx].getRefId() << std::endl;
                            if (thisScore > scoreArray[idx] && pairedSimilarFragments[jdx].getQueryId() < pairedSimilarFragments[idx].getQueryId()) {
                                scoreArray[idx] = thisScore;
                                prev[idx] = jdx;
                            }
                        } else if (pairedSimilarFragments[idx].getStrand() == NEGATIVE && pairedSimilarFragments[jdx].getQueryId() > pairedSimilarFragments[idx].getQueryId() &&
                                   pairedSimilarFragments[jdx].getQueryStartPos() > pairedSimilarFragments[idx].getQueryEndPos() &&
                                   pairedSimilarFragments[jdx].getRefEndPos() < pairedSimilarFragments[idx].getRefStartPos()

                                ) { // inversion
//                           int ref_del = pairedSimilarFragments[idx].getRefId() - pairedSimilarFragments[jdx].getRefId()-1;
//                           int query_del = pairedSimilarFragments[jdx].getQueryId() - pairedSimilarFragments[idx].getQueryId()-1;

//                            int ref_del = 0;
//                            for( int idi = pairedSimilarFragments[jdx].getRefId()+1; idi<pairedSimilarFragments[idx].getRefId(); ++idi){
//                                ref_del += (refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getRefEndPos() - refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getRefStartPos() + 1 );
//                            }
//                            int query_del =  0;
//                            for( int idi = pairedSimilarFragments[idx].getQueryId()+1; idi<pairedSimilarFragments[jdx].getQueryId(); ++idi){
//                                query_del += (queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryEndPos() - queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryStartPos() + 1 );
//                            }
/*
                            double ref_del = 0;
                            for( int idi = pairedSimilarFragments[jdx].getRefId()+1; idi<pairedSimilarFragments[idx].getRefId(); ++idi){
                                if( refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName().find("localAlignment") == std::string::npos ){
                                    ++ref_del;
                                }
                            }

                            double query_del =  0;
                            for( int idi = pairedSimilarFragments[idx].getQueryId()+1; idi<pairedSimilarFragments[jdx].getQueryId(); ++idi){
                                if( queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryGeneName().find("localAlignment") == std::string::npos ){
                                    ++query_del;
                                }
                            }
*/
                            double ref_del = 0;
//                            int smallerIdi = pairedSimilarFragments[jdx].getRefId();
//                            int largerIdi = pairedSimilarFragments[idx].getRefId();
//                            if ( smallerIdi > largerIdi ){
//                                int temp = smallerIdi;
//                                smallerIdi = largerIdi;
//                                largerIdi = temp;
//                            }
//                            for( int idi = smallerIdi+1; idi<largerIdi; ++idi){
                            for (int idi = pairedSimilarFragments[jdx].getRefId() + 1; idi < pairedSimilarFragments[idx].getRefId(); ++idi) {
                                if (refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName().find("localAlignment") == std::string::npos) {
                                    ++ref_del;

                                    if (refTimes.find(refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName()) != refTimes.end() && refTimes[refIndexMap[pairedSimilarFragments[idx].getRefChr()][idi].getReferenceGeneName()] >= refMaximumTimes) {
                                        ref_del += MAX_DIST_BETWEEN_MATCHES + MAX_DIST_BETWEEN_MATCHES;
                                    }
                                }
                            }
//
//                            smallerIdi = pairedSimilarFragments[jdx].getQueryId();
//                            largerIdi = pairedSimilarFragments[idx].getQueryId();
//                            if ( smallerIdi > largerIdi ){
//                                int temp = smallerIdi;
//                                smallerIdi = largerIdi;
//                                largerIdi = temp;
//                            }
                            double query_del = 0;
                            for (int idi = pairedSimilarFragments[idx].getQueryId() + 1; idi < pairedSimilarFragments[jdx].getQueryId(); ++idi) {
//                            for( int idi = smallerIdi+1; idi<largerIdi; ++idi){
                                if (queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryGeneName().find("localAlignment") == std::string::npos) {
                                    ++query_del;
//                                    std::cout << pairedSimilarFragments[idx].getQueryChr() << "\t" << queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getReferenceGeneName() << " pairedSimilarFragments[idx].getQueryId():" << pairedSimilarFragments[idx].getQueryId() << " pairedSimilarFragments[jdx].getQueryId():" << pairedSimilarFragments[jdx].getQueryId() <<" idi:" << idi << " queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryId():" << queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryId() << std::endl;
                                    assert(queryIndexMap[pairedSimilarFragments[idx].getQueryChr()][idi].getQueryId() == idi);
                                    if (queryTimes.find(pairedSimilarFragments[idx].getQueryChr()) != queryTimes.end()
                                        && queryTimes[pairedSimilarFragments[idx].getQueryChr()].find(idi) != queryTimes[pairedSimilarFragments[idx].getQueryChr()].end()
                                        && queryTimes[pairedSimilarFragments[idx].getQueryChr()][idi] >= queryMaximumTimes) {
                                        query_del += MAX_DIST_BETWEEN_MATCHES + MAX_DIST_BETWEEN_MATCHES;
                                    }
                                }
                            }

                            assert(ref_del >= 0);
                            assert(query_del >= 0);
                            double distance = (((ref_del + query_del) + std::abs(ref_del - query_del)) / (calculateIndelDistance));
//                            if( withNovelAnchros ){
//                                distance = (ref_del+query_del)/2;
//                            }
                            if (std::abs(ref_del) > MAX_DIST_BETWEEN_MATCHES || std::abs(query_del) > MAX_DIST_BETWEEN_MATCHES) {
                                // if this position is too large then the last node of j could not be i and the chain restart from j
                                break;
                            }
                            if (std::abs(std::abs(ref_del) - std::abs(query_del)) > MAX_DIST_BETWEEN_MATCHES) {
                                break;
                            }

                            double thisScore = thisIndexScore;
                            if (pairedSimilarFragments[idx].getRefId() != pairedSimilarFragments[jdx].getRefId()) {
                                thisScore += scoreArray[jdx];
                            }
                            if (distance > 0) {
                                if (pairedSimilarFragments[idx].getQueryGeneName().find("localAlignment") != std::string::npos || pairedSimilarFragments[jdx].getQueryGeneName().find("localAlignment") != std::string::npos) {
                                    thisScore = thisScore + INDEL_SCORE * distance;
                                } else {
                                    thisScore = thisScore + GAP_OPEN_PENALTY + INDEL_SCORE * distance;
                                }
                            }
//                            std::cout << "line 544 thisScore:" << thisScore << " distance:" << distance << " jdx:" << jdx << " idx:" << idx << " thisIndexScore:" << thisIndexScore << " scoreArray[jdx]:" << scoreArray[jdx] << " pairedSimilarFragments[idx].getRefId()" << pairedSimilarFragments[idx].getRefId() << " pairedSimilarFragments[jdx].getRefId():" << pairedSimilarFragments[jdx].getRefId() << std::endl;
                            if (thisScore > scoreArray[idx] &&
                                pairedSimilarFragments[jdx].getQueryId() > pairedSimilarFragments[idx].getQueryId()) {
                                scoreArray[idx] = thisScore;
                                prev[idx] = jdx;
                            }
                        }
                    }
                }
            }
        }
//        std::cout << "line 645" << std::endl;
        high.clear();
        for (int idx = 0; idx < n; ++idx) {
//            std::cout << scoreArray[idx] << "\t" << MIN_ALIGNMENT_SCORE << std::endl;
            if (scoreArray[idx] > MIN_ALIGNMENT_SCORE /*&& prev[idx] != -2*/) {
                Path p;
                p.index = idx;
                p.score = scoreArray[idx];
                high.push_back(p);
            }
        }
//        std::cout << "line 618" << std::endl;
        std::sort(high.begin(), high.end(), [](Path a, Path b) {
            return a.score > b.score;
        });
//        std::cout << "line 376 high.size(): " << high.size() << " highest value " << high[0].score << "highest index:" << high[0].index << std::endl;
        if (high.size() > 0) {
            done = false;
            i = 0;
            if (prev[high[i].index] != -2) { // only output the chain with highest score for each loop
//                std::cout << "line 574\t" << high[i].score << "\t" << pairedSimilarFragments[high[i].index].getReferenceGeneName() << std::endl;
                // trace back begin
//                double subMaxScore = 0.0;
//                int startIndex = high[i].index;
//                double currentScore=0.0;
//                for  (j = high[i].index;  prev[j] >= 0;  j = prev[j]) {
//                    std::cout << j << "\t";
//                    double scoreDifference = scoreArray[j]-scoreArray[prev[j]];
//                    currentScore+=scoreDifference;
//                    if( subMaxScore < currentScore ){
//                        subMaxScore=currentScore;
//                        startIndex=j;
//                    }
//                }
//                double scoreDifference = scoreArray[j];  // this one should be the score of the node itself
//                currentScore+=scoreDifference;
//                if( subMaxScore < currentScore ){
//                    subMaxScore=currentScore;
//                    startIndex=j;
//                }
//
//                assert(abs(subMaxScore-high[i].score)<0.01);
                // trace back end
                ans.clear();
                j = high[i].index;
                if (prev[j] >= 0) {
                    while (j >= 0) {
                        ans.push_back(j);
                        j = prev[j];
                    }
                } else {
                    ans.push_back(high[i].index);
                }
                std::vector<AlignmentMatch> chain;
                sortedOrthologPairChains.push_back(chain);
                reverse(ans.begin(), ans.end());
                int s = ans.size();

                uint32_t chainRefStart = std::numeric_limits<uint32_t>::max();
                uint32_t chainRefEnd = 0;
                uint32_t chainQueryStart = std::numeric_limits<uint32_t>::max();
                uint32_t chainQueryEnd = 0;
//                std::cout << "line 426" << std::endl;
                for (j = 0; j < s; j++) {
                    prev[ans[j]] = -2;
//                    std::cout << "line 429" << std::endl;
                    AlignmentMatch orthologPair2 = pairedSimilarFragments[ans[j]];
                    sortedOrthologPairChains[sortedOrthologPairChains.size() - 1].push_back(orthologPair2);

                    chainRefStart = chainRefStart < orthologPair2.getRefStartPos() ? chainRefStart : orthologPair2.getRefStartPos();
                    chainQueryStart = chainQueryStart < orthologPair2.getQueryStartPos() ? chainQueryStart : orthologPair2.getQueryStartPos();
                    chainRefEnd = chainRefEnd > orthologPair2.getRefEndPos() ? chainRefEnd : orthologPair2.getRefEndPos();
                    chainQueryEnd = chainQueryEnd > orthologPair2.getQueryEndPos() ? chainQueryEnd : orthologPair2.getQueryEndPos();
                }

                std::string lastRef = "";
                std::set<std::string> countedRefs;
                for (int ii = 0; ii < n; ii++) {
                    if (pairedSimilarFragments[ii].getRefChr() == pairedSimilarFragments[ans[0]].getRefChr() && (
                            (pairedSimilarFragments[ii].getRefStartPos() <= chainRefStart && chainRefStart <= pairedSimilarFragments[ii].getRefEndPos()) ||
                            (pairedSimilarFragments[ii].getRefStartPos() <= chainRefEnd && chainRefEnd <= pairedSimilarFragments[ii].getRefEndPos()) ||
                            (chainRefStart <= pairedSimilarFragments[ii].getRefStartPos() && pairedSimilarFragments[ii].getRefStartPos() <= chainRefEnd) ||
                            (chainRefStart <= pairedSimilarFragments[ii].getRefEndPos() && pairedSimilarFragments[ii].getRefEndPos() <= chainRefEnd))
                        /*&& pairedSimilarFragments[ii].getReferenceGeneName() != lastRef*/  ) {
                        lastRef = pairedSimilarFragments[ii].getReferenceGeneName();
                        if (countedRefs.find(pairedSimilarFragments[ii].getReferenceGeneName()) == countedRefs.end()) { // do not count it repeatedly
                            if (refTimes.find(pairedSimilarFragments[ii].getReferenceGeneName()) != refTimes.end()) {
                                refTimes[pairedSimilarFragments[ii].getReferenceGeneName()] =
                                        refTimes[pairedSimilarFragments[ii].getReferenceGeneName()] + 1;
                            } else {
                                refTimes[pairedSimilarFragments[ii].getReferenceGeneName()] = 1;
                            }
                            countedRefs.insert(pairedSimilarFragments[ii].getReferenceGeneName());
                        }
                    }
                }
                for (int ii = 0; ii < n; ii++) {
                    if (pairedSimilarFragments[ii].getQueryChr() == pairedSimilarFragments[ans[0]].getQueryChr() && (
                            (pairedSimilarFragments[ii].getQueryStartPos() <= chainQueryStart && chainQueryStart <= pairedSimilarFragments[ii].getQueryEndPos()) ||
                            (pairedSimilarFragments[ii].getQueryStartPos() <= chainQueryEnd && chainQueryEnd <= pairedSimilarFragments[ii].getQueryEndPos()) ||
                            (chainQueryStart <= pairedSimilarFragments[ii].getQueryStartPos() && pairedSimilarFragments[ii].getQueryStartPos() <= chainQueryEnd) ||
                            (chainQueryStart <= pairedSimilarFragments[ii].getQueryEndPos() && pairedSimilarFragments[ii].getQueryEndPos() <= chainQueryEnd))
                            ) {
                        if (queryTimes.find(pairedSimilarFragments[ii].getQueryChr()) == queryTimes.end()) {
                            std::map<int64_t, int64_t> t;
                            queryTimes[pairedSimilarFragments[ii].getQueryChr()] = t;
                        }

                        if (queryTimes[pairedSimilarFragments[ii].getQueryChr()].find(pairedSimilarFragments[ii].getQueryId()) != queryTimes[pairedSimilarFragments[ii].getQueryChr()].end()) {
                            queryTimes[pairedSimilarFragments[ii].getQueryChr()][pairedSimilarFragments[ii].getQueryId()] = queryTimes[pairedSimilarFragments[ii].getQueryChr()][pairedSimilarFragments[ii].getQueryId()] + 1;
                        } else {
                            queryTimes[pairedSimilarFragments[ii].getQueryChr()][pairedSimilarFragments[ii].getQueryId()] = 1;
                        }
                    }
                }

                for (int ii = 0; ii < n; ii++) { // this is to avoid nested chain
                    if (pairedSimilarFragments[ii].getQueryChr() == pairedSimilarFragments[ans[0]].getQueryChr() && (
                            (pairedSimilarFragments[ii].getQueryStartPos() <= chainQueryStart && chainQueryStart <= pairedSimilarFragments[ii].getQueryEndPos()) ||
                            (pairedSimilarFragments[ii].getQueryStartPos() <= chainQueryEnd && chainQueryEnd <= pairedSimilarFragments[ii].getQueryEndPos()) ||
                            (chainQueryStart <= pairedSimilarFragments[ii].getQueryStartPos() && pairedSimilarFragments[ii].getQueryStartPos() <= chainQueryEnd) ||
                            (chainQueryStart <= pairedSimilarFragments[ii].getQueryEndPos() && pairedSimilarFragments[ii].getQueryEndPos() <= chainQueryEnd))
                        && pairedSimilarFragments[ii].getRefChr() == pairedSimilarFragments[ans[0]].getRefChr() && (
                                (pairedSimilarFragments[ii].getRefStartPos() <= chainRefStart && chainRefStart <= pairedSimilarFragments[ii].getRefEndPos()) ||
                                (pairedSimilarFragments[ii].getRefStartPos() <= chainRefEnd && chainRefEnd <= pairedSimilarFragments[ii].getRefEndPos()) ||
                                (chainRefStart <= pairedSimilarFragments[ii].getRefStartPos() && pairedSimilarFragments[ii].getRefStartPos() <= chainRefEnd) ||
                                (chainRefStart <= pairedSimilarFragments[ii].getRefEndPos() && pairedSimilarFragments[ii].getRefEndPos() <= chainRefEnd))
                            ) {
                        prev[ii] = -2;
                    }
                }
            }
        }
        if (!done) {
            for (i = j = 0; i < n; i++) {
                if (prev[i] != -2 &&
                    (refTimes.find(pairedSimilarFragments[i].getReferenceGeneName()) == refTimes.end() || refTimes[pairedSimilarFragments[i].getReferenceGeneName()] < refMaximumTimes) &&
                    (queryTimes.find(pairedSimilarFragments[i].getQueryChr()) == queryTimes.end() ||
                     (queryTimes[pairedSimilarFragments[i].getQueryChr()].find(pairedSimilarFragments[i].getQueryId()) == queryTimes[pairedSimilarFragments[i].getQueryChr()].end()) ||
                     queryTimes[pairedSimilarFragments[i].getQueryChr()][pairedSimilarFragments[i].getQueryId()] < queryMaximumTimes)) {
                    if (i != j) { // those elements should be maintained for next loop
                        pairedSimilarFragments[j] = pairedSimilarFragments[i];
                    }
                    ++j;
                } else if (prev[i] != -2) {
//                    std::cout << "line 694, times" << pairedSimilarFragments[i].getReferenceGeneName() << "\t" << refTimes[pairedSimilarFragments[i].getReferenceGeneName()] << "\t" << pairedSimilarFragments[i].getQueryChr() << "\t" << pairedSimilarFragments[i].getQueryStartPos() << "\t" << pairedSimilarFragments[i].getQueryEndPos() << "\t" << queryTimes[pairedSimilarFragments[i].getQueryChr()][pairedSimilarFragments[i].getQueryId()] << std::endl;
                }
            }
            pairedSimilarFragments.resize(j);
//            std::cout << "line 703 size: " <<  pairedSimilarFragments.size() << std::endl;
        }
    } while (!done);
}
