//
// Created by Baoxing song on 20.10.18.
//

#include "deNovoGenomeVariantCalling.h"



void outputLocalAlignment(const int & chrWidth, const std::string & refFileName, const std::string & queryFileName,
                          std::vector<PairedSimilarFragment> pairedSimilarFragments0 , STRAND & strand,
                          const int & startRef, const int & startQuery, const std::string & refChr,
                          const std::string & queryChr, const std::string &refSeq, const std::string &querySeq, std::ofstream & ofile2,
                          std::map <std::string, Fasta> & refSequences, std::map <std::string, Fasta> & targetSequences ){

    for( int32_t i=0; i<pairedSimilarFragments0.size(); ++i ){
        if( strand == POSITIVE ){

            int32_t quePos = startQuery + pairedSimilarFragments0[i].getStart2() - 1;
            int32_t refPos = startRef + pairedSimilarFragments0[i].getStart1() - 1;
            int32_t refStart = pairedSimilarFragments0[i].getStart1();
            int32_t queryStart = pairedSimilarFragments0[i].getStart2();

            std::stringstream refAlign;
            std::stringstream queryAlign;

            std::vector<uint32_t> newCigar;
            for( int32_t j=0; j<pairedSimilarFragments0[i].getCigar().size(); ++j ){
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j]>>4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j]&0xf;
                if (  cigarType == 0 ){
                    refAlign <<  refSeq.substr(refStart-1, cigarLength);
                    queryAlign <<  querySeq.substr(queryStart-1, cigarLength);
                    refStart += cigarLength;
                    queryStart += cigarLength;
                }else if( cigarType == 1 ){
                    queryAlign <<  querySeq.substr(queryStart-1, cigarLength);
                    queryStart += cigarLength;
                    for ( int repeatI = 0; repeatI<cigarLength; ++repeatI ) {
                        refAlign << "-";
                    }
                }else if( cigarType == 2 ){
                    refAlign <<  refSeq.substr(refStart-1, cigarLength);
                    refStart += cigarLength;
                    for ( int repeatI = 0; repeatI<cigarLength; ++repeatI ) {
                        queryAlign << "-";
                    }
                }else{
                    std::cout << "denovogenomevariancalling cpp line 50 is not working properly, please update the source codce" << std::endl;
                }
            }

            ofile2 << "a\tscore=" << pairedSimilarFragments0[i].getScore() << std::endl;

            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
//            temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
            ofile2 << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << refPos-1 << "\t" << std::setw(9) << temp.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
//            temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
            ofile2 << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << quePos-1 << "\t" << std::setw(9) << temp.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  queryAlign.str() << std::endl;
            ofile2 << std::endl;

        } else {

            int32_t quePos = startQuery + querySeq.size() - 1 - pairedSimilarFragments0[i].getEnd2() + 1; // this is the position on the positive strand
            int32_t refPos = startRef + pairedSimilarFragments0[i].getStart1() - 1;
            int32_t refStart = pairedSimilarFragments0[i].getStart1();
            int32_t queryStart = pairedSimilarFragments0[i].getStart2();

            std::stringstream refAlign;
            std::stringstream queryAlign;

            std::vector<uint32_t> newCigar;
            for (int32_t j = 0; j < pairedSimilarFragments0[i].getCigar().size(); ++j) {
                uint32_t cigarLength = pairedSimilarFragments0[i].getCigar()[j] >> 4;
                uint32_t cigarType = pairedSimilarFragments0[i].getCigar()[j] & 0xf;
                if (cigarType == 0) {
                    refAlign << refSeq.substr(refStart - 1, cigarLength);
                    queryAlign << querySeq.substr(queryStart - 1, cigarLength);
                    refStart += cigarLength;
                    queryStart += cigarLength;
                } else if (cigarType == 1) {
                    queryAlign << querySeq.substr(queryStart - 1, cigarLength);
                    queryStart += cigarLength;
                    for (int repeatI = 0; repeatI < cigarLength; ++repeatI) {
                        refAlign << "-";
                    }
                } else if (cigarType == 2) {
                    refAlign << refSeq.substr(refStart - 1, cigarLength);
                    refStart += cigarLength;
                    for (int repeatI = 0; repeatI < cigarLength; ++repeatI) {
                        queryAlign << "-";
                    }
                } else {
                    std::cout
                            << "denovogenomevariancalling cpp line 50 is not working properly, please update the source codce"
                            << std::endl;
                }
            }

            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
//            temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
            ofile2 << "a\tscore=" << pairedSimilarFragments0[i].getScore() << std::endl;
            ofile2 << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right << std::setw(9) << refPos - 1 << "\t" << std::setw(9) << temp.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
//            temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
            ofile2 << "s\t" << std::left << std::setw(chrWidth) << queryFileName + "." + queryChr << "\t" << std::right << std::setw(9) << quePos - 1 << "\t" << std::setw(9) << temp.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" << queryAlign.str() << std::endl;
            ofile2 << std::endl;
        }
    }
}


void genomeAlignment( std::vector<std::vector<AlignmentMatch>> & alignmentMatchsMap,
                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                      const int32_t & widownWidth, const int32_t & wfaSize, const int32_t & wfaSize2,
                      const std::string & outPutMafFile, const std::string & outPutFragedFile, /*std::string & outPutLocalalignmentFile,*/
                      const int32_t & matchingScore, const  int32_t & mismatchingPenalty,
                      const int32_t & openGapPenalty1, const int32_t & extendGapPenalty1,
                      const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                      int32_t & seed_window_size, const int32_t & mini_cns_score, const int32_t & step_size,
                      const int32_t & matrix_boundary_distance, const  int32_t & scoreThreshold, const  int32_t & w, const  int32_t & xDrop,
                      const int32_t & min_wavefront_length, const int32_t & max_distance_threshold) {

    Scorei m(matchingScore, mismatchingPenalty);

    affine2p_penalties_t penalties = {
            .match = 0,
            .mismatch = -mismatchingPenalty+matchingScore,
            .gap_opening1 = -openGapPenalty1+matchingScore,
            .gap_extension1 = -extendGapPenalty1+matchingScore,
            .gap_opening2 = -openGapPenalty2+matchingScore,
            .gap_extension2 = -extendGapPenalty2+matchingScore,
    };
    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);


    bool outPutMaf = false;
    bool outPutFraged = false;
//    bool outPutLocalalignment = false;

    if ( outPutMafFile.size() > 0 ){
        outPutMaf = true;
    }
    if ( outPutFragedFile.size() > 0 ){
        outPutFraged = true;
    }
//    if ( outPutLocalalignmentFile.size() > 0 ){
//        outPutLocalalignment = true;
//    }

    std::map <std::string, Fasta> refSequences;
    readFastaFile(refFastaFilePath, refSequences);

    std::map <std::string, Fasta> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);

    int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    queryFileName = elems.back();

    for ( std::map <std::string, Fasta>::iterator itchr  =refSequences.begin(); itchr != refSequences.end(); ++itchr ){
        if ( (refFileName + "." + itchr->first).size() > chrWidth ){
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }
    for ( std::map <std::string, Fasta>::iterator itchr  =targetSequences.begin(); itchr != targetSequences.end(); ++itchr ){
        if ((queryFileName + "." + itchr->first).size() > chrWidth){
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ofragfile;
    std::ofstream  oLocalalignment;
    if( outPutMaf ){
        omaffile.open(outPutMafFile);
        omaffile << "##maf version=1" << std::endl;
    }

    if( outPutFraged ){
        ofragfile.open(outPutFragedFile);
        ofragfile << "##maf version=1" << std::endl;
    }
//    if(outPutLocalalignment){
//        oLocalalignment.open(outPutLocalalignmentFile);
//        oLocalalignment<< "##maf version=1" << std::endl;
//    }

    SequenceCharToUInt8 sequenceCharToUInt8;
    int32_t sizei = alignmentMatchsMap.size();
    bool checkResult = false;

    for(int32_t i=sizei-1; i>=0; --i){
        STRAND strand = alignmentMatchsMap[i][0].getStrand();

        if(  POSITIVE == strand ){
            size_t startRef = alignmentMatchsMap[i][0].getRefStartPos();
            size_t startQuery;
            size_t endRef;
            size_t endQuery;
            std::stringstream refAlign;
            std::stringstream queryAlign;
            std::string refChr = alignmentMatchsMap[i][0].getRefChr();
            std::string queryChr = alignmentMatchsMap[i][0].getQueryChr();
            int64_t alignmentScore = 0;
            startQuery = alignmentMatchsMap[i][0].getQueryStartPos();
            for( AlignmentMatch orthologPair : alignmentMatchsMap[i] ){

                if (  orthologPair.getRefStartPos()==startRef && orthologPair.getQueryStartPos()!=startQuery ) {
                    endQuery = orthologPair.getQueryStartPos()-1;
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery);
                    for ( int repeatI = 0; repeatI<querySeq.length(); ++repeatI ){
                        refAlign<<"-";
                    }
                    queryAlign << querySeq;
                    alignmentScore += openGapPenalty1 + extendGapPenalty1 * querySeq.length();
                }else if ( orthologPair.getRefStartPos()!=startRef && orthologPair.getQueryStartPos()==startQuery ){
                    endRef = orthologPair.getRefStartPos()-1;
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    refAlign<<refSeq;
                    for ( int repeatI = 0; repeatI<refSeq.length(); ++repeatI ) {
                        queryAlign << "-";
                    }
                    alignmentScore += openGapPenalty1 + extendGapPenalty1*refSeq.length();
                }else if ( orthologPair.getRefStartPos()==startRef && orthologPair.getQueryStartPos()==startQuery ){

                }else{
                    endRef = orthologPair.getRefStartPos()-1;
                    endQuery = orthologPair.getQueryStartPos()-1;

                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery);

//                    if ( outPutLocalalignment ){
//
//                        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
//                        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
//                        int32_t length1 = refSeq.length();
//
//
//                        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
//                        int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
//                        int32_t length2 = querySeq.length();
//
//                        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
//                                                                                                                           seq2, seq2_rev_com, length1, length2, seed_window_size,
//                                                                                                                           mini_cns_score, matrix_boundary_distance, openGapPenalty1,
//                                                                                                                           extendGapPenalty1, matchingScore,
//                                                                                                                           mismatchingPenalty, m, step_size, refSeq, querySeq,
//                                                                                                                           scoreThreshold, w, xDrop);
//
//                        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
//                        pairedSimilarFragments0 = pairedSimilarFragments;
//
//                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);
//
//
//                    }
                    /*if(false)*/{

                        std::string _alignment_q;
                        std::string _alignment_d;
                        int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);

                        alignmentScore += thiScore;
                        refAlign << _alignment_d;
                        queryAlign << _alignment_q;
                        if (checkResult){
                            std::string temp;

                            temp = _alignment_d;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(refSeq) == 0);

                            temp = _alignment_q;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(querySeq)==0);
                        }
                        if (outPutFraged  ){
                            ofragfile << "a\tscore=" << thiScore << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofragfile << std::endl;
                        }
                    }
                }
                {
                    startRef = orthologPair.getRefStartPos();
                    startQuery= orthologPair.getQueryStartPos();
                    endRef = orthologPair.getRefEndPos();
                    endQuery = orthologPair.getQueryEndPos();
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery);
/*
                    if ( outPutLocalalignment ){

                        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                        int32_t length1 = refSeq.length();

                        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                        int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                        int32_t length2 = querySeq.length();

                        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                           seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                           mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                           extendGapPenalty1, matchingScore,
                                                                                                                           mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                           scoreThreshold, w, xDrop);

                        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                        pairedSimilarFragments0 = pairedSimilarFragments;

                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);

                    }*/
                    /*if(false)*/{
                        std::string _alignment_q;
                        std::string _alignment_d;

                        int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth,  wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1,openGapPenalty2, extendGapPenalty2,  min_wavefront_length, max_distance_threshold, m);

                        alignmentScore += thiScore;
                        refAlign << _alignment_d;
                        queryAlign << _alignment_q;
                        if (checkResult) {
                            std::string temp;
                            temp = _alignment_d;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(refSeq) == 0);

                            temp = _alignment_q;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(querySeq) == 0);
                        }
                        if (outPutFraged ){
                            ofragfile << "a\tscore=" << thiScore << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofragfile << std::endl;
                        }
                    }
                }
                startRef = orthologPair.getRefEndPos()+1;
                startQuery = orthologPair.getQueryEndPos()+1;
            }

            /*if(false)*/{
                std::string refGenomerSequence = getSubsequence( refSequences, refChr, alignmentMatchsMap[i][0].getRefStartPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getRefEndPos());
                std::string queryGenomerSequence = getSubsequence( targetSequences, queryChr, alignmentMatchsMap[i][0].getQueryStartPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getQueryEndPos());

                if (checkResult) {
                    std::string temp = refAlign.str();
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
    //                temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
                    assert(temp.compare(refGenomerSequence)==0);

                    temp = queryAlign.str();
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
    //                temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
                    assert(temp.compare(queryGenomerSequence)==0);
                }
                omaffile << "a\tscore=" << alignmentScore << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][0].getRefStartPos()-1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][0].getQueryStartPos()-1 << "\t" << std::setw(9) << queryGenomerSequence.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  queryAlign.str() << std::endl;
                omaffile << std::endl;
            }
        } else {
            size_t startRef = alignmentMatchsMap[i][0].getRefStartPos();
            size_t startQuery;
            size_t endRef;
            size_t endQuery = alignmentMatchsMap[i][0].getQueryEndPos();
            std::stringstream refAlign;
            std::stringstream queryAlign;
            std::string refChr = alignmentMatchsMap[i][0].getRefChr();
            std::string queryChr = alignmentMatchsMap[i][0].getQueryChr();

            int64_t alignmentScore = 0;
            for( AlignmentMatch orthologPair : alignmentMatchsMap[i] ){

                if (  orthologPair.getRefStartPos()==startRef && orthologPair.getQueryEndPos()!=endQuery ) {
                    startQuery = orthologPair.getQueryEndPos()+1;
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery, strand);
                    for ( int repeatI = 0; repeatI<querySeq.length(); ++repeatI ){
                        refAlign<<"-";
                    }
                    queryAlign<<querySeq;
                    alignmentScore += openGapPenalty1 + extendGapPenalty1 * querySeq.length();
  //                  std::cout << "line 395" << std::endl;
                }else if ( orthologPair.getRefStartPos()!=startRef && orthologPair.getQueryEndPos()==endQuery ){
                    endRef = orthologPair.getRefStartPos()-1;
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    refAlign<<refSeq;
                    for ( int repeatI = 0; repeatI<refSeq.length(); ++repeatI ) {
                        queryAlign << "-";
                    }
                    alignmentScore += openGapPenalty1 + extendGapPenalty1*refSeq.length();
//                    std::cout << "line 404" << std::endl;
                }else if ( orthologPair.getRefStartPos()==startRef && orthologPair.getQueryEndPos()==endQuery ){

                }else{
                    endRef = orthologPair.getRefStartPos()-1;
                    startQuery = orthologPair.getQueryEndPos()+1;
//                    std::cout << "line 410 " << refChr << "\t" << startRef << "\t" << endRef << ":\t:" << queryChr << "\t" << startQuery << "\t" << endQuery << std::endl;
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery, strand);
//                    std::cout << "line 413 " << refSeq << std::endl;
//                    std::cout << "line 414 " << querySeq << std::endl;
/*
                    if ( outPutLocalalignment ){

                        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                        int32_t length1 = refSeq.length();


                        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                        int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                        int32_t length2 = querySeq.length();

                        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                           seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                           mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                           extendGapPenalty1, matchingScore,
                                                                                                                           mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                           scoreThreshold, w, xDrop);

                        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                        pairedSimilarFragments0 = pairedSimilarFragments;

                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);

                    }*/
                    /*if(false)*/{

                        std::string _alignment_q;
                        std::string _alignment_d;

                        //int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &affine_penalties, mm_allocator, widownWidth, wfaSize, wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, min_wavefront_length, max_distance_threshold, m, true  );
                        int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth, wfaSize,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);

                        alignmentScore += thiScore;
                        refAlign << _alignment_d;
                        queryAlign << _alignment_q;

                        if (checkResult) {
                            std::string temp;
                            temp = _alignment_d;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(refSeq) == 0);
                            temp = _alignment_q;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(querySeq) == 0);
                        }
                        if (outPutFraged  ){
                            ofragfile << "a\tscore=" << thiScore << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofragfile << std::endl;
                        }
                    }
                }
                {
                    startRef = orthologPair.getRefStartPos();
                    endQuery= orthologPair.getQueryEndPos();
                    endRef = orthologPair.getRefEndPos();
                    startQuery = orthologPair.getQueryStartPos();
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery, strand);
/*
                    if ( outPutLocalalignment ){

                        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                        int32_t length1 = refSeq.length();

                        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                        int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                        int32_t length2 = querySeq.length();

                        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                           seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                           mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                           extendGapPenalty1, matchingScore,
                                                                                                                           mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                           scoreThreshold, w, xDrop);

                        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                        pairedSimilarFragments0 = pairedSimilarFragments;

                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);
                    }*/
                    /*if(false)*/{
                        std::string _alignment_q;
                        std::string _alignment_d;

                        int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth, wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);

                        alignmentScore += thiScore;
                        refAlign << _alignment_d;
                        queryAlign << _alignment_q;
                        if (checkResult) {
                            std::string temp;
                            temp = _alignment_d;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(refSeq) == 0);

                            temp = _alignment_q;
                            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                            assert(temp.compare(querySeq) == 0);
                        }
                        if (outPutFraged ){
                            ofragfile << "a\tscore=" << thiScore << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofragfile << std::endl;
                        }
                    }
                }
                startRef = orthologPair.getRefEndPos()+1;
                endQuery= orthologPair.getQueryStartPos()-1;
            }
//            std::cout << refChr << " last align" << std::endl;


            /*if(false)*/{
                std::string refGenomerSequence = getSubsequence( refSequences, refChr, alignmentMatchsMap[i][0].getRefStartPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getRefEndPos());
                std::string queryGenomerSequence = getSubsequence( targetSequences, queryChr, alignmentMatchsMap[i][0].getQueryEndPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getQueryStartPos(), strand);

                if (checkResult) {
                    std::string temp = refAlign.str();
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
    //                temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
                    assert(temp.compare(refGenomerSequence)==0);

                    temp = queryAlign.str();
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
    //                temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
                    assert(temp.compare(queryGenomerSequence)==0);
                 }

                omaffile << "a\tscore=" << alignmentScore << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][0].getRefStartPos()-1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getQueryStartPos()-1 << "\t" << std::setw(9) << queryGenomerSequence.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  queryAlign.str() << std::endl;
                omaffile << std::endl;
            }

        }
    }
    if( outPutMaf ){
        omaffile.close();
    }
    if( outPutFraged ){
        ofragfile.close();
    }
//    if(outPutLocalalignment){
//        oLocalalignment.close();
//    }
    mm_allocator_delete(mm_allocator);
}

void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::map <std::string, Fasta> & refSequences){
    alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, chr, refSequences, 0);
}

void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::map <std::string, Fasta> & refSequences, int32_t refLetterNumber) {
    std::vector<Variant> sdiRecordsThisOne;
    alignmentToVcf(queryAlignSeq, refAlignSeq, sdiRecordsThisOne, chr, refSequences, refLetterNumber);
    for (std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin(); itVariant != sdiRecordsThisOne.end(); ++itVariant) {
        ovcffile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChromosome() << "_" << itVariant->getPosition() <<"\t"+ itVariant->getReference() << "\t" << itVariant->getAlternative() << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" << std::endl;
    }
}
void alignmentToVcf(std::string & queryAlignSeq, std::string & refAlignSeq, std::vector<Variant> & sdiRecordsThisOne, std::string chr, std::map <std::string, Fasta> & refSequences, int32_t refLetterNumber){

    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'U', 'T');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'R', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'Y', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'S', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'W', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'K', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'M', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'B', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'D', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'H', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'V', 'N');

    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'U', 'T');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'R', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'Y', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'S', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'W', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'K', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'M', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'B', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'D', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'H', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'V', 'N');

    FirstLastList sdiRecords;

    for (int ai = 0; ai < refAlignSeq.length(); ai++) {
        if (refAlignSeq[ai] != '-') {
            ++refLetterNumber;
        }
        if (refAlignSeq[ai] != queryAlignSeq[ai]) {
            if (queryAlignSeq[ai] == '-') {
                if (NULL != sdiRecords.getLast() &&
                    sdiRecords.getLast()->getMapSingleRecord().getPosition() - sdiRecords.getLast()->getMapSingleRecord().getChanginglength() == refLetterNumber
                    && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() < 0 &&
                    sdiRecords.getLast()->getMapSingleRecord().getAlternative().compare("-") == 0) {
                    int position = sdiRecords.getLast()->getMapSingleRecord().getPosition();
                    std::string ori = sdiRecords.getLast()->getMapSingleRecord().getReference() + refAlignSeq[ai];
                    std::string result = "-";
                    Variant mapSingleRecord = Variant(chr, position, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.deleteLast();
                    sdiRecords.insertLast(data);

                }else{
                    std::string ori(1, refAlignSeq[ai]);
                    std::string result = "-";
                    Variant mapSingleRecord = Variant(chr, refLetterNumber, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                }
            } else if (refAlignSeq[ai] == '-') {
                if (sdiRecords.getLast() == NULL) {
                    int position = refLetterNumber + 1;
                    std::string ori = "-";
                    std::string result(1, queryAlignSeq[ai]);
                    Variant mapSingleRecord = Variant(chr, position, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                } else {
                    if (NULL != sdiRecords.getLast() &&
                        sdiRecords.getLast()->getMapSingleRecord().getPosition() ==
                        (refLetterNumber + 1)
                        && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                        sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                        int position = refLetterNumber + 1;
                        std::string ori = "-";
                        std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + queryAlignSeq[ai];
                        Variant mapSingleRecord = Variant(chr, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.deleteLast();
                        sdiRecords.insertLast(data);
                    } else {
                        int position = refLetterNumber + 1;
                        std::string ori = "-";
                        std::string result(1, queryAlignSeq[ai]);
                        Variant mapSingleRecord = Variant(chr, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
                    }
                }
            } else {
                int position = refLetterNumber;
                std::string ori(1, refAlignSeq[ai]);
                std::string result(1, queryAlignSeq[ai]);
                Variant mapSingleRecord = Variant(chr, position, ori, result);
                Data *data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            }
        }
    }

    for (int runingCound = 0; runingCound < 2; ++runingCound) {
        if ((sdiRecords.getFirst() != NULL)
            && (sdiRecords.getFirst()->getNext() != NULL)) {
            Data *prevOne = sdiRecords.getFirst();
            Data *currOne = (sdiRecords.getFirst()->getNext());
            while ((currOne != NULL) && (NULL != currOne->getNext())) {
//                std::cout << "725" << std::endl;
                if (sdiRecords.getFirst() == currOne) {
                    prevOne = currOne;
                    currOne = prevOne->getNext();
                }
                if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                    prevOne->getMapSingleRecord().getChanginglength() < 0 &&
                    (prevOne->getMapSingleRecord().getPosition() + abs(prevOne->getMapSingleRecord().getChanginglength())) ==
                    currOne->getMapSingleRecord().getPosition()
                    && prevOne->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                    currOne->getMapSingleRecord().getAlternative().compare("-") == 0) { // merge two deletions

                    int position = prevOne->getMapSingleRecord().getPosition();
                    std::string ori = prevOne->getMapSingleRecord().getReference() + currOne->getMapSingleRecord().getReference();
                    std::string result = "-";
                    Variant mapSingleRecord2(chr, position, ori, result);
                    //delete prev one begin
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                    } else {
                        prevOne->getPrev()->setNext(currOne);
                        currOne->setPrev(prevOne->getPrev());
                        delete (prevOne);
                    } //delete prev one end
                    currOne->setMapSingleRecord(mapSingleRecord2);
                    prevOne = currOne->getPrev();
                    if (prevOne == NULL) {
                        prevOne = currOne;
                        currOne = prevOne->getNext();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() == 0 &&
                           0 == currOne->getMapSingleRecord().getReference().compare(currOne->getMapSingleRecord().getAlternative())) { // nonsense records
                    //delete current one
                    //std::cout << "820 delete prev" << std::endl;
                    prevOne->setNext(currOne->getNext());
                    currOne->getNext()->setPrev(prevOne);
                    delete (currOne);
                    currOne = prevOne->getNext();
                } else if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                           prevOne->getMapSingleRecord().getChanginglength() > 0
                           && currOne->getMapSingleRecord().getReference().compare(
                        prevOne->getMapSingleRecord().getAlternative()) == 0 &&
                           currOne->getMapSingleRecord().getPosition() ==
                           prevOne->getMapSingleRecord().getPosition()) { //delete one insertion and next reverse sence deletion
                    //delete current one and prev
                    //std::cout << "830 delete current one and prev" << std::endl;
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                        sdiRecords.deleteFirst();
                        prevOne = sdiRecords.getFirst();
                        currOne = (sdiRecords.getFirst()->getNext());
                    } else if (currOne == sdiRecords.getLast()) {
                        sdiRecords.deleteLast();
                        sdiRecords.deleteLast();
                        currOne = sdiRecords.getLast();
                        prevOne = currOne->getPrev();
                    } else {
                        currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                        currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                        Data *temp = currOne->getNext();
                        delete (currOne->getPrev());
                        delete (currOne);
                        currOne = temp;
                        prevOne = temp->getPrev();
                    }
                } else if (currOne->getMapSingleRecord().getChanginglength() > 0 &&
                           prevOne->getMapSingleRecord().getChanginglength() < 0
                           && currOne->getMapSingleRecord().getAlternative().compare(
                        prevOne->getMapSingleRecord().getReference()) == 0 &&
                           (currOne->getMapSingleRecord().getPosition() - 1) ==
                           prevOne->getMapSingleRecord().getPosition()) {
                    //delete current one and prev
                    //std::cout << "850 delete current one and prev" << std::endl;
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                        sdiRecords.deleteFirst();
                        prevOne = sdiRecords.getFirst();
                        currOne = (sdiRecords.getFirst()->getNext());
                    } else {
                        currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                        currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                        Data *temp = currOne->getNext();
                        delete (currOne->getPrev());
                        delete (currOne);
                        currOne = temp;
                        prevOne = temp->getPrev();
                    }
                } else {
                    prevOne = currOne;
                    currOne = prevOne->getNext();
                }//std::cout <<  (*itName) << ": link data structure end " << currOne->getMapSingleRecord().getPosition() << std::endl;
//                std::cout << "820" << std::endl;
            }
        }
    }

    if (sdiRecords.getFirst() != NULL) {
        Data *thisone = sdiRecords.getFirst();
        while (thisone != NULL) {
            sdiRecordsThisOne.push_back(thisone->getMapSingleRecord());
            thisone = (thisone->getNext());
        }
    }

    // clear RAM assigned by new Data() begin
    if (sdiRecords.getFirst() != NULL) {
        Data *thisone = sdiRecords.getFirst();
        while (thisone != NULL) {
            Data *tempData = thisone;
            thisone = (thisone->getNext());
            delete (tempData);
        }
    }// clear RAM assigned by new Data() end

    // transform link to vector and sort and merge nearby records begin
    bool ifChanged = true;
    while (ifChanged) { // link deletions together
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                sdiRecordsThisOne[j - 1].getChanginglength() < 0 &&
                (sdiRecordsThisOne[j - 1].getPosition() + abs(sdiRecordsThisOne[j - 1].getChanginglength())) == sdiRecordsThisOne[j].getPosition() &&
                sdiRecordsThisOne[j - 1].getAlternative().compare("-") == 0 &&
                sdiRecordsThisOne[j].getAlternative().compare("-") == 0) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = "-";
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }
    // transform link to vector and sort and merge nearby records end

    ifChanged = true;
    while (ifChanged && sdiRecordsThisOne.size() > 1) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j-1].getChanginglength() < 0 && sdiRecordsThisOne[j-1].getAlternative().compare("-") == 0 &&
                sdiRecordsThisOne[j].getChanginglength() > 0 && sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                ( (sdiRecordsThisOne[j-1].getPosition() + sdiRecordsThisOne[j-1].getReference().size()) == sdiRecordsThisOne[j].getPosition()) ){
                int position = sdiRecordsThisOne[j - 1].getPosition();

                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            }else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }

    // merge nearby indels begin
    ifChanged = true;
    while (ifChanged && sdiRecordsThisOne.size() > 1) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j].getChanginglength() < 0 &&  sdiRecordsThisOne[j].getAlternative().compare("-") == 0 &&
                ( //(sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition()-1)  && sdiRecordsThisOne[j - 1].getReference() != "-" ) ||
                (sdiRecordsThisOne[j - 1].getReference()[0]!='-' &&
                sdiRecordsThisOne[j - 1].getPosition() + sdiRecordsThisOne[j - 1].getReference().size() == sdiRecordsThisOne[j].getPosition() ) ) ){
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            }else if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                      sdiRecordsThisOne[j].getAlternative().compare("-") == 0 && (
                              (sdiRecordsThisOne[j - 1].getReference() != "-" && ( sdiRecordsThisOne[j - 1].getPosition() + (sdiRecordsThisOne[j - 1].getReference().size()) ) < (sdiRecordsThisOne[j].getPosition() ))
                        ||  (sdiRecordsThisOne[j - 1].getReference() == "-" &&  sdiRecordsThisOne[j - 1].getPosition() != (sdiRecordsThisOne[j].getPosition() ))    )
                      ){
                int position = sdiRecordsThisOne[j].getPosition()-1;
                std::string ori(1, refSequences[chr].getSequence()[position-1]);
                std::string result = ori;
                ori = ori + sdiRecordsThisOne[j].getReference();
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                       sdiRecordsThisOne[j].getAlternative().compare("-") == 0 &&  (sdiRecordsThisOne[j - 1].getReference() == "-" &&  sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition() )    )
                    ){
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                       sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                    sdiRecordsThisOne[j-1].getReference().compare("-") != 0 &&
                    ( sdiRecordsThisOne[j - 1].getPosition() + sdiRecordsThisOne[j - 1].getReference().size() ) == (sdiRecordsThisOne[j].getPosition())  ){
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();
                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            }else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                      sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                    sdiRecordsThisOne[j - 1].getReference() != "-" &&
                    ( sdiRecordsThisOne[j - 1].getPosition() + (sdiRecordsThisOne[j - 1].getReference().size()) ) < (sdiRecordsThisOne[j].getPosition()) ){
                int position = sdiRecordsThisOne[j].getPosition()-1;
                std::string ori(1, refSequences[chr].getSequence()[position-1]);
                std::string result = ori + sdiRecordsThisOne[j].getAlternative();
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }

    ifChanged = true;
    while (ifChanged && sdiRecordsThisOne.size() > 1) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j - 1].getReference().compare("-") != 0 &&
                sdiRecordsThisOne[j].getChanginglength() > 0 && sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                sdiRecordsThisOne[j - 1].getPosition() == sdiRecordsThisOne[j].getPosition()) {
                int position = sdiRecordsThisOne[j - 1].getPosition();

                std::string ori = sdiRecordsThisOne[j].getReference() + sdiRecordsThisOne[j - 1].getReference();
                std::string result = sdiRecordsThisOne[j].getAlternative() + sdiRecordsThisOne[j - 1].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
        sdiRecordsThisOne.shrink_to_fit();
    }

    if (sdiRecordsThisOne.size() >0 && sdiRecordsThisOne[0].getChanginglength() > 0 &&
        sdiRecordsThisOne[0].getReference().compare("-") == 0){
        if( sdiRecordsThisOne[0].getPosition() != 1 ){
            int position = sdiRecordsThisOne[0].getPosition()-1;
            std::string ori(1, refSequences[chr].getSequence()[position-1]);
            std::string result = ori + sdiRecordsThisOne[0].getAlternative();
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else if ( sdiRecordsThisOne[0].getPosition() == 1 && sdiRecordsThisOne.size() >1 && sdiRecordsThisOne[1].getPosition() == 1  ) {
            int position = 1;
            std::string ori = sdiRecordsThisOne[0].getReference() + sdiRecordsThisOne[1].getReference();
            std::string result = sdiRecordsThisOne[0].getAlternative() + sdiRecordsThisOne[1].getAlternative();

            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + 1);
            sdiRecordsThisOne.shrink_to_fit();
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else{
            int position = 1;
            std::string ori(1, refSequences[chr].getSequence()[0]);
            std::string result = sdiRecordsThisOne[0].getAlternative() + ori;
            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }
    }

    if (sdiRecordsThisOne.size() >0 && sdiRecordsThisOne[0].getChanginglength() < 0 &&
        sdiRecordsThisOne[0].getAlternative().compare("-") == 0){
        if( sdiRecordsThisOne[0].getPosition() != 1 ){
            int position = sdiRecordsThisOne[0].getPosition()-1;
            std::string ori(1, refSequences[chr].getSequence()[position-1]);
            std::string result = ori;
            ori += sdiRecordsThisOne[0].getReference();
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else if ( sdiRecordsThisOne[0].getPosition() == 1 && sdiRecordsThisOne.size() >1 && sdiRecordsThisOne[1].getPosition() == sdiRecordsThisOne[0].getPosition()+sdiRecordsThisOne[0].getReference().size()  ) {
            int position = 1;
            std::string ori = sdiRecordsThisOne[0].getReference() + sdiRecordsThisOne[1].getReference();
            std::string result = sdiRecordsThisOne[0].getAlternative() + sdiRecordsThisOne[1].getAlternative();

            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + 1);
            sdiRecordsThisOne.shrink_to_fit();
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else{
            int position = 1;
            std::string result(1, refSequences[chr].getSequence()[sdiRecordsThisOne[0].getReference().length()]);
            std::string ori = sdiRecordsThisOne[0].getReference() + result;
            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }
    }
}




void alignmentToVcfBackup(std::string & queryAlignSeq, std::string & refAlignSeq, std::ofstream & ovcffile, std::string chr, std::map <std::string, Fasta> & refSequences, int32_t refLetterNumber){

    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'U', 'T');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'R', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'Y', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'S', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'W', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'K', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'M', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'B', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'D', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'H', 'N');
    std::replace(refAlignSeq.begin(), refAlignSeq.end(), 'V', 'N');

    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'U', 'T');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'R', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'Y', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'S', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'W', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'K', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'M', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'B', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'D', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'H', 'N');
    std::replace(queryAlignSeq.begin(), queryAlignSeq.end(), 'V', 'N');

    FirstLastList sdiRecords;

    for (int ai = 0; ai < refAlignSeq.length(); ai++) {
        if (refAlignSeq[ai] != '-') {
            ++refLetterNumber;
        }
        if (refAlignSeq[ai] != queryAlignSeq[ai]) {
            if (queryAlignSeq[ai] == '-') {
                std::string ori(1, refAlignSeq[ai]);
                std::string result = "-";
                Variant mapSingleRecord = Variant(chr, refLetterNumber, ori, result);
                Data *data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            } else if (refAlignSeq[ai] == '-') {
                if (sdiRecords.getLast() == NULL) {
                    int position = refLetterNumber + 1;
                    std::string ori = "-";
                    std::string result(1, queryAlignSeq[ai]);
                    Variant mapSingleRecord = Variant(chr, position, ori, result);
                    Data *data = new Data(mapSingleRecord);
                    sdiRecords.insertLast(data);
                } else {
                    if (NULL != sdiRecords.getLast() &&
                        sdiRecords.getLast()->getMapSingleRecord().getPosition() ==
                        (refLetterNumber + 1)
                        && sdiRecords.getLast()->getMapSingleRecord().getChanginglength() > 0 &&
                        sdiRecords.getLast()->getMapSingleRecord().getReference().compare("-") == 0) {

                        int position = refLetterNumber + 1;
                        std::string ori = "-";
                        std::string result = sdiRecords.getLast()->getMapSingleRecord().getAlternative() + queryAlignSeq[ai];
                        Variant mapSingleRecord = Variant(chr, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.deleteLast();
                        sdiRecords.insertLast(data);
                    } else {
                        int position = refLetterNumber + 1;
                        std::string ori = "-";
                        std::string result(1, queryAlignSeq[ai]);
                        Variant mapSingleRecord = Variant(chr, position, ori, result);
                        Data *data = new Data(mapSingleRecord);
                        sdiRecords.insertLast(data);
                    }
                }
            } else {
                int position = refLetterNumber;
                std::string ori(1, refAlignSeq[ai]);
                std::string result(1, queryAlignSeq[ai]);
                Variant mapSingleRecord = Variant(chr, position, ori, result);
                Data *data = new Data(mapSingleRecord);
                sdiRecords.insertLast(data);
            }
        }
    }

    for (int runingCound = 0; runingCound < 2; ++runingCound) {
        if ((sdiRecords.getFirst() != NULL)
            && (sdiRecords.getFirst()->getNext() != NULL)) {
            Data *prevOne = sdiRecords.getFirst();
            Data *currOne = (sdiRecords.getFirst()->getNext());
            while ((currOne != NULL) && (NULL != currOne->getNext())) {
                //                        std::cout << "2109" << std::endl;
                if (sdiRecords.getFirst() == currOne) {
                    prevOne = currOne;
                    currOne = prevOne->getNext();
                }
                if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                    prevOne->getMapSingleRecord().getChanginglength() < 0 &&
                    (prevOne->getMapSingleRecord().getPosition() +
                     abs(prevOne->getMapSingleRecord().getChanginglength())) ==
                    currOne->getMapSingleRecord().getPosition()
                    && prevOne->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                    currOne->getMapSingleRecord().getAlternative().compare("-") == 0) { // merge two deletions

                    int position = prevOne->getMapSingleRecord().getPosition();
                    std::string ori = prevOne->getMapSingleRecord().getReference() +
                                      currOne->getMapSingleRecord().getReference();
                    std::string result = "-";
                    Variant mapSingleRecord2(chr, position, ori, result);
                    //delete prev one begin
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                    } else {
                        prevOne->getPrev()->setNext(currOne);
                        currOne->setPrev(prevOne->getPrev());
                        delete (prevOne);
                    } //delete prev one end
                    currOne->setMapSingleRecord(mapSingleRecord2);
                    prevOne = currOne->getPrev();
                    if (prevOne == NULL) {
                        prevOne = currOne;
                        currOne = prevOne->getNext();
                    }
                    //std::cout << "817 delete prev" << std::endl;
                } else if (currOne->getMapSingleRecord().getChanginglength() == 0 &&
                           0 == currOne->getMapSingleRecord().getReference().compare(currOne->getMapSingleRecord().getAlternative())) { // nonsense records
                    //delete current one
                    //std::cout << "820 delete prev" << std::endl;
                    prevOne->setNext(currOne->getNext());
                    currOne->getNext()->setPrev(prevOne);
                    delete (currOne);
                    currOne = prevOne->getNext();
                    //std::cout << "825 delete prev" << std::endl;
                } else if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                           prevOne->getMapSingleRecord().getChanginglength() > 0
                           && currOne->getMapSingleRecord().getReference().compare(
                        prevOne->getMapSingleRecord().getAlternative()) == 0 &&
                           currOne->getMapSingleRecord().getPosition() ==
                           prevOne->getMapSingleRecord().getPosition()) { //delete one insertion and next reverse sence deletion
                    //delete current one and prev
                    //std::cout << "830 delete current one and prev" << std::endl;
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                        sdiRecords.deleteFirst();
                        prevOne = sdiRecords.getFirst();
                        currOne = (sdiRecords.getFirst()->getNext());
                    } else if (currOne == sdiRecords.getLast()) {
                        sdiRecords.deleteLast();
                        sdiRecords.deleteLast();
                        currOne = sdiRecords.getLast();
                        prevOne = currOne->getPrev();
                    } else {
                        currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                        currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                        Data *temp = currOne->getNext();
                        delete (currOne->getPrev());
                        delete (currOne);
                        currOne = temp;
                        prevOne = temp->getPrev();
                    }
                    //std::cout << "844 delete current one and prev" << std::endl;
                } else if (currOne->getMapSingleRecord().getChanginglength() > 0 &&
                           prevOne->getMapSingleRecord().getChanginglength() < 0
                           && currOne->getMapSingleRecord().getAlternative().compare(
                        prevOne->getMapSingleRecord().getReference()) == 0 &&
                           (currOne->getMapSingleRecord().getPosition() - 1) ==
                           prevOne->getMapSingleRecord().getPosition()) {
                    //delete current one and prev
                    //std::cout << "850 delete current one and prev" << std::endl;
                    if (((currOne->getPrev())) == (sdiRecords.getFirst())) {
                        sdiRecords.deleteFirst();
                        sdiRecords.deleteFirst();
                        prevOne = sdiRecords.getFirst();
                        currOne = (sdiRecords.getFirst()->getNext());
                    } else {
                        currOne->getPrev()->getPrev()->setNext(currOne->getNext());
                        currOne->getNext()->setPrev(currOne->getPrev()->getPrev());
                        Data *temp = currOne->getNext();
                        delete (currOne->getPrev());
                        delete (currOne);
                        currOne = temp;
                        prevOne = temp->getPrev();
                    }
                    //                    std::cout << "865 delete current one and prev" << std::endl;
                } else {
                    prevOne = currOne;
                    currOne = prevOne->getNext();
                }//std::cout <<  (*itName) << ": link data structure end " << currOne->getMapSingleRecord().getPosition() << std::endl;
                //std::cout << "2009" << std::endl;
            }
        }
    }
    //end: merge link data structure

    //        std::cout << it0->first << " link data structure end" << std::endl;
    std::vector<Variant> sdiRecordsThisOne;
    if (sdiRecords.getFirst() != NULL) {
        Data *thisone = sdiRecords.getFirst();
        while (thisone != NULL) {
            sdiRecordsThisOne.push_back(thisone->getMapSingleRecord());
            thisone = (thisone->getNext());
        }
    }

    // clear RAM assigned by new Data() begin
    if (sdiRecords.getFirst() != NULL) {
        Data *thisone = sdiRecords.getFirst();
        while (thisone != NULL) {
            Data *tempData = thisone;
            thisone = (thisone->getNext());
            delete (tempData);
        }
    }// clear RAM assigned by new Data() end

//                std::cout << " transform link to vector and sort and merge nearby records begin" << std::endl;
    // transform link to vector and sort and merge nearby records begin
    bool ifChanged = true;
    while (ifChanged) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                sdiRecordsThisOne[j - 1].getChanginglength() < 0 &&
                (sdiRecordsThisOne[j - 1].getPosition() +
                 abs(sdiRecordsThisOne[j - 1].getChanginglength())) ==
                sdiRecordsThisOne[j].getPosition()
                && sdiRecordsThisOne[j - 1].getAlternative().compare("-") == 0 &&
                sdiRecordsThisOne[j].getAlternative().compare("-") == 0) {
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = "-";
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) ==
                       0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1;
             intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
        }
    }
    // transform link to vector and sort and merge nearby records end
//                std::cout << "transform link to vector and sort and merge nearby records end" << std::endl;


    ifChanged = true;
    while (ifChanged) {
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j-1].getChanginglength() < 0 &&
                sdiRecordsThisOne[j-1].getAlternative().compare("-") == 0 &&
                sdiRecordsThisOne[j].getChanginglength() >0 && sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                ( (sdiRecordsThisOne[j-1].getPosition() + sdiRecordsThisOne[j-1].getReference().size()) == sdiRecordsThisOne[j].getPosition()) ){

                int position = sdiRecordsThisOne[j - 1].getPosition();

                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            }else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }

        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
            sdiRecordsThisOne.shrink_to_fit();
        }
//                    std::cout << "merge nearby indels end 1078" << std::endl;
    }



//                std::cout << "merge nearby indels begin" << std::endl;
    // merge nearby indels begin
    ifChanged = true;
    while (ifChanged) {
//                    std::cout << "merge nearby indels begin 1013" << std::endl;
        std::sort(sdiRecordsThisOne.begin(), sdiRecordsThisOne.end());
        ifChanged = false;
        std::vector<int> sdiRecordsToRomove;
        int oldSize = sdiRecordsThisOne.size();
        for (int j = 1; j < oldSize; j++) {
            if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                sdiRecordsThisOne[j].getAlternative().compare("-") == 0 &&
                (sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition()-1) ||
                 (sdiRecordsThisOne[j - 1].getReference()[0]!='-' && sdiRecordsThisOne[j - 1].getPosition() + sdiRecordsThisOne[j - 1].getReference().size() == sdiRecordsThisOne[j].getPosition() ) ) ){
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            }else if (sdiRecordsThisOne[j].getChanginglength() < 0 &&
                      sdiRecordsThisOne[j].getAlternative().compare("-") == 0 &&
                      sdiRecordsThisOne[j - 1].getPosition() != (sdiRecordsThisOne[j].getPosition()-1) ){

                int position = sdiRecordsThisOne[j].getPosition()-1;
                std::string ori(1, refSequences[chr].getSequence()[position-1]);
                std::string result = ori;
                ori += sdiRecordsThisOne[j].getReference();
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                       sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                       sdiRecordsThisOne[j - 1].getPosition() == (sdiRecordsThisOne[j].getPosition()-1) ){
                int position = sdiRecordsThisOne[j - 1].getPosition();
                std::string ori = sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                std::string result = sdiRecordsThisOne[j - 1].getAlternative() + sdiRecordsThisOne[j].getAlternative();

                ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
                result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsToRomove.push_back(j - 1);
                sdiRecordsThisOne[j] = mapSingleRecord2;
                j++;
                ifChanged = true;
            }else if (sdiRecordsThisOne[j].getChanginglength() > 0 &&
                      sdiRecordsThisOne[j].getReference().compare("-") == 0 &&
                      sdiRecordsThisOne[j - 1].getPosition() != (sdiRecordsThisOne[j].getPosition()-1) ){
                int position = sdiRecordsThisOne[j].getPosition()-1;
                std::string ori(1, refSequences[chr].getSequence()[position-1]);
                std::string result = ori + sdiRecordsThisOne[j].getAlternative();
                Variant mapSingleRecord2(chr, position, ori, result);
                sdiRecordsThisOne[j] = mapSingleRecord2;
            } else if (sdiRecordsThisOne[j].getReference().compare(sdiRecordsThisOne[j].getAlternative()) == 0) {
                sdiRecordsToRomove.push_back(j); // it does not affect sorting
            }
        }
//                    std::cout << "merge nearby indels end 1048" << std::endl;
        for (int intTpRomoveIndex = sdiRecordsToRomove.size() - 1; intTpRomoveIndex >= 0; --intTpRomoveIndex) {
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + sdiRecordsToRomove[intTpRomoveIndex]);
            sdiRecordsThisOne.shrink_to_fit();
        }
//                    std::cout << "merge nearby indels end 1078" << std::endl;
    }
//                std::cout << "merge nearby indels end" << std::endl;

    if (sdiRecordsThisOne.size() >0 && sdiRecordsThisOne[0].getChanginglength() > 0 &&
        sdiRecordsThisOne[0].getReference().compare("-") == 0){
        if( sdiRecordsThisOne[0].getPosition() != 1 ){
            int position = sdiRecordsThisOne[0].getPosition()-1;
            std::string ori(1, refSequences[chr].getSequence()[position-1]);
            std::string result = ori + sdiRecordsThisOne[0].getAlternative();
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else if ( sdiRecordsThisOne[0].getPosition() == 1 && sdiRecordsThisOne.size() >1 && sdiRecordsThisOne[1].getPosition() == 1  ) {
            int position = 1;
            std::string ori = sdiRecordsThisOne[0].getReference() + sdiRecordsThisOne[1].getReference();
            std::string result = sdiRecordsThisOne[0].getAlternative() + sdiRecordsThisOne[1].getAlternative();

            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + 1);
            sdiRecordsThisOne.shrink_to_fit();
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else{
            int position = 1;
            std::string ori(1, refSequences[chr].getSequence()[0]);
            std::string result = sdiRecordsThisOne[0].getAlternative() + ori;
            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }
    }
    if (sdiRecordsThisOne.size() >0 && sdiRecordsThisOne[0].getChanginglength() < 0 &&
        sdiRecordsThisOne[0].getAlternative().compare("-") == 0){
        if( sdiRecordsThisOne[0].getPosition() != 1 ){
            int position = sdiRecordsThisOne[0].getPosition()-1;
            std::string ori(1, refSequences[chr].getSequence()[position-1]);
            std::string result = ori;
            ori += sdiRecordsThisOne[0].getReference();
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else if ( sdiRecordsThisOne[0].getPosition() == 1 && sdiRecordsThisOne.size() >1 && sdiRecordsThisOne[1].getPosition() == sdiRecordsThisOne[0].getPosition()+sdiRecordsThisOne[0].getReference().size()  ) {
            int position = 1;
            std::string ori = sdiRecordsThisOne[0].getReference() + sdiRecordsThisOne[1].getReference();
            std::string result = sdiRecordsThisOne[0].getAlternative() + sdiRecordsThisOne[1].getAlternative();

            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());

            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne.erase(sdiRecordsThisOne.begin() + 1);
            sdiRecordsThisOne.shrink_to_fit();
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }else{
            int position = 1;
            std::string result(1, refSequences[chr].getSequence()[sdiRecordsThisOne[0].getReference().length()]);
            std::string ori = sdiRecordsThisOne[0].getReference() + ori;
            ori.erase(std::remove(ori.begin(), ori.end(), '-'), ori.end());
            result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
            Variant mapSingleRecord2(chr, position, ori, result);
            sdiRecordsThisOne[0] = mapSingleRecord2;
        }
    }

    for (std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin();
         itVariant != sdiRecordsThisOne.end(); ++itVariant) {
        ovcffile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChromosome() << "_" << itVariant->getPosition() <<"\t"+
        itVariant->getReference() << "\t" << itVariant->getAlternative() << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" << std::endl;
    }
}


void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                                      const int32_t & widownWidth, const int32_t & wfaSize, const int32_t & wfaSize2, const std::string & outPutMafFile, const std::string & outPutVcfFile,
                                      const std::string & outPutFragedFile, /*std::string & outPutLocalalignmentFile, */const int32_t & matchingScore, const  int32_t & mismatchingPenalty,
                                      const  int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const  int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
                                      const int32_t & min_wavefront_length, const int32_t & max_distance_threshold, int32_t & seed_window_size, const int32_t & mini_cns_score, const int32_t & step_size,
                                      const int32_t & matrix_boundary_distance, const int32_t & scoreThreshold, const  int32_t & w, const  int32_t & xDrop) {

    Scorei m(matchingScore, mismatchingPenalty);

    affine2p_penalties_t penalties = {
            .match = 0,
            .mismatch = -mismatchingPenalty+matchingScore,
            .gap_opening1 = -openGapPenalty1+matchingScore,
            .gap_extension1 = -extendGapPenalty1+matchingScore,
            .gap_opening2 = -openGapPenalty2+matchingScore,
            .gap_extension2 = -extendGapPenalty2+matchingScore,
    };
    std::cout << "match:" << penalties.match << std::endl;
    std::cout << "mismatch:" << penalties.mismatch << std::endl;
    std::cout << "gap_opening1:" << penalties.gap_opening1 << std::endl;
    std::cout << "gap_extension1:" << penalties.gap_extension1 << std::endl;
    std::cout << "gap_opening2:" << penalties.gap_opening2 << std::endl;
    std::cout << "gap_extension2:" << penalties.gap_extension2 << std::endl;
    std::cout << "line 1501" << std::endl;
    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_512M);

    bool outPutMaf = false;
    bool outPutVcf = false;
    bool outPutFraged = false;
//    bool outPutLocalalignment = false;

    if ( outPutMafFile.size() > 0 ){
        outPutMaf = true;
    }
    if ( outPutVcfFile.size() > 0 ){
        outPutVcf = true;
    }
    if ( outPutFragedFile.size() > 0 ){
        outPutFraged = true;
    }
//    if ( outPutLocalalignmentFile.size() > 0 ){
//        outPutLocalalignment = true;
//    }

    std::map <std::string, Fasta> refSequences;
    readFastaFile(refFastaFilePath, refSequences);

    std::map <std::string, Fasta> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);

    int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    queryFileName = elems.back();

    for ( std::map <std::string, Fasta>::iterator itchr  =refSequences.begin(); itchr != refSequences.end(); ++itchr ){
        if ( (refFileName + "." + itchr->first).size() > chrWidth ){
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }
    for ( std::map <std::string, Fasta>::iterator itchr  =targetSequences.begin(); itchr != targetSequences.end(); ++itchr ){
        if ((queryFileName + "." + itchr->first).size() > chrWidth){
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ovcffile;
    std::ofstream ofragfile;
    std::ofstream  oLocalalignment;
    if( outPutMaf ){
        omaffile.open(outPutMafFile);
        omaffile << "##maf version=1" << std::endl;
    }

    SequenceCharToUInt8 sequenceCharToUInt8;

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon));
    if( ltm->tm_mday < 10 ){
        filedate = filedate + "0" + std::to_string((ltm->tm_mday));
    }else{
        filedate = filedate + std::to_string((ltm->tm_mday));
    }

    if( outPutVcf ){
        ovcffile.open(outPutVcfFile);
        ovcffile << "##fileformat=VCFv4.3" << std::endl;
        ovcffile << "##fileDate=" << filedate << std::endl;
        ovcffile << "##source=proali" << std::endl;
        ovcffile <<"##reference=" << refFileName << std::endl;
        ovcffile <<"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
        ovcffile <<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
        ovcffile <<"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
        ovcffile <<"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
        ovcffile <<"##FILTER=<ID=q30,Description=\"Quality below 30\">" << std::endl;
        std::string accession = queryFileName;
        accession = songStrReplaceAll(accession, ".fasta", "");
        accession = songStrReplaceAll(accession, ".fa", "");
        accession.erase(std::remove(accession.begin(), accession.end(), ' '), accession.end());
        ovcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << accession << std::endl;
    }

    if( outPutFraged ){
        ofragfile.open(outPutFragedFile);
        ofragfile << "##maf version=1" << std::endl;
    }
//    if(outPutLocalalignment){
//        oLocalalignment.open(outPutLocalalignmentFile);
//        oLocalalignment<< "##maf version=1" << std::endl;
//    }

    for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it0=alignmentMatchsMap.begin(); it0 != alignmentMatchsMap.end(); ++it0 ) {
        if( it0->second.size() > 0 ){
            size_t startRef = 1;
            size_t startQuery = 1;
            size_t endRef;
            size_t endQuery;
            std::stringstream refAlign;
            std::stringstream queryAlign;

            std::string refChr = it0->second[0].getRefChr();
            std::string queryChr = it0->second[0].getQueryChr();
            std::cout << refChr << " align start" << std::endl;

            int64_t alignmentScore = 0;
            STRAND lastStrand = POSITIVE;
            AlignmentMatch lastAlignmentMatch;

            int mafRefStart = 0;
            int mafQueryStart = 0;
            std::string mafStrand = "+";
            bool hasInversion = false;

            for (AlignmentMatch alignmentMatch : it0->second) {
                if( alignmentMatch.getStrand() == NEGATIVE ){
                    hasInversion = true;
                }

                if( lastStrand == POSITIVE && alignmentMatch.getStrand() == POSITIVE ){
//                    std::cout << "line 1625" << std::endl;
                    if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() != startQuery) {
                        endQuery = alignmentMatch.getQueryStartPos() - 1;
                        std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);
                        for (int repeatI = 0; repeatI < querySeq.length(); ++repeatI) {
                            refAlign << "-";
                        }
                        queryAlign << querySeq;
                        alignmentScore -= openGapPenalty1 + extendGapPenalty1*querySeq.length();
                    } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryStartPos() == startQuery) {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                        refAlign << refSeq;
//                        std::cout << "line 621 " <<  startRef << "\t" << endRef << std::endl;
                        for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                            queryAlign << "-";
                        }
                        alignmentScore -= openGapPenalty1 + extendGapPenalty1*refSeq.length();
                    } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() == startQuery) {

                    } else {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        endQuery = alignmentMatch.getQueryStartPos() - 1;

                        std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                        std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);

                        std::string _alignment_q;
                        std::string _alignment_d;

                        int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth, wfaSize, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);
                        std::string temp;

                        temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq) == 0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq) == 0);

                        alignmentScore += thiScore;
                        refAlign << _alignment_d;
                        queryAlign << _alignment_q;

                        //assert(temp.compare(querySeq) == 0);

                        if (outPutFraged  ){
                            ofragfile << "a\tscore=" << thiScore << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofragfile << std::endl;
                        }
/*
                        if ( outPutLocalalignment ){

                            int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                            int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                            int32_t length1 = refSeq.length();


                            int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                            int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                            int32_t length2 = querySeq.length();

                            std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                               seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                               mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                               extendGapPenalty1, matchingScore,
                                                                                                                               mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                               scoreThreshold, w, xDrop);

                            std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                            pairedSimilarFragments0 = pairedSimilarFragments;

                            outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, lastStrand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);
                        }
                        */
                    }
                    mafStrand="+";
                }else if ( lastStrand == NEGATIVE && alignmentMatch.getStrand() == NEGATIVE
                             && lastAlignmentMatch.getRefStartPos() > alignmentMatch.getRefEndPos()
                             && lastAlignmentMatch.getQueryEndPos() < alignmentMatch.getQueryStartPos()
                         ){
                    if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() != startQuery) {
                        endQuery = alignmentMatch.getQueryEndPos() + 1;
                        std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery, alignmentMatch.getStrand());
                        for (int repeatI = 0; repeatI < querySeq.length(); ++repeatI) {
                            refAlign << "-";
                        }
                        queryAlign << querySeq;
                        alignmentScore -= openGapPenalty1 + extendGapPenalty1*querySeq.length();
                    } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryEndPos() == startQuery) {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                        refAlign << refSeq;
                        for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                            queryAlign << "-";
                        }
                        alignmentScore -= openGapPenalty1 + extendGapPenalty1*refSeq.length();
                    } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() == startQuery) {

                    } else {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        endQuery = alignmentMatch.getQueryEndPos() + 1;

                        std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                        std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                        std::string _alignment_q;
                        std::string _alignment_d;

                        int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth, wfaSize,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);
                        std::string temp;

                        temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq) == 0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq) == 0);

                        alignmentScore += thiScore;
                        refAlign << _alignment_d;
                        queryAlign << _alignment_q;

                        if (outPutFraged ){
                            ofragfile << "a\tscore=" << thiScore << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << endQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofragfile << std::endl;
                        }
                        /*
                        if ( outPutLocalalignment ){

                            int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                            int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                            int32_t length1 = refSeq.length();


                            int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                            int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                            int32_t length2 = querySeq.length();

                            std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                               seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                               mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                               extendGapPenalty1, matchingScore,
                                                                                                                               mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                               scoreThreshold, w, xDrop);

                            std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                            pairedSimilarFragments0 = pairedSimilarFragments;

                            outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, lastStrand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);
                        }
                         */
                    }
                    if ( mafQueryStart > alignmentMatch.getQueryStartPos() ){
                        mafQueryStart = alignmentMatch.getQueryStartPos();
                    }
                }else{
//                    std::cout << "line 710" << std::endl;
                    std::string temp1 = refAlign.str();
                    std::string temp2 = queryAlign.str();
                    if( outPutMaf && temp1.size()>0 ) {
                        temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
                        temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());
                        omaffile << "a\tscore=" << alignmentScore << std::endl;
                        omaffile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right
                                 << std::setw(9) << mafRefStart << "\t" << std::setw(9)
                                 << temp1.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t"
                                 << refAlign.str() << std::endl;

                        omaffile << "s\t" << std::left << std::setw(chrWidth) << queryFileName + "." + queryChr << "\t" << std::right
                                 << std::setw(9) << mafQueryStart << "\t" << std::setw(9)
                                 << temp2.size() << "\t" << mafStrand << "\t"<< targetSequences[queryChr].getSequence().size() << "\t"
                                 << queryAlign.str() << std::endl;
                        omaffile << std::endl;
                    }
                    alignmentScore = 0;
                    refAlign.str(std::string());
                    queryAlign.str(std::string());

                    mafRefStart = alignmentMatch.getRefStartPos();
                    mafQueryStart = alignmentMatch.getQueryStartPos();
                }
                {
//                    std::cout << "line 751" << std::endl;
                    if( POSITIVE == alignmentMatch.getStrand() ){
                        mafStrand = "+";
                    }else{
                        mafStrand = "-";
                    }
                    startRef = alignmentMatch.getRefStartPos();
                    startQuery = alignmentMatch.getQueryStartPos();
                    endRef = alignmentMatch.getRefEndPos();
                    endQuery = alignmentMatch.getQueryEndPos();
                    std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                    std::string _alignment_q;
                    std::string _alignment_d;

                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth,  wfaSize2, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);

                    std::string temp;

                    temp = _alignment_d;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(refSeq) == 0);
                    temp = _alignment_q;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(querySeq) == 0);

                    alignmentScore += thiScore;
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    if (outPutFraged){
                        ofragfile << "a\tscore=" << thiScore << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                        ofragfile << std::endl;
                    }
                    /*
                    if ( outPutLocalalignment ){

                        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                        int32_t length1 = refSeq.length();


                        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                        int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                        int32_t length2 = querySeq.length();

                        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                           seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                           mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                           extendGapPenalty1, matchingScore,
                                                                                                                           mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                           scoreThreshold, w, xDrop);

                        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                        pairedSimilarFragments0 = pairedSimilarFragments;
                        STRAND thisStrand = alignmentMatch.getStrand();
                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, thisStrand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);
                    }
                     */
                }
                startRef = alignmentMatch.getRefEndPos() + 1;
                startQuery = alignmentMatch.getQueryEndPos() + 1;
                if(alignmentMatch.getStrand() == NEGATIVE){
                    startQuery = alignmentMatch.getQueryStartPos() - 1;
                    std::cout << "line 775" << std::endl;
                }
                lastStrand = alignmentMatch.getStrand();
                lastAlignmentMatch = alignmentMatch;
            }
            if( !hasInversion ){

                std::cout << refChr << " last align" << std::endl;
                endRef = refSequences[refChr].getSequence().length();
                endQuery = targetSequences[queryChr].getSequence().length();
                if (startRef > endRef && startQuery <= endQuery) {
                    std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);
                    for (int repeatI = 0; repeatI < querySeq.length(); ++repeatI) {
                        refAlign << "-";
                    }
                    queryAlign << querySeq;
                    alignmentScore += -6 + -2*querySeq.length();
                } else if (startRef <= endRef && startQuery > endQuery) {
                    std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                    refAlign << refSeq;
                    for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                        queryAlign << "-";
                    }
                    alignmentScore += -6 + -2*refSeq.length();

                } else if (startRef > endRef && startQuery > endQuery) {

                } else {
                    std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);

                    std::string _alignment_q;
                    std::string _alignment_d;

                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, &penalties, mm_allocator, widownWidth, wfaSize,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2, min_wavefront_length, max_distance_threshold, m);
                    std::string temp;

                    temp = _alignment_d;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(refSeq) == 0);
                    temp = _alignment_q;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(querySeq) == 0);

                    alignmentScore += thiScore;
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    if (outPutFraged ){
                        ofragfile << "a\tscore=" << thiScore << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                        ofragfile << std::endl;
                    }/*
                    if ( outPutLocalalignment ){

                        int8_t * seq1 = sequenceCharToUInt8.seq_to_int8(refSeq);
                        int8_t * seq1_rev_com = sequenceCharToUInt8.rev_comp(seq1, refSeq.length());
                        int32_t length1 = refSeq.length();


                        int8_t * seq2 = sequenceCharToUInt8.seq_to_int8(querySeq);
                        int8_t * seq2_rev_com = sequenceCharToUInt8.rev_comp(seq2, querySeq.length());
                        int32_t length2 = querySeq.length();

                        std::vector<PairedSimilarFragment> pairedSimilarFragments0 = findSimilarFragmentsForPairedSequence(seq1, seq1_rev_com,
                                                                                                                           seq2, seq2_rev_com, length1, length2, seed_window_size,
                                                                                                                           mini_cns_score, matrix_boundary_distance, openGapPenalty1,
                                                                                                                           extendGapPenalty1, matchingScore,
                                                                                                                           mismatchingPenalty, m, step_size, refSeq, querySeq,
                                                                                                                           scoreThreshold, w, xDrop);

                        std::vector<PairedSimilarFragment> pairedSimilarFragments = syntenic(pairedSimilarFragments0);
                        pairedSimilarFragments0 = pairedSimilarFragments;
                        STRAND thisStrand = POSITIVE;
                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, thisStrand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, oLocalalignment, refSequences, targetSequences);
                    }*/
                }

                std::string temp = refAlign.str();
                temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                assert(temp.compare(refSequences[refChr].getSequence()) == 0);

                temp = queryAlign.str();
                temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
//                temp.erase(std::remove(temp.begin(), temp.end(), '\0'), temp.end());
                assert(temp.compare(targetSequences[queryChr].getSequence()) == 0);
            }

            std::cout << refChr << " align done" << std::endl;


            if( outPutMaf ) {
                std::string temp1 = refAlign.str();
                std::string temp2 = queryAlign.str();
                temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
                temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());
                omaffile << "a\tscore=" << alignmentScore << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right
                      << std::setw(9) << mafRefStart << "\t" << std::setw(9)
                      << temp1.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t"
                      << refAlign.str() << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth) << queryFileName + "." + queryChr << "\t" << std::right
                     << std::setw(9) << mafQueryStart << "\t" << std::setw(9)
                     << temp2.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t"
                     << queryAlign.str() << std::endl;
                omaffile << std::endl;
            }

            if (outPutVcf) {
                std::string queryAlignSeq = queryAlign.str();
                std::string refAlignSeq = refAlign.str();
                alignmentToVcf(queryAlignSeq, refAlignSeq, ovcffile, it0->first, refSequences);
            }
        }
    }

    if( outPutMaf ){
        omaffile.close();
    }
    if( outPutVcf ){
        ovcffile.close();
    }
    if( outPutFraged ){
        ofragfile.close();
    }
//    if(outPutLocalalignment){
//        oLocalalignment.close();
//    }
    mm_allocator_delete(mm_allocator);
}
