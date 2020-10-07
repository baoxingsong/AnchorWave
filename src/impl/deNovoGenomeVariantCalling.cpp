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

            ofile2 << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << refPos-1 << "\t" << std::setw(9) << temp.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
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

            ofile2 << "a\tscore=" << pairedSimilarFragments0[i].getScore() << std::endl;
            ofile2 << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right << std::setw(9) << refPos - 1 << "\t" << std::setw(9) << temp.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            ofile2 << "s\t" << std::left << std::setw(chrWidth) << queryFileName + "." + queryChr << "\t" << std::right << std::setw(9) << quePos - 1 << "\t" << std::setw(9) << temp.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" << queryAlign.str() << std::endl;
            ofile2 << std::endl;
        }
    }
}




void genomeAlignment( std::vector<std::vector<OrthologPair2>> & alignmentMatchsMap,
                      const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                      const size_t & widownWidth, const std::string & outPutFilePath, bool & outPutAlignmentForEachInterval,
                      const bool & localAlignment, const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
                      const int32_t & extendGapPenalty1, int32_t & seed_window_size, const int32_t & mini_cns_score, const int32_t & step_size,
                      const int32_t & matrix_boundary_distance, const  int32_t & scoreThreshold, const  int32_t & w, const  int32_t & xDrop) {


    Scorei m(matchingScore, mismatchingPenalty);

    std::ofstream ofile;
    ofile.open(outPutFilePath+".maf");
    ofile << "##maf version=1 scoring=PROALI" << std::endl;


    std::map <std::string, Fasta> refSequences;
    readFastaFile(refFastaFilePath, refSequences);

    std::map <std::string, Fasta> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);
    int32_t sizei = alignmentMatchsMap.size();

    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    std::string refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    std::string queryFileName = elems.back();
    int chrWidth = 4;
    for ( std::map <std::string, Fasta>::iterator itchr  =refSequences.begin(); itchr != refSequences.end(); ++itchr ){
        if ( (refFileName + "." + itchr->first).size() > chrWidth){
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }

    for ( std::map <std::string, Fasta>::iterator itchr  =targetSequences.begin(); itchr != targetSequences.end(); ++itchr ){
        if ( (queryFileName + "." + itchr->first).size() > chrWidth ){
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    SequenceCharToUInt8 sequenceCharToUInt8;

    for(int32_t i=sizei-1; i>=0; --i){
//        int32_t sizej = alignmentMatchsMap[i].size();
        size_t startRef = alignmentMatchsMap[i][0].getRefStartPos();
        size_t startQuery;
        STRAND strand = alignmentMatchsMap[i][0].getStrand();
        size_t endRef;
        size_t endQuery;
        std::stringstream refAlign;
        std::stringstream queryAlign;
        std::string refChr = alignmentMatchsMap[i][0].getRefChr();
        std::string queryChr = alignmentMatchsMap[i][0].getQueryChr();
//        std::cout << refChr << " align start" << std::endl;
        if(  POSITIVE == strand ){
            int64_t alignmentScore = 0;
            startQuery = alignmentMatchsMap[i][0].getQueryStartPos();
            for( OrthologPair2 orthologPair : alignmentMatchsMap[i] ){
//                std::cout  << orthologPair.getRefChr() << "\t"
//                           << orthologPair.getRefStartPos() << "\t"
//                           << orthologPair.getRefEndPos() << "\t"
//                           << orthologPair.getQueryChr() << "\t"
//                           << orthologPair.getQueryStartPos() << "\t"
//                           << orthologPair.getQueryEndPos() << "\t"
//                           << orthologPair.getStrand() << "\t"
//                           << orthologPair.getReferenceGeneName() << "\t"
//                           << orthologPair.getQueryGeneName() << std::endl;
                if (  orthologPair.getRefStartPos()==startRef && orthologPair.getQueryStartPos()!=startQuery ) {
                    endQuery = orthologPair.getQueryStartPos()-1;
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery);
                    for ( int repeatI = 0; repeatI<querySeq.length(); ++repeatI ){
                        refAlign<<"-";
                    }
                    queryAlign<<querySeq;
                    alignmentScore += -6 + -2*querySeq.length();
                }else if ( orthologPair.getRefStartPos()!=startRef && orthologPair.getQueryStartPos()==startQuery ){
                    endRef = orthologPair.getRefStartPos()-1;
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    refAlign<<refSeq;
                    for ( int repeatI = 0; repeatI<refSeq.length(); ++repeatI ) {
                        queryAlign << "-";
                    }
                    alignmentScore += -6 + -2*refSeq.length();
                }else if ( orthologPair.getRefStartPos()==startRef && orthologPair.getQueryStartPos()==startQuery ){

                }else{
                    endRef = orthologPair.getRefStartPos()-1;
                    endQuery = orthologPair.getQueryStartPos()-1;

                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery);

                    if ( localAlignment ){

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


                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, ofile, refSequences, targetSequences);

                    }else{

                        std::string _alignment_q;
                        std::string _alignment_d;
                        int64_t thiScore = alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1  );
                        alignmentScore += thiScore;
                        refAlign<<_alignment_d;
                        queryAlign<<_alignment_q;

                        std::string temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq)==0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq)==0);

                        if (outPutAlignmentForEachInterval && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                            ofile << "a\tscore=" << thiScore << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofile << std::endl;
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

                    if ( localAlignment ){

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

                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, ofile, refSequences, targetSequences);

                    }else{
                        std::string _alignment_q;
                        std::string _alignment_d;
                        int64_t thiScore =alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1  );
                        alignmentScore += thiScore;
                        refAlign<<_alignment_d;
                        queryAlign<<_alignment_q;

                        std::string temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq)==0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq)==0);

                        if (outPutAlignmentForEachInterval && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                            ofile << "a\tscore=" << thiScore << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofile << std::endl;

                        }
                    }
                }
                startRef = orthologPair.getRefEndPos()+1;
                startQuery= orthologPair.getQueryEndPos()+1;
            }
//            std::cout << refChr << " last align" << std::endl;

            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());

            std::string refGenomerSequence = getSubsequence( refSequences, refChr, alignmentMatchsMap[i][0].getRefStartPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getRefEndPos());
            //assert(temp.compare(refGenomerSequence)==0);

            std::string queryGenomerSequence = getSubsequence( targetSequences, queryChr, alignmentMatchsMap[i][0].getQueryStartPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getQueryEndPos());

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            //assert(temp.compare(queryGenomerSequence)==0);
            if( !outPutAlignmentForEachInterval && !localAlignment  ){
                ofile << "a\tscore=" << alignmentScore << std::endl;

                ofile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][0].getRefStartPos()-1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;
                ofile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][0].getQueryStartPos()-1 << "\t" << std::setw(9) << queryGenomerSequence.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  queryAlign.str() << std::endl;
                ofile << std::endl;
            }
        } else {
            int64_t alignmentScore = 0;
            endQuery = alignmentMatchsMap[i][0].getQueryEndPos();
            for( OrthologPair2 orthologPair : alignmentMatchsMap[i] ){
//
//                std::cout  << orthologPair.getRefChr() << "\t"
//                           << orthologPair.getRefStartPos() << "\t"
//                           << orthologPair.getRefEndPos() << "\t"
//                           << orthologPair.getQueryChr() << "\t"
//                           << orthologPair.getQueryStartPos() << "\t"
//                           << orthologPair.getQueryEndPos() << "\t"
//                           << orthologPair.getStrand() << "\t"
//                           << orthologPair.getReferenceGeneName() << "\t"
//                           << orthologPair.getQueryGeneName() << std::endl;

                if (  orthologPair.getRefStartPos()==startRef && orthologPair.getQueryEndPos()!=endQuery ) {
                    startQuery = orthologPair.getQueryEndPos()+1;
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery, strand);
                    for ( int repeatI = 0; repeatI<querySeq.length(); ++repeatI ){
                        refAlign<<"-";
                    }
                    queryAlign<<querySeq;
                    alignmentScore += -6 + -2*querySeq.length();
                }else if ( orthologPair.getRefStartPos()!=startRef && orthologPair.getQueryEndPos()==endQuery ){
                    endRef = orthologPair.getRefStartPos()-1;
                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    refAlign<<refSeq;
                    for ( int repeatI = 0; repeatI<refSeq.length(); ++repeatI ) {
                        queryAlign << "-";
                    }
                    alignmentScore += -6 + -2*refSeq.length();
                }else if ( orthologPair.getRefStartPos()==startRef && orthologPair.getQueryEndPos()==endQuery ){

                }else{
                    endRef = orthologPair.getRefStartPos()-1;
                    startQuery = orthologPair.getQueryEndPos()+1;

                    std::string refSeq = getSubsequence( refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence( targetSequences, queryChr, startQuery, endQuery, strand);




                    if ( localAlignment ){

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

                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, ofile, refSequences, targetSequences);

                    }else{
                        //                ofile << refSeq.length() << "\t" << querySeq.length() << std::endl;

                        std::string _alignment_q;
                        std::string _alignment_d;
                        int64_t thiScore =alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1  );
                        alignmentScore += thiScore;
                        refAlign<<_alignment_d;
                        queryAlign<<_alignment_q;

                        std::string temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq)==0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq)==0);
                        if (outPutAlignmentForEachInterval && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                            ofile << "a\tscore=" << thiScore << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofile << std::endl;
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

                    if ( localAlignment ){

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

                        outputLocalAlignment(chrWidth, refFileName, queryFileName, pairedSimilarFragments0, strand, startRef, startQuery, refChr, queryChr,refSeq, querySeq, ofile, refSequences, targetSequences);
                    }else{
                        std::string _alignment_q;
                        std::string _alignment_d;
                        int64_t thiScore = alignSlidingWindow( querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1 );
                        alignmentScore += thiScore;
                        refAlign<<_alignment_d;
                        queryAlign<<_alignment_q;

                        std::string temp = _alignment_d;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(refSeq)==0);
                        temp = _alignment_q;
                        temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                        assert(temp.compare(querySeq)==0);

                        if (outPutAlignmentForEachInterval && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                            ofile << "a\tscore=" << thiScore << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                            ofile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                            ofile << std::endl;
                        }
                    }
                }
                startRef = orthologPair.getRefEndPos()+1;
                endQuery= orthologPair.getQueryStartPos()-1;
            }
//            std::cout << refChr << " last align" << std::endl;

            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());

            std::string refGenomerSequence = getSubsequence( refSequences, refChr, alignmentMatchsMap[i][0].getRefStartPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getRefEndPos());
            //assert(temp.compare(refGenomerSequence)==0);

            std::string queryGenomerSequence = getSubsequence( targetSequences, queryChr, alignmentMatchsMap[i][0].getQueryEndPos(), alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getQueryStartPos(), strand);
            if ( !outPutAlignmentForEachInterval && !localAlignment ){
                ofile << "a\tscore=" << alignmentScore << std::endl;
                ofile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][0].getRefStartPos()-1 << "\t" << std::setw(9) << refGenomerSequence.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << refAlign.str() << std::endl;
                ofile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << alignmentMatchsMap[i][alignmentMatchsMap[i].size()-1].getQueryStartPos()-1 << "\t" << std::setw(9) << queryGenomerSequence.size() << "\t-\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  queryAlign.str() << std::endl;
                ofile << std::endl;
            }
            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
            //assert(temp.compare(queryGenomerSequence)==0);
        }
    }
    ofile.close();

}

void deNovoGenomeVariantCalling( std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                 const std::string & refFastaFilePath, const std::string & targetFastaFilePath,
                                 const size_t & widownWidth, const std::string & outPutMafFile, const std::string & outPutVcfFile,
                                 const std::string & outPutFragedFile, const int32_t & matchingScore, const  int32_t & mismatchingPenalty, const  int32_t & openGapPenalty1,
                                 const int32_t & extendGapPenalty1) {
    bool outPutMaf = false;
    bool outPutVcf = false;
    bool outPutFraged = false;

    if ( outPutMafFile.size() > 0 ){
        outPutMaf = true;
    }
    if ( outPutVcfFile.size() > 0 ){
        outPutVcf = true;
    }
    if ( outPutFragedFile.size() > 0 ){
        outPutFraged = true;
    }

    std::map <std::string, Fasta> refSequences;
    readFastaFile(refFastaFilePath, refSequences);

    std::map <std::string, Fasta> targetSequences;
    readFastaFile(targetFastaFilePath, targetSequences);

    int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    if ( outPutMaf || outPutFraged ){
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
    }

    std::ofstream omaffile;
    std::ofstream ovcffile;
    std::ofstream ofragfile;
    if( outPutMaf ){
        omaffile.open(outPutMafFile);
    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon)) + std::to_string((ltm->tm_mday));


    if( outPutVcf ){
        ovcffile.open(outPutVcfFile);
        ovcffile << "##fileformat=VCFv4.3" << std::endl;
        ovcffile << "##fileDate=" << filedate << std::endl;
        ovcffile << "##source=proali" << std::endl;
        ovcffile <<"##reference=" << refFileName << std::endl;
        ovcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER" << std::endl;
    }

    if( outPutFraged ){
        ofragfile.open(outPutFragedFile);
    }

    for ( std::map<std::string, std::vector<AlignmentMatch>>::iterator it0=alignmentMatchsMap.begin(); it0!=alignmentMatchsMap.end(); ++it0 ) {
        if( it0->second.size() > 0 ){
            size_t startRef = 1;
            size_t startQuery = 1;
            size_t endRef;
            size_t endQuery;
            std::stringstream refAlign;
            std::stringstream queryAlign;

            std::string refChr = it0->second[0].getDatabaseChr();
            std::string queryChr = it0->second[0].getQueryChr();
            std::cout << refChr << " align start" << std::endl;

            int64_t alignmentScore = 0;

            for (AlignmentMatch alignmentMatch : it0->second) {

                refChr = alignmentMatch.getDatabaseChr();
                queryChr = alignmentMatch.getQueryChr();
                if (alignmentMatch.getDatabaseStart() == startRef && alignmentMatch.getQueryStart() != startQuery) {
                    endQuery = alignmentMatch.getQueryStart() - 1;
                    std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);
                    for (int repeatI = 0; repeatI < querySeq.length(); ++repeatI) {
                        refAlign << "-";
                    }
                    queryAlign << querySeq;
                    alignmentScore += -6 + -2*querySeq.length();
                } else if (alignmentMatch.getDatabaseStart() != startRef && alignmentMatch.getQueryStart() == startQuery) {
                    endRef = alignmentMatch.getDatabaseStart() - 1;
                    std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                    refAlign << refSeq;
                    for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                        queryAlign << "-";
                    }
                    alignmentScore += -6 + -2*refSeq.length();

                } else if (alignmentMatch.getDatabaseStart() == startRef && alignmentMatch.getQueryStart() == startQuery) {

                } else {
                    endRef = alignmentMatch.getDatabaseStart() - 1;
                    endQuery = alignmentMatch.getQueryStart() - 1;

                    std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);

                    std::string _alignment_q;
                    std::string _alignment_d;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1 );
                    alignmentScore += thiScore;

                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    std::string temp = _alignment_d;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(refSeq) == 0);
                    temp = _alignment_q;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(querySeq) == 0);

                    if (outPutFraged && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                        ofragfile << "a\tscore=" << thiScore << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                        ofragfile << std::endl;
                    }
                }
                {
                    startRef = alignmentMatch.getDatabaseStart();
                    startQuery = alignmentMatch.getQueryStart();
                    endRef = alignmentMatch.getDatabaseEnd();
                    endQuery = alignmentMatch.getQueryEnd();
                    std::string refSeq = getSubsequence(refSequences, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence(targetSequences, queryChr, startQuery, endQuery);

    //                ofile << refSeq.length() << "\t" << querySeq.length() << std::endl;

                    std::string _alignment_q;
                    std::string _alignment_d;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1 );
                    alignmentScore += thiScore;
                    refAlign << _alignment_d;
                    queryAlign << _alignment_q;

                    std::string temp = _alignment_d;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(refSeq) == 0);
                    temp = _alignment_q;
                    temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                    assert(temp.compare(querySeq) == 0);

                    if (outPutFraged && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                        ofragfile << "a\tscore=" << thiScore << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                        ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                        ofragfile << std::endl;
                    }

                }
                startRef = alignmentMatch.getDatabaseEnd() + 1;
                startQuery = alignmentMatch.getQueryEnd() + 1;
    //            ++alignmentMatchIndex;
            }
            std::cout << refChr << " last align" << std::endl;
            endRef = refSequences[refChr].getSequence().length();
            endQuery = targetSequences[queryChr].getSequence().length();
    //        std::cout << refChr << " last line 789" << std::endl;
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
                //std::cout << refChr << " last line 811" << std::endl << refSeq << std::endl << querySeq << std::endl ;
                int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, widownWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1 );
                alignmentScore += thiScore;
                refAlign << _alignment_d;
                queryAlign << _alignment_q;

                std::string temp = _alignment_d;
                temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                assert(temp.compare(refSeq) == 0);
                temp = _alignment_q;
                temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());
                assert(temp.compare(querySeq) == 0);

                if (outPutFraged && ((refSeq.size() <=widownWidth || querySeq.size() <=widownWidth)) && refSeq.size() <=(2*widownWidth) && querySeq.size()<=(2*widownWidth) ){
                    ofragfile << "a\tscore=" << thiScore << std::endl;
                    ofragfile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right <<  std::setw(9) << startRef-1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t" << _alignment_d << std::endl;
                    ofragfile << "s\t" << std::left << std::setw(chrWidth)  << queryFileName + "." + queryChr << "\t" << std::right <<  std::setw(9) << startQuery-1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << targetSequences[queryChr].getSequence().size() << "\t" <<  _alignment_q << std::endl;
                    ofragfile << std::endl;
                }
            }
    //        std::cout << refChr << " last line 822" << std::endl;


            std::string temp = refAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());


            assert(temp.compare(refSequences[refChr].getSequence()) == 0);

            temp = queryAlign.str();
            temp.erase(std::remove(temp.begin(), temp.end(), '-'), temp.end());

            assert(temp.compare(targetSequences[queryChr].getSequence()) == 0);

            std::cout << refChr << " align done" << std::endl;

            if( outPutMaf ) {

                omaffile << "a\tscore=" << alignmentScore << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth) << refFileName + "." + refChr << "\t" << std::right
                      << std::setw(9) << 0 << "\t" << std::setw(9)
                      << refSequences[refChr].getSequence().size() << "\t+\t" << refSequences[refChr].getSequence().size() << "\t"
                      << refAlign.str() << std::endl;
                omaffile << "s\t" << std::left << std::setw(chrWidth) << queryFileName + "." + queryChr << "\t"
                      << std::right << std::setw(9)
                      << 0 << "\t"
                      << std::setw(9) << 0 << "\t+\t"
                      << targetSequences[queryChr].getSequence().size() << "\t" << queryAlign.str() << std::endl;
                omaffile << std::endl;
            }


            if (outPutVcf) {
                std::cout << std::endl;
                std::string queryAlignSeq = queryAlign.str();
                std::string refAlignSeq = refAlign.str();
                FirstLastList sdiRecords;
                int refLetterNumber = 0;
                for (int ai = 0; ai < refAlignSeq.length(); ai++) {
                    if (refAlignSeq[ai] != '-') {
                        ++refLetterNumber;
                    }
                    if (refAlignSeq[ai] != queryAlignSeq[ai]) {
                        if (queryAlignSeq[ai] == '-') {
                            std::string ori(1, refAlignSeq[ai]);
                            std::string result = "-";
                            Variant mapSingleRecord = Variant(it0->first, refLetterNumber, ori, result);
                            Data *data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
                        } else if (refAlignSeq[ai] == '-') {
                            if (sdiRecords.getLast() == NULL) {
                                int position = refLetterNumber + 1;
                                std::string ori = "-";
                                std::string result(1, queryAlignSeq[ai]);
                                Variant mapSingleRecord = Variant(it0->first, position, ori, result);
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
                                    std::string result =
                                            sdiRecords.getLast()->getMapSingleRecord().getAlternative() +
                                            queryAlignSeq[ai];
                                    Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                                    Data *data = new Data(mapSingleRecord);
                                    sdiRecords.deleteLast();
                                    sdiRecords.insertLast(data);
                                } else {
                                    int position = refLetterNumber + 1;
                                    std::string ori = "-";
                                    std::string result(1, queryAlignSeq[ai]);
                                    Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                                    Data *data = new Data(mapSingleRecord);
                                    sdiRecords.insertLast(data);
                                }
                            }
                        } else {
                            int position = refLetterNumber;
                            std::string ori(1, refAlignSeq[ai]);
                            std::string result(1, queryAlignSeq[ai]);
                            Variant mapSingleRecord = Variant(it0->first, position, ori, result);
                            Data *data = new Data(mapSingleRecord);
                            sdiRecords.insertLast(data);
                        }
                    }
                }

                for (int runingCound = 0; runingCound < 2; ++runingCound) {
                    std::cout << "round " << runingCound << std::endl;
                    if ((sdiRecords.getFirst() != NULL)
                        && (sdiRecords.getFirst()->getNext() != NULL)) {
                        Data *prevOne = sdiRecords.getFirst();
                        Data *currOne = (sdiRecords.getFirst()->getNext());
                        //                    std::cout << "2107" << std::endl;
                        while ((currOne != NULL) && (NULL != currOne->getNext())) {
                            //                        std::cout << "2109" << std::endl;
                            if (sdiRecords.getFirst() == currOne) {
                                prevOne = currOne;
                                currOne = prevOne->getNext();
                            }
                            //std::cout << "793" << std::endl;
                            if (currOne->getMapSingleRecord().getChanginglength() < 0 &&
                                prevOne->getMapSingleRecord().getChanginglength() < 0 &&
                                (prevOne->getMapSingleRecord().getPosition() +
                                 abs(prevOne->getMapSingleRecord().getChanginglength())) ==
                                currOne->getMapSingleRecord().getPosition()
                                && prevOne->getMapSingleRecord().getAlternative().compare("-") == 0 &&
                                currOne->getMapSingleRecord().getAlternative().compare("-") ==
                                0) { // merge to deletions
                                //std::cout << "797 delete prev" << std::endl;
                                int position = prevOne->getMapSingleRecord().getPosition();
                                std::string ori = prevOne->getMapSingleRecord().getReference() +
                                                  currOne->getMapSingleRecord().getReference();
                                std::string result = "-";
                                Variant mapSingleRecord2(it0->first, position, ori, result);
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
                            } else if (currOne->getMapSingleRecord().getChanginglength() == 0 && 0 ==
                                                                                                 currOne->getMapSingleRecord().getReference().compare(
                                                                                                         currOne->getMapSingleRecord().getAlternative())) { // nonsense records
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

                //        std::cout << " begin to sort" << std::endl;
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
                            std::string ori =
                                    sdiRecordsThisOne[j - 1].getReference() + sdiRecordsThisOne[j].getReference();
                            std::string result = "-";
                            Variant mapSingleRecord2(it0->first, position, ori, result);
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

                for (std::vector<Variant>::iterator itVariant = sdiRecordsThisOne.begin();
                     itVariant != sdiRecordsThisOne.end(); ++itVariant) {
                    ovcffile << itVariant->getChromosome() << "\t" << itVariant->getPosition() << "\t" << itVariant->getChromosome() << "_" << itVariant->getPosition() <<"\t"+
                             itVariant->getReference() << "\t" << itVariant->getAlternative() << "\t50\tPASS" << std::endl;
                }
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
}
