//
// Created by song on 8/4/18.
//

#include "TransferGffWithNucmerResult.h"

void splitCIGAR( std::string & cigarString, std::vector<std::string> & cigarElems) {
    std::regex reg("([0-9]+[MIDNSHPX=])");
    std::smatch match;
    while( regex_search(cigarString, match, reg) ){
        cigarElems.push_back( match[1] );
        cigarString = match.suffix().str();
    }
}

void readSam(std::vector<AlignmentMatch> & alignmentMatchsMapT, std::ifstream & infile, std::map<std::string, Transcript> & transcriptHashMap){

    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;

    std::string databaseChr;
    size_t databaseStart;
    size_t databaseEnd;
    std::string queryChr;
    size_t queryStart;
    size_t queryEnd;

    while (std::getline(infile, line)){ // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if( line[0] != '@' ){ //ignore the header
            elems.clear();
            split(line, delim, elems);
            queryStart=stoi(elems[3]);
//            if(databaseStart>0 && transcriptHashMap.find(elems[0])!=transcriptHashMap.end() ){ // ignore those none mapping records
            queryChr=elems[2];
            if( queryChr.compare("*") != 0 && transcriptHashMap.find(elems[0])!=transcriptHashMap.end() ){ // ignore those none mapping records
                //std::cout << "begain to analysis " << line << std::endl;
                databaseChr=transcriptHashMap[elems[0]].getChromeSomeName();

                databaseStart = transcriptHashMap[elems[0]].getPStart();
                databaseEnd = transcriptHashMap[elems[0]].getPEnd();
//                databaseEnd=databaseStart;
                queryEnd=queryStart;
                std::vector<std::string> cigarElems;
                splitCIGAR(elems[5], cigarElems);
//                std::cout << "cigar split done" << std::endl;
                int headClipping = 0;
                int tailClipping = 0;
                for(int i=0; i<cigarElems.size(); ++i) {
                    std::string cVal = cigarElems[i];
                    char cLetter = cVal[cVal.length() - 1];
                    int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                    if( i == cigarElems.size()-1 && (cLetter == 'H' || cLetter == 'S' ) ){ // ignore the last soft/hard clipping
                        tailClipping += cLen;
                        continue;
                    }

                    switch (cLetter) {
                        case 'H':
                            headClipping += cLen;
//                            databaseEnd += cLen;
                            break;
                        case 'S':
                            headClipping += cLen;
//                            databaseEnd += cLen;
                            break;
                        case 'M':
                            queryEnd += cLen;
//                            databaseEnd += cLen;
                            break;
                        case '=':
                            queryEnd += cLen;
//                            databaseEnd += cLen;
                            break;
                        case 'X':
                            queryEnd += cLen;
//                            databaseEnd += cLen;
                            break;
                        case 'I':
//                            databaseEnd += cLen;
                            break;
                        case 'D':
                            queryEnd += cLen;
                            break;
                        case 'N':
                            queryEnd += cLen;
                            break;
                        case 'P':
                            break;
                        default:
                            std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                            break;
                    }
                }

                --queryEnd;
                int samFlag = stoi(elems[1]);

                std::map<int32_t, int32_t> positionsMap;
                size_t cdsSequenceLength = 0;
                if (transcriptHashMap[elems[0]].getStrand()==POSITIVE ){
                    int32_t cdsPosition = 0;
                    int32_t chromosomePosition = transcriptHashMap[elems[0]].getPStart()-1;
                    for ( int32_t cdsIndex = 0; cdsIndex < transcriptHashMap[elems[0]].getCdsVector().size(); ++cdsIndex ){
                        //std::cout << transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart() << "\t" << transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd() << std::endl;
                        for ( int32_t i=transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart(); i<=transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd(); ++i ){
                            cdsPosition++;
                            chromosomePosition++;
                            cdsSequenceLength++;
                            positionsMap[cdsPosition]=chromosomePosition;
                        }
                        if ( cdsIndex < transcriptHashMap[elems[0]].getCdsVector().size()-1 ){ // for intron
                            for ( int32_t i=transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd()+1; i<transcriptHashMap[elems[0]].getCdsVector()[cdsIndex+1].getStart(); ++i ){
                                chromosomePosition++;
                            }
                        }
                    }
                    assert( chromosomePosition==transcriptHashMap[elems[0]].getPEnd());
                }else{
                    int32_t cdsPosition = 0;
                    int32_t chromosomePosition = transcriptHashMap[elems[0]].getPEnd() + 1;
                    for ( int32_t cdsIndex = transcriptHashMap[elems[0]].getCdsVector().size()-1; cdsIndex >=0 ; --cdsIndex ){
                        for ( int32_t i=transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getEnd(); i>=transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart(); --i ){
                            cdsPosition++;
                            chromosomePosition--;
                            cdsSequenceLength++;
                            positionsMap[cdsPosition]=chromosomePosition;
                        }
                        if ( cdsIndex >0 ){ // for intron
                            for ( int32_t i=transcriptHashMap[elems[0]].getCdsVector()[cdsIndex].getStart()-1; i>transcriptHashMap[elems[0]].getCdsVector()[cdsIndex-1].getEnd(); --i ){
                                chromosomePosition--;
                            }
                        }
                    }
                    assert(transcriptHashMap[elems[0]].getPStart() == chromosomePosition);
                }

                if ( 0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==POSITIVE ){
                    databaseStart = positionsMap[headClipping+1];
                    databaseEnd   = positionsMap[cdsSequenceLength-tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if ( 0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==NEGATIVE ){
                    databaseEnd    = positionsMap[headClipping+1];
                    databaseStart  = positionsMap[cdsSequenceLength-tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if ( 0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==POSITIVE ) {
                    databaseEnd = positionsMap[cdsSequenceLength - headClipping];
                    databaseStart   = positionsMap[1 + tailClipping];
                    assert(databaseStart < databaseEnd);
                } else if ( 0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==NEGATIVE ) {
                    databaseStart = positionsMap[cdsSequenceLength - headClipping];
                    databaseEnd   = positionsMap[1 + tailClipping];
                    assert(databaseStart < databaseEnd);
                }

                double thisScore = 1.0 - (tailClipping+headClipping)/(double)cdsSequenceLength;
                if( (0 == samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==POSITIVE)
                    || (0 != samFlag % 32 && transcriptHashMap[elems[0]].getStrand()==NEGATIVE) ){

                    AlignmentMatch orthologPair(databaseChr, queryChr,
                                                databaseStart, databaseEnd, queryStart, queryEnd, thisScore, POSITIVE, elems[0], elems[0] );
                    alignmentMatchsMapT.push_back(orthologPair);
                }else{
                    AlignmentMatch orthologPair(databaseChr, queryChr,
                                                databaseStart, databaseEnd, queryStart, queryEnd, thisScore, NEGATIVE, elems[0], elems[0] );
                    alignmentMatchsMapT.push_back(orthologPair);
                }
            }
        }
    }

}


void setupAnchorsWithSpliceAlignmentResult( const std::string & gffFilePath, const std::string & samFile, std::map<std::string, std::vector<AlignmentMatch>> & alignmentMatchsMap, double & inversion_PENALTY, double & MIN_ALIGNMENT_SCORE, bool & considerInversion){
    std::ifstream infile(samFile);
    if( ! infile.good()){
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit (1);
    }

    std::map<std::string, std::vector<std::string> > geneNameMap; // key is chromosome and value is gene names
    std::map<std::string, Gene> geneHashMap;  // key is gene name, value is a gene structure
    std::map<std::string, Transcript> transcriptHashMap; // key is transcript name, value is a transcript structure
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);

    std::vector<AlignmentMatch> alignmentMatchsMapT0;

    readSam( alignmentMatchsMapT0, infile, transcriptHashMap);
    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;
    for( AlignmentMatch orthologPair2 : alignmentMatchsMapT0 ){
        if(  alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end() ){
            alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
        }
        alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(orthologPair2);
    }
    bool keepTandemDuplication = false;
    for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it=alignmentMatchsMapT.begin(); it != alignmentMatchsMapT.end(); ++it){
        myAlignmentMatchSort( it->second, inversion_PENALTY,  MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);
        std::vector<AlignmentMatch> sortedAlignmentMatchs;
        if( it->second.size()>1 ) {
            longestPath(it->second, sortedAlignmentMatchs, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
        }else{
            sortedAlignmentMatchs = it->second;
        }
        std::vector<unsigned long> keepIndexs;
        for(unsigned long i=0; i<sortedAlignmentMatchs.size(); ++i ){
            if(sortedAlignmentMatchs[i].getStrand() == POSITIVE && sortedAlignmentMatchs[i].getQueryChr().compare(sortedAlignmentMatchs[i].getRefChr())==0){
                if ( keepIndexs.size() == 0 ){
                    keepIndexs.push_back(i);
                }else if ( sortedAlignmentMatchs[keepIndexs[keepIndexs.size()-1]].getQueryEndPos() < sortedAlignmentMatchs[i].getQueryStartPos() &&
                           sortedAlignmentMatchs[keepIndexs[keepIndexs.size()-1]].getRefEndPos() < sortedAlignmentMatchs[i].getRefStartPos() ){
                    keepIndexs.push_back(i);
                }
            }else if(sortedAlignmentMatchs[i].getStrand() == NEGATIVE && sortedAlignmentMatchs[i].getQueryChr().compare(sortedAlignmentMatchs[i].getRefChr())==0 && considerInversion){
                if ( keepIndexs.size() == 0 ){
                    keepIndexs.push_back(i);
                }else if ( ( sortedAlignmentMatchs[keepIndexs[keepIndexs.size()-1]].getRefStartPos() > sortedAlignmentMatchs[i].getRefEndPos() || sortedAlignmentMatchs[i-1].getStrand() == POSITIVE  ) &&
                sortedAlignmentMatchs[keepIndexs[keepIndexs.size()-1]].getQueryEndPos() < sortedAlignmentMatchs[i].getQueryStartPos() ){
                    keepIndexs.push_back(i);
                }
            }
        }
        alignmentMatchsMap[it->first] = std::vector<AlignmentMatch>();
        for (unsigned long j : keepIndexs){
            alignmentMatchsMap[it->first].push_back(sortedAlignmentMatchs[j]);
        }
    }
}


void setupAnchorsWithSpliceAlignmentResultQuota( const std::string & gffFilePath, const std::string & samFile, std::vector<std::vector<AlignmentMatch>> & alignmentMatchsMap,
                                                 double & INDEL_SCORE, double & GAP_OPEN_PENALTY, double & MIN_ALIGNMENT_SCORE,
                                                 int & MAX_DIST_BETWEEN_MATCHES, int & refMaximumTimes, int & queryMaximumTimes,
                                                  double & calculateIndelDistance){
    std::ifstream infile(samFile);
    if( ! infile.good()){
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit (1);
    }

    std::map<std::string, std::vector<std::string> > geneNameMap; // key is chromosome and value is gene names
    std::map<std::string, Gene> geneHashMap;  // key is gene name, value is a gene structure
    std::map<std::string, Transcript> transcriptHashMap; // key is transcript name, value is a transcript structure
    readGffFileWithEveryThing(gffFilePath, geneNameMap, geneHashMap, transcriptHashMap);

    std::vector<AlignmentMatch> alignmentMatchsMapT;

    readSam( alignmentMatchsMapT, infile, transcriptHashMap);

    std::cout << "reading sam file done. size:" << alignmentMatchsMapT.size() << std::endl;
    // begin setting index, they are necessary in the longest path approach
    std::map<std::string, int> queryIndex;
    myOrthologPairsSortQueryQuota( alignmentMatchsMapT);
    for ( size_t ii=0; ii<alignmentMatchsMapT.size() ;++ii){
        if( queryIndex.find( alignmentMatchsMapT[ii].getQueryChr() ) == queryIndex.end() ){
            queryIndex[alignmentMatchsMapT[ii].getQueryChr()] = 0;
        }else{
            queryIndex[alignmentMatchsMapT[ii].getQueryChr()] = queryIndex[alignmentMatchsMapT[ii].getQueryChr()] + 1;
        }
        alignmentMatchsMapT[ii].setQueryId(queryIndex[alignmentMatchsMapT[ii].getQueryChr()]);
    }
    std::cout << "query id setting done"<< std::endl;


    myOrthologPairsSortQuota( alignmentMatchsMapT);
    std::map<std::string, int> refIndex;

    for ( size_t ii=0; ii<alignmentMatchsMapT.size() ;++ii){
        if( refIndex.find(  alignmentMatchsMapT[ii].getRefChr() ) == refIndex.end() ){
            refIndex[alignmentMatchsMapT[ii].getRefChr()] = 0;
        }else{
            refIndex[alignmentMatchsMapT[ii].getRefChr()] = refIndex[alignmentMatchsMapT[ii].getRefChr()] + 1;
        }
        alignmentMatchsMapT[ii].setRefId(refIndex[alignmentMatchsMapT[ii].getRefChr()]);
    }
    // index setting end
    std::cout << "reference id setting done"<< std::endl;


    longestPathQuotav2 (alignmentMatchsMapT, alignmentMatchsMap, INDEL_SCORE, GAP_OPEN_PENALTY,
                      MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES, refMaximumTimes, queryMaximumTimes, calculateIndelDistance);
    for ( size_t i =0; i<alignmentMatchsMap.size(); ++i){
        std::vector<size_t> keepIndexs;
        std::set<size_t> keepIndexset;
        for ( size_t j =0; j<alignmentMatchsMap[i].size(); ++j){
            if(alignmentMatchsMap[i][j].getStrand() == POSITIVE){
                if ( keepIndexs.size() == 0 ){
                    keepIndexs.push_back(j);
                    keepIndexset.insert(j);
                } else if ( alignmentMatchsMap[i][keepIndexs[keepIndexs.size()-1]].getQueryEndPos() < alignmentMatchsMap[i][j].getQueryStartPos() &&
                        alignmentMatchsMap[i][keepIndexs[keepIndexs.size()-1]].getRefEndPos() < alignmentMatchsMap[i][j].getRefStartPos() ){
                    keepIndexs.push_back(j);
                    keepIndexset.insert(j);
                }
            } else {
                if ( keepIndexs.size() == 0 ){
                    keepIndexs.push_back(j);
                    keepIndexset.insert(j);
                } else if ( alignmentMatchsMap[i][keepIndexs[keepIndexs.size()-1]].getQueryStartPos() > alignmentMatchsMap[i][j].getQueryEndPos() &&
                           alignmentMatchsMap[i][keepIndexs[keepIndexs.size()-1]].getRefEndPos() < alignmentMatchsMap[i][j].getRefStartPos() ){
                    keepIndexs.push_back(j);
                    keepIndexset.insert(j);
                }
            }
        }
        std::vector<size_t> toRemoves;
        int this_size = alignmentMatchsMap[i].size();
        for ( int jj = this_size-1; jj>=0; jj-- ) {
            if ( keepIndexset.find(jj) == keepIndexset.end() ){
                alignmentMatchsMap[i].erase(alignmentMatchsMap[i].begin()+jj);
            }
        }
    }
}
