//
// Created by baoxing on 10/10/17.
//

#include "checkOrfState.h"


void checkOrfState( Transcript& targetTranscript, std::map<std::string, Fasta>& targetGenome,
                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix, const int& minIntron){
    std::stringstream metaInformation;
    bool orfShift = false;
    if( ifSpliceSitesOk(targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix) ){
        metaInformation << "spliceSitesConserved";
    }else{
        metaInformation << "spliceSitesDestroyed";
        orfShift = true;
    }
    if( ifIntronEnoughLarge(targetTranscript, minIntron) ){
        metaInformation << "_intronsLargeEnough";
    }else{
        metaInformation << "_intronsNotLargeEnough";
        orfShift = true;
    }
    std::string cdsSequenceString = targetTranscript.getCdsSequence();
    if (cdsSequenceString.length() < 3) {
        metaInformation << "_exonLengthLessThan3";
        orfShift = true;
    } else {
        metaInformation << "_exonLengthMoreThan3";
        if (ifLengthDivisibleByThree(cdsSequenceString)) {
            metaInformation << "_exonLengthIsDivisibleBy3";
        } else {
            metaInformation << "_exonLengthIsNotMultipleOf3";
            orfShift = true;
        }
        if (ifNewStopCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_prematureStopCodon";
            orfShift = true;
        } else {
            metaInformation << "_noPrematureStopCodon";
        }
        if (ifEndWithStopCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_endWithStopCodon";
        } else {
            metaInformation << "_notEndWithStopCodon";
            orfShift = true;
        }
        if (ifStartWithStartCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_startWithStartCodon";
        } else {
            metaInformation << "_notStartWithStartCodon";
            orfShift = true;
        }
    }
    if( !orfShift ){
        metaInformation << "_ConservedFunction";
    }
    metaInformation << " location:" << targetTranscript.getPStart() << "-" << targetTranscript.getPEnd();
    if( POSITIVE ==  targetTranscript.getStrand()){
        metaInformation << "_positive";
    }else{
        metaInformation << "_negative";
    }
    std::string tempMetaInformation = metaInformation.str();
    targetTranscript.setMetaInformation(tempMetaInformation);
    targetTranscript.setIfOrfShift(orfShift);
}

bool checkSpliceSites( std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    s1t = gtIUPACcodesTranslation(s1t);
    s2t = agIUPACcodesTranslation(s2t);
    std::vector<std::string>& dornors = nucleotideCodeSubstitutionMatrix.getDornors();
    std::vector<std::string>& acceptors = nucleotideCodeSubstitutionMatrix.getAcceptors();

    for( size_t size =0; size<dornors.size(); ++size ){
        std::string dornor=dornors[size];
        std::string acceptor=acceptors[size];
        std::set<std::string> allDornors = nucleotideCodeSubstitutionMatrix.getAllPossiableDornors(dornor);
        std::set<std::string> allAcceptors = nucleotideCodeSubstitutionMatrix.getAllPossiableAcceptors(acceptor);
        if( allDornors.find(s1t)!=allDornors.end() && allAcceptors.find(s2t)!=allAcceptors.end() ){
            return true;
        }
    }
//    std::cerr << s1t << " " << s2t <<std::endl;
    return false;
}

bool checkSpliceSites(std::string& s1, std::string& s2, std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    s2 = agIUPACcodesTranslation(s2);
    s1 = gtIUPACcodesTranslation(s1);
    s2t = agIUPACcodesTranslation(s2t);
    s1t = gtIUPACcodesTranslation(s1t);

    if ( checkSpliceSites(s1t, s2t, nucleotideCodeSubstitutionMatrix) ) {
        return true;
    } else {
        std::set<std::string> allS1s;
        nucleotideCodeSubstitutionMatrix.getAllPossibleWithIupac(s1, allS1s);
        if (allS1s.find(s1t) == allS1s.end()) {
            std::cout << "problem splice site: " << s1t << "   " << s2t << std::endl;
            return false;
        }
        std::set<std::string> allS2s;
        nucleotideCodeSubstitutionMatrix.getAllPossibleWithIupac(s2, allS2s);
        if (allS2s.find(s2t) == allS2s.end() ) {
            //std::cout << s1t << "  " << s2t << std::endl;
            return false;
        }
    }
    return true;
}




std::string agIUPACcodesTranslation(std::string& ag) {
    std::string agString = ag;
    if ( agString.length()==2 && ('A' == ag[0] || 'R' == ag[0] || 'W' == ag[0]
                                  || 'M' == ag[0] || 'D' == ag[0]
                                  || 'H' == ag[0] || 'V' == ag[0] || 'N' == ag[0])
         && ('G' == ag[1] || 'K' == ag[1]
             || 'R' == ag[1] || 'S' == ag[1]
             || 'B' == ag[1] || 'D' == ag[1]
             || 'V' == ag[1] || 'N' == ag[1])) {
        agString = "AG";
    }
    return agString;
}

std::string gtIUPACcodesTranslation(std::string& gt) {
    std::string gtString = gt;
    if ( gtString.length()==2 && ('G' == gtString[0] || 'K' == gtString[0]
                                  || 'R' == gtString[0] || 'S' == gtString[0]
                                  || 'B' == gtString[0] || 'D' == gtString[0]
                                  || 'V' == gtString[0] || 'N' == gtString[0])
         && ('T' == gtString[1] || 'K' == gtString[1]
             || 'Y' == gtString[1]
             || 'W' == gtString[1]
             || 'U' == gtString[1]
             || 'B' == gtString[1]
             || 'D' == gtString[1]
             || 'H' == gtString[1] || 'N' == gtString[1])) {
        gtString = "GT";
    }
    return gtString;
}

bool ifLengthDivisibleByThree(std::string& sequence){
    return 0 == sequence.length() % 3;
}

bool ifNewStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    for (size_t j = 0; j < cdsSequence.length() - 3; j += 3) {
        std::string threeNaInFrame = cdsSequence.substr(j,3);
        if (nucleotideCodeSubstitutionMatrix.getMustStopCodons().find(threeNaInFrame) !=
            nucleotideCodeSubstitutionMatrix.getMustStopCodons().end()) {// new stop codon
            //std::cout << "443 " << cdsSequence << " " << threeNaInFrame << " " << j << std::endl;
            return true;
        }
    }
    return false;
}

bool ifEndWithStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    std::string threeNaInFrame = cdsSequence.substr(cdsSequence.length() - 3, 3);
    //std::cout << "472 " << threeNaInFrame << std::endl << cdsSequence << std::endl;
    return nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().find(threeNaInFrame)
           != nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().end();
}

bool ifStartWithStartCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    std::string threeNaInFrame = cdsSequence.substr(0, 3);
    //std::cout << "479 " << threeNaInFrame << std::endl << cdsSequence << std::endl;
    return nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().find(threeNaInFrame)
           != nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().end();
}


bool ifSpliceSitesOk(Transcript& targetTranscript, Transcript &referenceTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     std::map<std::string, Fasta>& referenceGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    bool ifSelectTaur10 = true;

    std::string chrName = targetTranscript.getChromeSomeName();

    for (size_t i = 1; i < referenceTranscript.getCdsVector().size(); i++) {
        int le;
        int ts;
        std::string s1;
        std::string s2;
        int let;
        int tst;
        std::string s1t;
        std::string s2t;
        if (referenceTranscript.getStrand() == POSITIVE) {
            STRAND thisStrand = POSITIVE;

            le = referenceTranscript.getCdsVector()[i-1].getEnd();
            ts = referenceTranscript.getCdsVector()[i].getStart();
            int thisStart01 = le + 1;
            int thisEnd01 = le + 2;
            int thisStart02 = ts - 2;
            int thisEnd02 = ts - 1;
            s1 = getSubsequence(referenceGenome, chrName, thisStart01, thisEnd01, thisStrand);
            s2 = getSubsequence(referenceGenome, chrName, thisStart02, thisEnd02, thisStrand);

            let = targetTranscript.getCdsVector()[i-1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            int thisStart1 = let + 1;
            int thisEnd1 = let + 2;
            int thisStart2 = tst - 2;
            int thisEnd2 = tst - 1;
            s1t = getSubsequence(targetGenome, chrName, thisStart1, thisEnd1, thisStrand);
            s2t = getSubsequence(targetGenome, chrName, thisStart2, thisEnd2, thisStrand);
        } else {
            STRAND thisStrand = NEGATIVE;

            le = referenceTranscript.getCdsVector()[i-1].getEnd();
            ts = referenceTranscript.getCdsVector()[i].getStart();
            int thisStart01 = ts - 2;
            int thisEnd01 = ts - 1;
            int thisStart02 = le + 1;
            int thisEnd02 = le + 2;
            s1 = getSubsequence(referenceGenome, chrName, thisStart01, thisEnd01, thisStrand);
            s2 = getSubsequence(referenceGenome, chrName, thisStart02, thisEnd02, thisStrand);

            let = targetTranscript.getCdsVector()[i - 1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            int thisStart1 = tst - 2;
            int thisEnd1 = tst - 1;
            int thisStart2 = let + 1;
            int thisEnd2 = let + 2;
            s1t = getSubsequence(targetGenome, chrName, thisStart1, thisEnd1, thisStrand);
            s2t = getSubsequence(targetGenome, chrName, thisStart2, thisEnd2, thisStrand);
        }

        if( ! checkSpliceSites(s1, s2, s1t, s2t, nucleotideCodeSubstitutionMatrix)  ){
            return false;
        }
    }
    return ifSelectTaur10;
}

bool ifSpliceSitesOk(Transcript& targetTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {


    std::string chrName = targetTranscript.getChromeSomeName();

    for (size_t i = 1; i < targetTranscript.getCdsVector().size(); i++) {
        int let;
        int tst;
        std::string s1t;
        std::string s2t;

        if (targetTranscript.getStrand() == POSITIVE) {
            let = targetTranscript.getCdsVector()[i-1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            int thisStart1 = let + 1;
            int thisEnd1 = let + 2;

            int thisStart2 = tst - 2;
            int thisEnd2 = tst - 1;
            STRAND thisStrand = POSITIVE;
            s1t = getSubsequence(targetGenome, chrName, thisStart1, thisEnd1, thisStrand);
            s2t = getSubsequence(targetGenome, chrName, thisStart2, thisEnd2, thisStrand);
        } else {
            let = targetTranscript.getCdsVector()[i - 1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            int thisStart1 = tst - 2;
            int thisEnd1 = tst - 1;

            int thisStart2 = let + 1;
            int thisEnd2 = let + 2;
            STRAND thisStrand = NEGATIVE;
            s1t = getSubsequence(targetGenome, chrName, thisStart1, thisEnd1, thisStrand);
            s2t = getSubsequence(targetGenome, chrName, thisStart2, thisEnd2, thisStrand);
            //std::cerr << let << "\t" << tst << "\t" << s1t << "\t" << s2t << std::endl;
        }

        if( ! checkSpliceSites( s1t, s2t, nucleotideCodeSubstitutionMatrix)  ){
            return false;
        }
    }
    return true;
}

bool ifIntronEnoughLarge(Transcript& targetTranscript, const int& minIntron) {
    for (size_t i = 1; i < targetTranscript.getCdsVector().size(); i++) {
        if( (targetTranscript.getCdsVector()[i].getStart() - targetTranscript.getCdsVector()[i-1].getEnd() - 1) < minIntron ){
            return false;
        }
    }
    return true;
}

void checkOrfPversionBin(int const  index, int const & binSize, std::vector<Transcript>   & transcripts, std::map<std::string, Fasta> & genome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,  int const & minIntron, std::atomic_int & number_of_runing_threads){
    int checked = 0;
    for ( size_t i=index; checked<binSize && i<(transcripts.size()); ++checked, ++i ){
        TranscriptUpdateCdsInformation(transcripts[i], genome);
        checkOrfState(transcripts[i], genome, nucleotideCodeSubstitutionMatrix, minIntron);
    }
    --number_of_runing_threads;
}


void checkOrfPversion(Transcript& transcript, std::map<std::string, Fasta>& genome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,  int& minIntron, std::atomic_int & number_of_runing_threads){
    TranscriptUpdateCdsInformation(transcript, genome);
    checkOrfState(transcript, genome, nucleotideCodeSubstitutionMatrix, minIntron);
    --number_of_runing_threads;
}