//
// Created by song on 8/19/18.
//

#include "CheckAndUpdateTranscriptsEnds.h"


void updateTranscriptEnds(Transcript & transcript){
    if( POSITIVE == transcript.getStrand() ){
        if( transcript.getCdsVector().size()>0 ){
            transcript.setPEnd(transcript.getPEnd() + 3);
            transcript.getCdsVector()[transcript.getCdsVector().size()-1].setEnd(transcript.getPEnd() + 3);
        }
        if( transcript.getThreePrimerUtr().size() > 0 ){
            transcript.getThreePrimerUtr()[0].setStart(transcript.getThreePrimerUtr()[0].getStart()+3);
        }
    } else {
        if( transcript.getCdsVector().size()>0 ) {
            transcript.setPStart(transcript.getPStart() - 3);
            transcript.getCdsVector()[0].setStart(transcript.getPStart() - 3);
        }
        if( transcript.getThreePrimerUtr().size() > 0 ){
            transcript.getThreePrimerUtr()[transcript.getThreePrimerUtr().size()-1].setEnd(transcript.getThreePrimerUtr()[transcript.getThreePrimerUtr().size()-1].getEnd()-3);
        }
    }
}

void CheckAndUpdateTranscriptsEnds(std::map<std::string, Transcript> & transcripts, std::map<std::string, Fasta>& sequences,
     NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix, const int & minIntron){
    size_t checkedNumber = 0;
    size_t caseNumber = 0;
    for( std::map<std::string, Transcript>::iterator it = transcripts.begin();
        it != transcripts.end(); ++it){
        checkedNumber++;
        TranscriptUpdateCdsInformation(it->second, sequences);
        checkOrfState( it->second, sequences, nucleotideCodeSubstitutionMatrix, minIntron);
        if( it->second.getIfOrfShift() && it->second.getMetaInformation().find("notEndWithStopCodon")!=std::string::npos ){
            ++caseNumber;
        }
        if( checkedNumber >= 100 ){
            break;
        }
    }
    if( caseNumber > 95 ){
        for( std::map<std::string, Transcript>::iterator it = transcripts.begin();
             it != transcripts.end(); ++it){
            updateTranscriptEnds(it->second);
        }
    }
}

void CheckAndUpdateTranscriptsEnds(std::map<std::string, std::vector<Transcript> > & transcriptHashSet, std::map<std::string, Fasta>& sequences,
                                   NucleotideCodeSubstitutionMatrix & nucleotideCodeSubstitutionMatrix, const int & minIntron){
    size_t checkedNumber = 0;
    size_t caseNumber = 0;
    for( std::map<std::string, std::vector<Transcript>>::iterator it = transcriptHashSet.begin();
         it != transcriptHashSet.end(); ++it){
        for( std::vector<Transcript>::iterator it1=it->second.begin(); it1!=it->second.end(); ++it1 ){
            checkedNumber++;
            TranscriptUpdateCdsInformation(*it1, sequences);
            checkOrfState( *it1, sequences, nucleotideCodeSubstitutionMatrix, minIntron);
            if( it1->getIfOrfShift() && it1->getMetaInformation().find("notEndWithStopCodon")!=std::string::npos ){
                ++caseNumber;
            }
            if( checkedNumber >= 100 ){
                break;
            }
        }
    }
    if( caseNumber > 95 ){
        for( std::map<std::string, std::vector<Transcript>>::iterator it = transcriptHashSet.begin();
             it != transcriptHashSet.end(); ++it){
            for( std::vector<Transcript>::iterator it1=it->second.begin(); it1!=it->second.end(); ++it1 ){
                updateTranscriptEnds(*it1);
            }
        }
    }
}
