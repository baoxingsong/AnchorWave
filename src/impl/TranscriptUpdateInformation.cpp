//
// Created by baoxing on 10/10/17.
//

#include "TranscriptUpdateInformation.h"

void TranscriptUpdateCdsInformation(Transcript & transcript, std::map<std::string, Fasta>& genome) {
    if(transcript.getCdsVector().size()>0){
        transcript.updateInforCds();

        STRAND strand = transcript.getStrand();

        std::string genomeSequence = getSubsequence(genome, transcript.getChromeSomeName(),  transcript.getPStart(), transcript.getPEnd(), strand);
        transcript.setGeneomeSequence( genomeSequence );

        std::string reserveString;
        reserveString.reserve(transcript.getGeneomeSequence().size());
        std::stringstream cdsss(reserveString);
        if( POSITIVE == strand ){
            for( size_t i=0; i < transcript.getCdsVector().size(); i++ ){
                int start = transcript.getCdsVector()[i].getStart();
                int end = transcript.getCdsVector()[i].getEnd();
                cdsss << getSubsequence(genome, transcript.getChromeSomeName(), start, end, strand);
            }
        }else{
            for( size_t i=transcript.getCdsVector().size(); i>0; i-- ){
                int start = transcript.getCdsVector()[i-1].getStart();
                int end = transcript.getCdsVector()[i-1].getEnd();
                cdsss << getSubsequence(genome, transcript.getChromeSomeName(), start, end, strand);
            }
        }
        std::string cdsSequence = cdsss.str();
        transcript.setCdsSequence( cdsSequence );
    }
}
