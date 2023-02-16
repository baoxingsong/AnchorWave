//
// Created by baoxing on 10/10/17.
//

#include "TranscriptUpdateInformation.h"

void TranscriptUpdateCdsInformation(Transcript &transcript, std::map<std::string, std::tuple<std::string, long, long, int> > &genome) {
    if (transcript.getCdsVector().size() > 0) {
        transcript.updateInforCds();

        STRAND strand = transcript.getStrand();

        std::string sq = getSubsequence2(genome, transcript.getChromeSomeName(), transcript.getPStart(), transcript.getPEnd(), strand);
        transcript.setGeneomeSequence(sq);

        std::string reserveString;
        reserveString.reserve(transcript.getGeneomeSequence().size());
        std::stringstream ss(reserveString);

        if (POSITIVE == strand) {
            for (size_t i = 0; i < transcript.getCdsVector().size(); i++) {
                int start = transcript.getCdsVector()[i].getStart();
                int end = transcript.getCdsVector()[i].getEnd();

                std::string s2 = getSubsequence2(genome, transcript.getChromeSomeName(), start, end, strand);
                ss << s2;
            }
        } else {
            for (size_t i = transcript.getCdsVector().size(); i > 0; i--) {
                int start = transcript.getCdsVector()[i - 1].getStart();
                int end = transcript.getCdsVector()[i - 1].getEnd();

                std::string s2 = getSubsequence2(genome, transcript.getChromeSomeName(), start, end, strand);
                ss << s2;
            }
        }

        std::string cdsSequence = ss.str();
        transcript.setCdsSequence(cdsSequence);
    }
}
