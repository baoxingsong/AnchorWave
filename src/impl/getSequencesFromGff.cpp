//
// Created by baoxing on 10/10/17.
//

#include "getSequencesFromGff.h"
#include "CheckAndUpdateTranscriptsEnds.h"

void getSequences(const std::string &gffFile, const std::string &genomeFile,
                  const std::string &outputCdsSequences, std::map<std::string, std::string> &parameters, const int &minExon, const bool &exonModel) {
    std::string regex = get_parameters("cdsParentRegex", parameters);
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;

    std::map<std::string, std::string> genome;
    readFastaFile(genomeFile, genome);

    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    if (exonModel) {
        readGffFile_exon(gffFile, transcriptHashSet, regex, minExon);
    } else {
        readGffFile(gffFile, transcriptHashSet, regex, minExon);
    }
    std::set<std::string> toRemoveChrs;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it) {
        if (genome.find(it->first) == genome.end()) {
            toRemoveChrs.insert(it->first);
        }
    }
    for (std::string chr: toRemoveChrs) {
        transcriptHashSet.erase(chr);
    }
    if (exonModel) {

    } else {
        // CheckAndUpdateTranscriptsEnds(transcriptHashSet, genome, nucleotideCodeSubstitutionMatrix);
    }
//    CheckAndUpdateTranscriptsEnds( transcriptHashSet, genome, nucleotideCodeSubstitutionMatrix);

    std::map<std::string, std::string> transcript_to_gene_map;
    get_transcript_to_gene_map_from_gff(gffFile, transcript_to_gene_map);

    std::map<std::string, std::string> usedSeq;
    std::set<std::string> geneBlackList;
    std::map<std::string, std::string> seqToOutPut;

    for (std::map<std::string, std::vector<Transcript> >::iterator it1 = transcriptHashSet.begin(); it1 != transcriptHashSet.end(); ++it1) {
        if (genome.find(it1->first) != genome.end()) {
            for (std::vector<Transcript>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
                if (transcript_to_gene_map.find((*it2).getName()) != transcript_to_gene_map.end()) {
                    TranscriptUpdateCdsInformation((*it2), genome);
                    std::string cdsSequence = (*it2).getCdsSequence();
                    std::string transcriptName = (*it2).getName();
                    if (usedSeq.find(cdsSequence) == usedSeq.end()) {
                        if (seqToOutPut.find(transcript_to_gene_map[(*it2).getName()]) == seqToOutPut.end()) {
                            seqToOutPut[transcript_to_gene_map[(*it2).getName()]] = cdsSequence;
                        } else if (seqToOutPut[transcript_to_gene_map[(*it2).getName()]].length() < cdsSequence.length()) {
                            seqToOutPut[transcript_to_gene_map[(*it2).getName()]] = cdsSequence;
                        }
                        usedSeq[cdsSequence] = transcriptName;
                    } else {
                        geneBlackList.insert(transcript_to_gene_map[(*it2).getName()]);
                        geneBlackList.insert(transcript_to_gene_map[usedSeq[cdsSequence]]);
                    }
                }
            }
        }
    }

    std::ofstream oCfile;
    oCfile.open(outputCdsSequences);
    for (std::map<std::string, std::string>::iterator it = seqToOutPut.begin(); it != seqToOutPut.end(); ++it) {
        if (geneBlackList.find(it->first) == geneBlackList.end()) {
            oCfile << ">" << usedSeq[it->second] << " " << it->first << std::endl;
            oCfile << it->second << std::endl;
        }
    }
    oCfile.close();
}
