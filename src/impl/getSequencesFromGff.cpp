//
// Created by baoxing on 10/10/17.
//

#include "getSequencesFromGff.h"

void getSequences(const std::string &gffFile, const std::string &genomeFile, const std::string &outputCdsSequences, const int &minExon, const bool &exonModel) {

    std::map<std::string, std::tuple<std::string, long, long, int> > genome;
    readFastaFile(genomeFile, genome);

    std::map<std::string, std::vector<Transcript> > map_v_ts;
    if (exonModel) {
        readGffFile(gffFile, map_v_ts, "exon", minExon);
    } else {
        readGffFile(gffFile, map_v_ts, "CDS", minExon);
    }

    std::set<std::string> set_rm_chr;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = map_v_ts.begin(); it != map_v_ts.end(); ++it) {
        if (genome.find(it->first) == genome.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::string chr: set_rm_chr) {
        map_v_ts.erase(chr);
    }

    std::map<std::string, std::string> map_ts_to_gene; // map: transcript to gene
    get_map_from_gff(gffFile, map_ts_to_gene);

    std::map<std::string, std::string> map_used;
    std::set<std::string> geneBlackList;
    std::map<std::string, std::string> seqToOutPut;

    for (std::map<std::string, std::vector<Transcript> >::iterator it1 = map_v_ts.begin(); it1 != map_v_ts.end(); ++it1) {
        if (genome.find(it1->first) != genome.end()) {
            for (std::vector<Transcript>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
                if (map_ts_to_gene.find((*it2).getName()) != map_ts_to_gene.end()) {
                    TranscriptUpdateCdsInformation((*it2), genome);
                    std::string seq_cds = (*it2).getCdsSequence();
                    std::string name_transcript = (*it2).getName();

                    if (map_used.find(seq_cds) == map_used.end()) {
                        if (seqToOutPut.find(map_ts_to_gene[(*it2).getName()]) == seqToOutPut.end()) {
                            seqToOutPut[map_ts_to_gene[(*it2).getName()]] = seq_cds;
                        }
                        else if (seqToOutPut[map_ts_to_gene[(*it2).getName()]].length() < seq_cds.length()) {
                            seqToOutPut[map_ts_to_gene[(*it2).getName()]] = seq_cds;
                        }

                        map_used[seq_cds] = name_transcript;
                    } else {
                        geneBlackList.insert(map_ts_to_gene[(*it2).getName()]);
                        geneBlackList.insert(map_ts_to_gene[map_used[seq_cds]]);
                    }
                }
            }
        }
    }

    std::ofstream oCfile;
    oCfile.open(outputCdsSequences);
    for (std::map<std::string, std::string>::iterator it = seqToOutPut.begin(); it != seqToOutPut.end(); ++it) {
        if (geneBlackList.find(it->first) == geneBlackList.end()) {
            oCfile << ">" << map_used[it->second] << " " << it->first << std::endl;
            oCfile << it->second << std::endl; // does it work?
        }
    }
    oCfile.close();
}
