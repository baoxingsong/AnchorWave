//
// Created by baoxing on 10/10/17.
//

#include "getSequencesFromGff.h"
#include "CheckAndUpdateTranscriptsEnds.h"

void getSequences(const std::string& gffFile, const std::string& genomeFile,
                  const std::string& outputCdsSequences, std::map<std::string, std::string>& parameters, const int & minExon){
    std::string regex = get_parameters("cdsParentRegex", parameters);
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix(parameters);

    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    readGffFile (gffFile, transcriptHashSet, regex, minExon);
    std::map<std::string, Fasta> genome;
    readFastaFile(genomeFile, genome);
    CheckAndUpdateTranscriptsEnds( transcriptHashSet, genome, nucleotideCodeSubstitutionMatrix, minExon);

    std::map<std::string, std::string > transcript_to_gene_map;

    get_transcript_to_gene_map_from_gff (gffFile, transcript_to_gene_map);

    std::ofstream oCfile;
    oCfile.open(outputCdsSequences);
    std::set<std::string> usedGenes;

    for( std::map<std::string, std::vector<Transcript> >::iterator it1=transcriptHashSet.begin();
            it1!=transcriptHashSet.end(); ++it1 ){
        if( genome.find(it1->first) != genome.end() ){
            for ( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2 ) {
                if( transcript_to_gene_map.find((*it2).getName()) != transcript_to_gene_map.end() && usedGenes.find( transcript_to_gene_map[(*it2).getName()] ) == usedGenes.end() ){
                    TranscriptUpdateCdsInformation((*it2), genome);
                    //checkOrfState((*it2), genome, nucleotideCodeSubstitutionMatrix, minExon);
                    std::string cdsSequence = (*it2).getCdsSequence();
                    oCfile << ">" << (*it2).getName() << " " << transcript_to_gene_map[(*it2).getName()] << std::endl;
                    oCfile << (*it2).getCdsSequence() << std::endl;
                    usedGenes.insert( transcript_to_gene_map[(*it2).getName()] );
                }
            }
        }
    }
    oCfile.close();
}
