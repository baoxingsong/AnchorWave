//
// Created by Baoxing song on 08.08.18.
// this function here does not read real GFF correctly, it ignores short CDS.
// It is designed to do so.
//

#include "readGffFileWithEverything.h"
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, std::vector<std::string> > & geneNameMap,
                                std::map<std::string, Gene > & geneHashMap,
                                std::map<std::string, Transcript> & transcriptHashMap, const int & minExon){
    std::set<std::string> cdsParentRegex;
    std::set<std::string> transcriptParentRegex;
    std::set<std::string> transcriptIdRegex;
    std::set<std::string> geneIdRegex;
    //tair10 begin
    cdsParentRegex.insert("Parent=([\\s\\S]*?)[;,]");//CDS
    cdsParentRegex.insert("Parent=([\\s\\S^,^;]*?)$");//exon
    transcriptParentRegex.insert("Parent=([\\s\\S]*?)[;,]");
    transcriptIdRegex.insert("ID=([\\s\\S]*?)[;,]");
    //geneIdRegex.insert("ID=([\\s\\S]*?)[;,;]");
    geneIdRegex.insert("ID=([-_0-9:a-zA-Z.]*?)[;,]");
    geneIdRegex.insert("ID=([-_0-9:a-zA-Z.]*?)$"); //melon
    // tair10 end

    //ensembl begin
    cdsParentRegex.insert("transcript_id\\s*\"([\\s\\S]*?)\"[;,]");
    transcriptParentRegex.insert("gene_id\\s*\"([\\s\\S]*?)\"[;,]");
    transcriptIdRegex.insert("transcript_id\\s*\"([\\s\\S]*?)\"[;,]");
    geneIdRegex.insert("gene_id\\s*\"([\\s\\S]*?)\"[;,]");
    //ensembl end

    /* //begin
    cdsParentRegex.insert("transcript_id \"([\\s\\S]*?)\"[;,]");
    transcriptParentRegex.insert("");
    transcriptIdRegex.insert();
    geneIdRegex.insert();
    end */

    transcriptParentRegex.insert("Parent=([-_0-9:a-zA-Z.]*?)$");

    std::set<std::string> transcriptNickNames;
    transcriptNickNames.insert("transcript");
    transcriptNickNames.insert("pseudogenic_transcript");
    transcriptNickNames.insert("mRNA");
    transcriptNickNames.insert("miRNA");
    transcriptNickNames.insert("tRNA");
    transcriptNickNames.insert("ncRNA");
    transcriptNickNames.insert("mRNA_TE_gene");
    transcriptNickNames.insert("rRNA");
    transcriptNickNames.insert("snoRNA");
    transcriptNickNames.insert("snRNA");
    transcriptNickNames.insert("lincRNA");

    std::set<std::string> geneNickNames;
    geneNickNames.insert("gene");
    geneNickNames.insert("pseudogene");
    geneNickNames.insert("transposable_element_gene");
    geneNickNames.insert("lincRNA_gene");
    geneNickNames.insert("tRNA_gene");
    /*

    std::set<std::string> transcriptNickNames;
    transcriptNickNames.insert("transcript");
    transcriptNickNames.insert("pseudogenic_transcript");
    transcriptNickNames.insert("mRNA");
    transcriptNickNames.insert("miRNA");
    transcriptNickNames.insert("tRNA");
    transcriptNickNames.insert("ncRNA");
    transcriptNickNames.insert("mRNA_TE_gene");
    transcriptNickNames.insert("rRNA");
    transcriptNickNames.insert("snoRNA");
    transcriptNickNames.insert("snRNA");
    transcriptNickNames.insert("lincRNA");

    std::set<std::string> geneNickNames;
    geneNickNames.insert("gene");
    geneNickNames.insert("pseudogene");
    geneNickNames.insert("transposable_element_gene");
    geneNickNames.insert("tRNA_gene");
    geneNickNames.insert("lincRNA_gene");
    geneNickNames.insert("miRNA_gene");
*/
    std::set<std::string> ignoreTypes;
//    ignoreTypes.insert("pseudogene");
//    ignoreTypes.insert("pseudogenic_transcript");
//    ignoreTypes.insert("pseudogenic_exon");
//    ignoreTypes.insert("transposable_element_gene");
//    ignoreTypes.insert("mRNA_TE_gene");
    ignoreTypes.insert("start_codon");
    ignoreTypes.insert("stop_codon");
    ignoreTypes.insert("protein");
    ignoreTypes.insert("chromosome");
    ignoreTypes.insert("transposable_element");
    ignoreTypes.insert("transposon_fragment");
    ignoreTypes.insert("contig");
    ignoreTypes.insert("miRNA_gene");
    readGffFileWithEveryThing (filePath, geneNameMap, geneHashMap, transcriptHashMap,
                               cdsParentRegex, transcriptParentRegex, transcriptIdRegex, geneIdRegex,
                               transcriptNickNames, geneNickNames, ignoreTypes, minExon);
}

void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, std::vector<std::string> > & geneNameMap,
                                std::map<std::string, Gene > & geneHashMap,
                                std::map<std::string, Transcript> & transcriptHashMap,
        std::set<std::string> & cdsParentRegex, std::set<std::string> & transcriptParentRegex,
        std::set<std::string> & transcriptIdRegex, std::set<std::string> & geneIdRegex, std::set<std::string> & transcriptNickNames,
        std::set<std::string> & geneNickNames, std::set<std::string> & ignoreTypes, const int & minExon){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    std::vector<std::regex> regCdsParents;
    for( std::string cds : cdsParentRegex ){
        std::regex regCdsParent(cds);
        regCdsParents.push_back(regCdsParent);
    }
    std::vector<std::regex> regTranscriptIds;
    for( std::string transcript : transcriptIdRegex ){
        std::regex regTranscript(transcript);
        regTranscriptIds.push_back(regTranscript);
    }
    std::vector<std::regex> regTranscriptParents;
    for( std::string transcript : transcriptParentRegex ){
        std::regex regTranscript(transcript);
        regTranscriptParents.push_back(regTranscript);
    }
    std::vector<std::regex> regGeneIds;
    for( std::string gene : geneIdRegex ){
        std::regex regGene(gene);
        regGeneIds.push_back(regGene);
    }
    char splim='\t';
    std::string line;
    std::smatch match;
    std::string chromosomeName;
    std::string score;
    std::string source;
    std::string lastColumnInformation;
    bool matched;
    while (std::getline(infile, line)){
        if(line[0]=='#' || line.size()<9){
        }else {
            std::vector<std::string> elemetns;
            split(line, splim, elemetns);
            if(elemetns.size()==9){
                int start = stoi(elemetns[3]);
                int end = stoi(elemetns[4]);
                if (start > end) {
                    int temp = start;
                    start = end;
                    end = temp;
                }
                STRAND strand;
                if (elemetns[6].compare("-") == 0) {
                    strand = NEGATIVE;
                } else {
                    strand = POSITIVE;
                }
                chromosomeName = elemetns[0];
                source = elemetns[1];
                score = elemetns[5];
                lastColumnInformation = elemetns[8];
                if (geneNickNames.find(elemetns[2]) != geneNickNames.end()) {
                    //                std::cout << "line 119" << std::endl;
                    matched = false;
                    for (std::regex regGene : regGeneIds) {
                        regex_search(lastColumnInformation, match, regGene);
                        if (!match.empty()) {
                            std::string geneId = match[1];
                            if (geneHashMap.find(geneId) != geneHashMap.end()) {
                            } else {
                                Gene gene1(geneId, chromosomeName, strand);
                                geneHashMap[geneId] = gene1;
                            }
//                            std::cout << "line 170 " << geneId << std::endl;
                            geneHashMap[geneId].setStart(start);
                            geneHashMap[geneId].setEnd(end);
                            geneHashMap[geneId].setSource(source);
                            geneHashMap[geneId].setScore(score);
                            geneHashMap[geneId].setLastColumnInformation(lastColumnInformation);
                            matched = true;
                            break; //jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout
                                << "the gene record does not fellow the regular expression provided for finding gene ID, please check"
                                << std::endl << line << std::endl;
                    }
                } else if (transcriptNickNames.find(elemetns[2]) != transcriptNickNames.end()) {
                    //                std::cout << "line 141" << std::endl;
                    matched = false;
                    std::string transcriptId = "";
                    for (std::regex regTranscriptId : regTranscriptIds) {
                        regex_search(lastColumnInformation, match, regTranscriptId);
                        if (!match.empty()) {
                            transcriptId = match[1];
                            matched = true;
                            break; //jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout
                                << "the transcript record does not fellow the regular expression provided for finding transcript ID, please check"
                                << std::endl << line << std::endl;
                        continue; //goto next line in the gff file
                    }
                    //                std::cout << "line 156" << std::endl;
                    matched = false;
                    for (std::regex regTranscriptParent : regTranscriptParents) {
                        regex_search(lastColumnInformation, match, regTranscriptParent);
                        if (!match.empty()) {
                            std::string geneId = match[1];
                            if (transcriptHashMap.find(transcriptId) != transcriptHashMap.end()) {
                            } else {
                                Transcript transcript1(transcriptId, elemetns[0], strand);
                                transcriptHashMap[transcriptId] = transcript1;
                            }
                            transcriptHashMap[transcriptId].setStart(start);
                            transcriptHashMap[transcriptId].setEnd(end);
                            transcriptHashMap[transcriptId].setSource(source);
                            transcriptHashMap[transcriptId].setScore(score);
                            transcriptHashMap[transcriptId].setLastColumnInformation(lastColumnInformation);
                            transcriptHashMap[transcriptId].setType(elemetns[2]);
                            //std::cout << "transcript  " << transcriptId << std::endl;
                            if (geneHashMap.find(geneId) != geneHashMap.end()) {

                            } else {
                                //std::cout << "there is something wrong with the gff file, maybe could not generate correct result" << std::endl;
                                Gene gene1(geneId, chromosomeName, strand);
                                geneHashMap[geneId] = gene1;
                            }
                            geneHashMap[geneId].addTranscript(transcriptHashMap[transcriptId]);
                            //std::cout << "line 178 very good very good " << geneHashMap[geneId].getName() << " " << transcriptHashMap[transcriptId].getName() << std::endl;
                            matched = true;
                            break; // jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout
                                << "the transcript record does not fellow the regular expression provided for finding gene ID, please check"
                                << std::endl << line << std::endl;
                    }
                } else if (elemetns[2].compare("CDS") == 0 && (end-start+1) >= minExon ) {
                    //                std::cout << "line 188" << std::endl;
                    matched = false;
                    for (std::regex regCdsParent : regCdsParents) {
                        regex_search(lastColumnInformation, match, regCdsParent);
                        if (!match.empty()) {
                            std::string transcriptId = match[1];
                            if (transcriptHashMap.find(transcriptId) != transcriptHashMap.end()) {
                            } else {
                                Transcript transcript1(transcriptId, chromosomeName, strand);
                                transcript1.setSource(elemetns[3]);
                                transcriptHashMap[transcriptId] = transcript1;
                            }
                            GenomeBasicFeature cds(start, end, score, elemetns[7], lastColumnInformation);
                            cds.setType(elemetns[2]);
                            transcriptHashMap[transcriptId].addCds(cds);
                            //                        std::cout << transcriptId << "adding cds" << std::endl;
                            matched = true;
                            break; // jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout
                                << "the CDS record does not fellow the regular expression provided, please check"
                                << std::endl << line << std::endl;
                    }
                } else if (elemetns[2].compare("exon") == 0 || elemetns[2].compare("pseudogenic_exon") == 0) {
                    //                std::cout << "line 208 " << line << std::endl;
                    matched = false;
                    for (std::regex regCdsParent : regCdsParents) {
                        regex_search(lastColumnInformation, match, regCdsParent);
                        if (!match.empty()) {
                            std::string transcriptId = match[1];
                            if (transcriptHashMap.find(transcriptId) != transcriptHashMap.end()) {
                            } else {
                                Transcript transcript1(transcriptId, chromosomeName, strand);
                                transcript1.setSource(source);
                                transcriptHashMap[transcriptId] = transcript1;
                            }
                            GenomeBasicFeature exon(start, end, score, elemetns[7], lastColumnInformation);
                            exon.setType(elemetns[2]);
                            transcriptHashMap[transcriptId].addExon(exon);
                            //                        std::cout << transcriptId << "adding exon" << std::endl;
                            matched = true;
                            break; // jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout << "the exon record does not fellow the regular expression provided, please check"
                                  << std::endl << line << std::endl;
                    }
                } else if (elemetns[2].compare("five_prime_utr") == 0 || elemetns[2].compare("five_prime_UTR") == 0) {
                    //                std::cout << "line 231" << std::endl;
                    matched = false;
                    for (std::regex regCdsParent : regCdsParents) {
                        regex_search(lastColumnInformation, match, regCdsParent);
                        if (!match.empty()) {
                            std::string transcriptId = match[1];
                            if (transcriptHashMap.find(transcriptId) != transcriptHashMap.end()) {
                            } else {
                                Transcript transcript1(transcriptId, chromosomeName, strand);
                                transcript1.setSource(source);
                                transcriptHashMap[transcriptId] = transcript1;
                            }
                            GenomeBasicFeature five_prime_utr(start, end, score, elemetns[7], lastColumnInformation);
                            five_prime_utr.setType(elemetns[2]);
                            transcriptHashMap[transcriptId].addFivePrimerUtr(five_prime_utr);
                            matched = true;
                            break; // jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout
                                << "the five_prime_utr record does not fellow the regular expression provided, please check"
                                << std::endl << line << std::endl;
                    }
                } else if (elemetns[2].compare("three_prime_utr") == 0 || elemetns[2].compare("three_prime_UTR") == 0) {
                    //                std::cout << "line 255" << std::endl;
                    matched = false;
                    for (std::regex regCdsParent : regCdsParents) {
                        regex_search(lastColumnInformation, match, regCdsParent);
                        if (!match.empty()) {
                            std::string transcriptId = match[1];
                            if (transcriptHashMap.find(transcriptId) != transcriptHashMap.end()) {
                            } else {
                                Transcript transcript1(transcriptId, chromosomeName, strand);
                                transcript1.setSource(source);
                                transcriptHashMap[transcriptId] = transcript1;
                            }
                            GenomeBasicFeature three_prime_utr(start, end, score, elemetns[7], lastColumnInformation);
                            three_prime_utr.setType(elemetns[2]);
                            transcriptHashMap[transcriptId].addThreePrimerUtr(three_prime_utr);
                            matched = true;
                            break; // jump out for loop
                        }
                    }
                    if (!matched) {
                        std::cout
                                << "the three_prime_utr record does not fellow the regular expression provided, please check"
                                << std::endl << line << std::endl;
                    }
                } else if (ignoreTypes.find(elemetns[2]) != ignoreTypes.end()) { // ignore those elements

                } else {
//                    std::cout << "we could not analysis the line in the gff/gtf file: " << line << std::endl;
                }
            }
        }
    }
    std::map<std::string, std::vector<Gene>> geneHashSet;
    for (std::map<std::string, Gene>::iterator it=geneHashMap.begin(); it!=geneHashMap.end(); ++it){
        if(geneHashSet.find(it->second.getChromeSomeName()) == geneHashSet.end() ){
            geneHashSet[it->second.getChromeSomeName()]=std::vector<Gene>();
        }
//        std::cout << it->second.getName() << "\t" << it->second.getTranscriptVector().size() << "\t" << std::endl;
        for( std::vector<std::string>::iterator transcript= it->second.getTranscriptVector().begin();
             transcript != it->second.getTranscriptVector().end();
             ++transcript){
            if( transcriptHashMap[*transcript].getCdsVector().size()>0 ){
                transcriptHashMap[*transcript].updateInforCds();
            }
        }
        geneHashSet[it->second.getChromeSomeName()].push_back(it->second);
        //geneNameMap[it->second.getChromeSomeName()].push_back(it->second.getName());
    }
    // the purpose of geneHashSet if to keep the IDs in geneNameMap in order
    // this is very important for the gff record transformation
    for (std::map<std::string, std::vector<Gene>>::iterator it=geneHashSet.begin(); it!=geneHashSet.end(); ++it){
        std::sort(it->second.begin(), it->second.end(), [](Gene a, Gene b) {
            return a.getStart() < b.getStart();
        });
        if(geneNameMap.find(it->first) == geneNameMap.end() ){
            geneNameMap[it->first]=std::vector<std::string>();
        }
        for( Gene gene : it->second ){
            geneNameMap[it->first].push_back(gene.getName());
        }
    }
}
