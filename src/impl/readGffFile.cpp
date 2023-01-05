//
// Created by baoxing on 10/10/17.
//

#include "readGffFile.h"

void readGffFile(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const int &minExon) {
    std::string cdsParentRegex = "([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
    readGffFile(filePath, transcriptHashSet, cdsParentRegex, minExon);
}

int32_t gffmin(const int32_t &a, const int32_t &b) {
    if (a < b) {
        return a;
    }
    return b;
}
void get_transcript_to_gene_map_from_gff(const std::string &filePath, std::map <std::string, std::string> &transcript_to_gene_map) {
    std::vector <std::string> transcriptParentRegex;
    transcriptParentRegex.push_back("ID=(\\S+?);.*Parent=(\\S+?);");
    transcriptParentRegex.push_back("ID=(\\S+?);.*Parent=(\\S+?)$");
    transcriptParentRegex.push_back("ID=(\\S+?);.*geneID=(\\S+?)$");

    std::vector <std::regex> regTranscriptParents;
    for (std::string transcript: transcriptParentRegex) {
        std::regex regTranscript(transcript);
        regTranscriptParents.push_back(regTranscript);
    }

    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line[0] == '#') {
            continue;
        }
        std::string subline = line.substr(0, gffmin(10000, line.size()));
        std::smatch match;
        for (size_t i = 0; i < regTranscriptParents.size(); i++) {
            std::regex reg = regTranscriptParents[i];
            while (regex_search(subline, match, reg)) {
                std::string transcript_id = match[1];
                std::string gene_id = match[2];
                transcript_to_gene_map[transcript_id] = gene_id;

                if (i == 0) {
                    line = match.suffix();
                } else {
                    break;
                }
            }
        }
    }

    while (std::getline(infile, line)) {
        if (line[0] == '#') {
            continue;
        }

        std::smatch match;
        std::regex reg("Parent=(\\S+?);ID=(\\S+?);");
        std::string subline = line.substr(0, gffmin(10000, line.size()));
        while (regex_search(subline, match, reg)) {
            std::string transcript_id = match[2];
            std::string gene_id = match[1];
            transcript_to_gene_map[transcript_id] = gene_id;
            line = match.suffix();
        }

        std::regex reg2("Parent=(\\S+?);ID=(\\S+?)$");
        while (regex_search(line, match, reg2)) {
            std::string transcript_id = match[2];
            std::string gene_id = match[1];
            transcript_to_gene_map[transcript_id] = gene_id;
            line = match.suffix();
        }
    }
}


void readGffFile(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
    std::map <std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::regex reg("^(\\S*)\t([\\s\\S]*)\tCDS\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t" + cdsParentRegex);
    std::string line;
    while (std::getline(infile, line)) {
        std::smatch match;
        std::string subline = line.substr(0, gffmin(10000, line.size()));
        regex_search(subline, match, reg);
        //std::cout << line << std::endl;
        if (match.empty() || line[0] == '#' || line.size() < 9) {
        } else {
            //std::cout << line << std::endl;
            int start = stoi(match[3]);
            int end = stoi(match[4]);
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }
            if ((end - start + 1) >= minExon) {
                std::string information = match[9];
                //std::cout << information << std::endl;
                if (transcriptHashMap.find(information) != transcriptHashMap.end()) {
                } else {
                    std::string chromosomeName = match[1];
                    STRAND strand;
                    if (match[6].compare("-") == 0) {
                        strand = NEGATIVE;
                    } else {
                        strand = POSITIVE;
                    }
                    Transcript transcript1(information, chromosomeName, strand);
                    transcriptHashMap[information] = transcript1;
                }
                GenomeBasicFeature cds(start, end);
                //cds.setTranscript(transcriptHashMap[information]);
                transcriptHashMap[information].addCds(cds);
            }
        }
    }

    for (std::map<std::string, Transcript>::iterator it = transcriptHashMap.begin(); it != transcriptHashMap.end(); ++it) {
        if (transcriptHashSet.find(it->second.getChromeSomeName()) == transcriptHashSet.end()) {
            transcriptHashSet[it->second.getChromeSomeName()] = std::vector<Transcript>();
        }

        it->second.updateInforCds();
        transcriptHashSet[it->second.getChromeSomeName()].push_back(it->second);
    }

    for (std::map < std::string, std::vector < Transcript >> ::iterator it = transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](Transcript a, Transcript b) {
            return a.getPStart() < b.getPStart();
        });
    }
}

void readGffFile_exon(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
    std::map <std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }
    std::regex reg("^(\\S*)\t([\\s\\S]*)\texon\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t" + cdsParentRegex);
    std::string line;
    while (std::getline(infile, line)) {
        std::smatch match;
        std::string subline = line.substr(0, gffmin(10000, line.size()));
        regex_search(subline, match, reg);

        if (match.empty() || line[0] == '#' || line.size() < 9) {
        } else {
            int start = stoi(match[3]);
            int end = stoi(match[4]);
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }
            //std::cout << start << "\t" << end << std::endl;
            if ((end - start + 1) >= minExon) {
                std::string information = match[9];
                if (transcriptHashMap.find(information) != transcriptHashMap.end()) {
                } else {
                    std::string chromosomeName = match[1];
                    STRAND strand;
                    if (match[6].compare("-") == 0) {
                        strand = NEGATIVE;
                    } else {
                        strand = POSITIVE;
                    }
                    Transcript transcript1(information, chromosomeName, strand);
                    transcriptHashMap[information] = transcript1;
                }
                GenomeBasicFeature cds(start, end);
                //cds.setTranscript(transcriptHashMap[information]);
                transcriptHashMap[information].addCds(cds);
//                std::cout << information << "\t" << start << "\t" << end << std::endl;
            }
        }
    }
    infile.close();

    for (std::map<std::string, Transcript>::iterator it = transcriptHashMap.begin(); it != transcriptHashMap.end(); ++it) {
        if (transcriptHashSet.find(it->second.getChromeSomeName()) == transcriptHashSet.end()) {
            transcriptHashSet[it->second.getChromeSomeName()] = std::vector<Transcript>();
        }

        it->second.updateInforCds();
        transcriptHashSet[it->second.getChromeSomeName()].push_back(it->second);
    }

    for (std::map < std::string, std::vector < Transcript >> ::iterator it = transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](Transcript a, Transcript b) {
            return a.getPStart() < b.getPStart();
        });
    }
}
