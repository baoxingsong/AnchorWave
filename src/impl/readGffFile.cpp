//
// Created by baoxing on 10/10/17.
//

#include "readGffFile.h"

void get_map_from_gff(const std::string &filePath, std::map <std::string, std::string> &map_transcript_to_gene) {
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.size() < 24 || line[0] == '#') {
            continue;
        }

        std::string id;
        size_t p_id_b = line.find("ID=");
        if (p_id_b != std::string::npos) {
            std::string line_s = line.substr(p_id_b + 3);
            size_t p_id_e = line_s.find(";");
            id = line_s.substr(0, p_id_e);
        } else {
            continue;
        }

        size_t p_p_b = line.find("Parent=");
        if (p_p_b != std::string::npos) {
            std::string line_s = line.substr(p_p_b + 7);
            size_t p_p_e = line_s.find(";");
            std::string parent = line_s.substr(0, p_p_e);
            map_transcript_to_gene[id] = parent;
            continue;
        }

        size_t p_g_b = line.find("geneID=");
        if (p_g_b != std::string::npos) {
            std::string line_s = line.substr(p_g_b + 7);
            size_t p_g_e = line_s.find(";");
            std::string gene_id = line_s.substr(0, p_g_e);
            map_transcript_to_gene[id] = gene_id;
        }
    }

    infile.close();
}

void readGffFile(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &type, const int &minExon) {
    std::map <std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // --line like :    chr1	NAM	CDS	34722	35318	.	+	0	ID=Zm00001eb000010_P001;Parent=Zm00001eb000010_T001;protein_id=Zm00001eb000010_P001
        size_t pos_type = line.find(type); // type like "CDS","exon" ...
        if (pos_type == std::string::npos || line.size() < 24 || line[0] == '#') {
            continue;
        }

        int size = line.size();
        char name[size]; // name
        char c_3[size];  // type
        int start, end;  // start, end
        char c_6[1];     // strand
        char c_8[size];  // information

        int ret = sscanf(line.c_str(), "%s%*s%s%d%d%*s%s%*s%s", name, c_3, &start, &end, c_6, c_8);
        if (ret == 6 && std::string(c_3) == type) {
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }

            if ((end - start + 1) >= minExon) {
                std::string info = std::string(c_8);
                size_t p_p_b = info.find("Parent=");
                if (p_p_b == std::string::npos)
                    continue;
                info = info.substr(p_p_b + 7);
                size_t p_p_e = info.find(";");
                info = info.substr(0, p_p_e);

                // process for the type like : chr1	TAIR10	CDS	3760	3913	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
                size_t p_p_b_pre = info.find(",");
                info = info.substr(0, p_p_b_pre);

                if (transcriptHashMap.find(info) != transcriptHashMap.end()) {
                } else {
                    std::string chromosomeName = std::string(name);

                    STRAND strand;
                    if (c_6[0] == '-') {
                        strand = NEGATIVE;
                    } else {
                        strand = POSITIVE;
                    }

                    Transcript transcript1(info, chromosomeName, strand);
                    transcriptHashMap[info] = transcript1;
                }

                GenomeBasicFeature cds(start, end);
                transcriptHashMap[info].addCds(cds);
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

void readGffFile1(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
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
        regex_search(line, match, reg);

        if (match.empty() || line[0] == '#' || line.size() < 9) {
        } else {
            int start = stoi(match[3]);
            int end = stoi(match[4]);
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }
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
        regex_search(line, match, reg);

        if (match.empty() || line[0] == '#' || line.size() < 9) {
        } else {
            int start = stoi(match[3]);
            int end = stoi(match[4]);
            if (start > end) {
                int temp = start;
                start = end;
                end = temp;
            }

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
                transcriptHashMap[information].addCds(cds);
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