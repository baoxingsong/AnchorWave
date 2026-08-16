//
// Created by baoxing on 10/10/17.
//

#include "readGffFile.h"
#include "../io/CompressedInput.h"

#include <array>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <utility>

namespace {

std::string readableAnnotationPath(const std::string &filePath) {
    try {
        return anchorwave::io::materializeInputFile(filePath);
    } catch (const std::exception &error) {
        std::cerr << error.what() << std::endl;
        exit(1);
    }
}

bool splitGffColumns(const std::string &line, std::array<std::string, 9> &columns) {
    size_t begin = 0;
    for (size_t column = 0; column < 8; ++column) {
        const size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            return false;
        }
        columns[column] = line.substr(begin, tab - begin);
        begin = tab + 1;
    }
    columns[8] = line.substr(begin);
    return true;
}

bool parseGffCoordinate(const std::string &text, int &coordinate) {
    if (text.empty()) {
        return false;
    }
    char *end = NULL;
    errno = 0;
    const long value = std::strtol(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0' || value < 1 || value > INT_MAX) {
        return false;
    }
    coordinate = static_cast<int>(value);
    return true;
}

std::string gffAttributeValue(const std::string &attributes, const std::string &key) {
    size_t begin = 0;
    while (begin < attributes.size()) {
        size_t end = attributes.find(';', begin);
        if (end == std::string::npos) {
            end = attributes.size();
        }
        size_t fieldBegin = begin;
        while (fieldBegin < end && (attributes[fieldBegin] == ' ' || attributes[fieldBegin] == '\t')) {
            ++fieldBegin;
        }
        if (attributes.compare(fieldBegin, key.size(), key) == 0) {
            size_t valueBegin = fieldBegin + key.size();
            if (valueBegin < end && attributes[valueBegin] == '=') {
                ++valueBegin;
            } else if (valueBegin < end &&
                       (attributes[valueBegin] == ' ' || attributes[valueBegin] == '\t')) {
                while (valueBegin < end &&
                       (attributes[valueBegin] == ' ' || attributes[valueBegin] == '\t')) {
                    ++valueBegin;
                }
            } else {
                begin = end + 1;
                continue;
            }
            if (valueBegin < end && attributes[valueBegin] == '"') {
                ++valueBegin;
                if (end > valueBegin && attributes[end - 1] == '"') {
                    --end;
                }
            }
            return attributes.substr(valueBegin, end - valueBegin);
        }
        begin = end + 1;
    }
    return std::string();
}

void readParentedFeatures(const std::string &filePath,
                          std::map<std::string, std::vector<Transcript> > &transcriptHashSet,
                          const std::string &featureType,
                          const int minExon) {
    std::map<std::string, Transcript> transcriptHashMap;
    std::ifstream infile(readableAnnotationPath(filePath));
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::array<std::string, 9> columns;
        if (!splitGffColumns(line, columns) || columns[2] != featureType) {
            continue;
        }

        int start = 0;
        int end = 0;
        if (!parseGffCoordinate(columns[3], start) || !parseGffCoordinate(columns[4], end)) {
            continue;
        }
        if (start > end) {
            std::swap(start, end);
        }
        if ((end - start + 1) < minExon) {
            continue;
        }

        std::string parent = gffAttributeValue(columns[8], "Parent");
        if (parent.empty()) {
            // Accept conventional GTF as advertised by the command-line interface.
            parent = gffAttributeValue(columns[8], "transcript_id");
        }
        const size_t comma = parent.find(',');
        if (comma != std::string::npos) {
            parent.erase(comma);
        }
        if (parent.empty()) {
            continue;
        }

        if (transcriptHashMap.find(parent) == transcriptHashMap.end()) {
            const STRAND strand = columns[6] == "-" ? NEGATIVE : POSITIVE;
            transcriptHashMap[parent] = Transcript(parent, columns[0], strand);
        }
        transcriptHashMap[parent].addCds(GenomeBasicFeature(start, end));
    }

    for (std::map<std::string, Transcript>::iterator it = transcriptHashMap.begin();
         it != transcriptHashMap.end(); ++it) {
        it->second.updateInforCds();
        const std::string chromosome = it->second.getChromeSomeName();
        transcriptHashSet[chromosome].push_back(std::move(it->second));
    }
    for (std::map<std::string, std::vector<Transcript> >::iterator it = transcriptHashSet.begin();
         it != transcriptHashSet.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](const Transcript &a, const Transcript &b) {
            return a.getPStart() < b.getPStart();
        });
    }
}

}  // namespace

void get_map_from_gff(const std::string &filePath, std::map <std::string, std::string> &map_transcript_to_gene) {
    std::ifstream infile(readableAnnotationPath(filePath));
    if (!infile.good()) {
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::array<std::string, 9> columns;
        if (!splitGffColumns(line, columns)) {
            continue;
        }

        const std::string id = gffAttributeValue(columns[8], "ID");
        std::string parent = gffAttributeValue(columns[8], "Parent");
        const size_t comma = parent.find(',');
        if (comma != std::string::npos) {
            parent.erase(comma);
        }
        if (!id.empty() && !parent.empty()) {
            map_transcript_to_gene[id] = parent;
            continue;
        }

        const std::string geneId = gffAttributeValue(columns[8], "geneID");
        if (!id.empty() && !geneId.empty()) {
            map_transcript_to_gene[id] = geneId;
            continue;
        }

        const std::string transcriptId = gffAttributeValue(columns[8], "transcript_id");
        const std::string gtfGeneId = gffAttributeValue(columns[8], "gene_id");
        if (!transcriptId.empty() && !gtfGeneId.empty()) {
            map_transcript_to_gene[transcriptId] = gtfGeneId;
        }
    }

    infile.close();
}

void readGffFile(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &type, const int &minExon) {
    std::map <std::string, Transcript> transcriptHashMap;
    std::ifstream infile(readableAnnotationPath(filePath));
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
        const std::string chromosome = it->second.getChromeSomeName();
        transcriptHashSet[chromosome].push_back(std::move(it->second));
    }

    for (std::map < std::string, std::vector < Transcript >> ::iterator it = transcriptHashSet.begin(); it != transcriptHashSet.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](const Transcript &a, const Transcript &b) {
            return a.getPStart() < b.getPStart();
        });
    }
}

void readGffFile1(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
    (void)cdsParentRegex;
    readParentedFeatures(filePath, transcriptHashSet, "CDS", minExon);
}
void readGffFile_exon(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon) {
    (void)cdsParentRegex;
    readParentedFeatures(filePath, transcriptHashSet, "exon", minExon);
}
