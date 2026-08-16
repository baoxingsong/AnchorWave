//
// Created by baoxing on 10/10/17.
//

#include "readFastaFile.h"
#include "../io/CompressedInput.h"

void readFastaFile(const std::string &filePath, std::map<std::string, std::tuple<std::string, long, long, int> > &map) {
    std::string readablePath;
    try {
        readablePath = anchorwave::io::materializeInputFile(filePath);
    } catch (const std::exception &error) {
        std::cerr << error.what() << std::endl;
        exit(1);
    }

    std::ifstream infile(readablePath);
    if (!infile.good()) {
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit(1);
    }

    std::string name = "";
    std::string line;

    long offset = 0;
    size_t size = 0;
    size_t line_bases = 0;
    size_t line_last = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) {
            ++offset;
            continue;
        }
        if (line[0] == '>') {
            if (name.size() > 0) {
                map[name] = std::make_tuple(readablePath, size, offset, line_bases);

                if (line_bases > 0) {
                    if(line_last == line_bases) {
                        offset += size + size/line_bases;
                    }
                    else {
                        offset += size + size/line_bases + 1;
                    }
                }
            }

            size_t p_sp = line.find_first_of(" \t");
            if(p_sp != std::string::npos) {
                name = line.substr(1, p_sp - 1);
            }
            else {
                name = line.substr(1);
            }

            size = 0;
            line_bases = 0;
            offset += line.size() + 1;
        } else {
            size += line.size();
            line_last = line.size();

            if(line.size() > line_bases) {
                line_bases = line.size();
            }
        }
    }

    if (name.size() > 0) {
        map[name] = std::make_tuple(readablePath, size, offset, line_bases);
    }

    infile.close();
}
