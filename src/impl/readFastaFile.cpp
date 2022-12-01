//
// Created by baoxing on 10/10/17.
//

#include "readFastaFile.h"

const char replaceCharArray[] = {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};
std::set<char> set_replaceChar(replaceCharArray, replaceCharArray + 10);

char char_replace(char str) {
    if (str == 'A' || str == 'C' || str == 'T' || str == 'G') {
        return str;
    }

    if (str >= 'a' && str <= 'z')
        str = str - 32; // toupper, 32 = 'a' - 'A'

    if (set_replaceChar.find(str) != set_replaceChar.end()) {
        return 'N';
    }

    if (str == 'U') {
        return 'T';
    }

    return str;
}

void readFastaFile(const std::string &filePath, std::map<std::string, std::string> &sequences) {
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit(1);
    }

    std::string name = "";
    std::stringstream sequencestream;
    std::string line = "";
    while (std::getline(infile, line)) {
        if (line[0] == '>') { // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if (name.size() > 0) {
                std::string sequence = sequencestream.str();
                transform(sequence.begin(), sequence.end(), sequence.begin(), char_replace);
                sequences[name] = sequence;
            }

            name = line.substr(1, line.find(" ", 0) - 1);
            if (name[0] == '>') {
                name = name.substr(1, name.find("\t", 0) - 1);
            } else {
                name = name.substr(0, name.find("\t", 0) - 1);
            }
            sequencestream.str(std::string());
        } else {
            sequencestream << line;
        }
    }

    if (name.size() > 0) {
        std::string sequence = sequencestream.str();
        transform(sequence.begin(), sequence.end(), sequence.begin(), char_replace);

        if (sequence.size() > 0) {
            sequences[name] = sequence;
        }
    }

    infile.close();
}

void readFastaFileWorkWithIUPACcode(const std::string &filePath, std::map<std::string, std::string> &sequences) {
    std::ifstream infile(filePath);
    if (!infile.good()) {
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit(1);
    }

    std::string name = "";
    std::stringstream sequencestream;
    std::string line = "";
    while (std::getline(infile, line)) {
        if (line[0] == '>') { // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if (name.size() > 0) {
                std::string sequence = sequencestream.str();
                transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                sequences[name] = sequence;
            }

            name = line.substr(1, line.find(" ", 0) - 1);
            if (name[0] == '>') {
                name = name.substr(1, name.find("\t", 0) - 1);
            } else {
                name = name.substr(0, name.find("\t", 0) - 1);
            }
            sequencestream.str(std::string());
        } else {
            sequencestream << line;
        }
    }

    if (name.size() > 0) {
        std::string sequence = sequencestream.str();
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        if (sequence.size() > 0) {
            sequences[name] = sequence;
        }
    }

    infile.close();
}
