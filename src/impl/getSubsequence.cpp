//
// Created by baoxing on 10/10/17.
//

#include "getSubsequence.h"

#include <cerrno>
#include <cstring>

const char replaceCharArray[] = {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};
std::set<char> set_replaceChar(replaceCharArray, replaceCharArray + 10);

char char_replace(char c) {
    if (c == 'A' || c == 'C' || c == 'T' || c == 'G') {
        return c;
    }

    if (c >= 'a' && c <= 'z') {
        c = c - 32;
    }
    if (set_replaceChar.find(c) != set_replaceChar.end()) {
        return 'N';
    }
    if (c == 'U') {
        return 'T';
    }
    return c;
}

namespace {

using FastaEntry = std::tuple<std::string, long, long, int>;
using FastaIndex = std::map<std::string, FastaEntry>;

const FastaEntry *findEntry(const FastaIndex &index,
                            const std::string &sequenceName) {
    const auto found = index.find(sequenceName);
    if (found == index.end()) {
        std::cerr << "could not find " << sequenceName << std::endl;
        return nullptr;
    }
    return &found->second;
}

std::string readInterval(const FastaEntry &entry, int fd,
                         int requestedStart, int requestedEnd) {
    long sequenceSize;
    long offset;
    int lineBases;
    std::tie(std::ignore, sequenceSize, offset, lineBases) = entry;
    if (sequenceSize <= 0 || lineBases <= 0) {
        return std::string();
    }

    int start = requestedStart;
    int end = requestedEnd;
    if (start > end) {
        std::swap(start, end);
    }
    start = std::max(1, start);
    if (start > sequenceSize) {
        start = static_cast<int>(sequenceSize);
        end = start;
    } else if (end > sequenceSize) {
        end = static_cast<int>(sequenceSize);
    }

    const int expectedBases = end - start + 1;
    const int physicalStart = start + (start - 1) / lineBases;
    const int physicalEnd = end + (end - 1) / lineBases;
    // Preserve the legacy one-byte look-ahead. It covers a line terminator at
    // the right edge; the string is then compacted and bounded to the exact
    // requested number of bases.
    const std::size_t bytesToRead =
            static_cast<std::size_t>(physicalEnd - physicalStart + 2);
    std::string sequence(bytesToRead, '\0');
    const ssize_t bytesRead = pread(
            fd, &sequence[0], bytesToRead, offset + physicalStart - 1);
    if (bytesRead < 0) {
        throw std::runtime_error(
                std::string("pread failed: ") + std::strerror(errno));
    }
    sequence.resize(static_cast<std::size_t>(bytesRead));

    const auto compactEnd = std::remove_if(
            sequence.begin(), sequence.end(),
            [](char base) { return base == '\n' || base == '\r'; });
    sequence.erase(compactEnd, sequence.end());
    if (sequence.size() > static_cast<std::size_t>(expectedBases)) {
        sequence.resize(static_cast<std::size_t>(expectedBases));
    }
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                   char_replace);
    return sequence;
}

std::string readIntervalOpeningFile(const FastaEntry &entry,
                                    int start, int end) {
    const std::string &path = std::get<0>(entry);
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error(
                std::string("cannot open FASTA file ") + path + ": " +
                std::strerror(errno));
    }
    try {
        std::string sequence = readInterval(entry, fd, start, end);
        close(fd);
        return sequence;
    } catch (...) {
        close(fd);
        throw;
    }
}

}  // namespace

std::string getSubsequence2(const FastaIndex &index,
                            const std::string &sequenceName,
                            const int &start, const int &end) {
    const FastaEntry *entry = findEntry(index, sequenceName);
    return entry == nullptr
           ? std::string()
           : readIntervalOpeningFile(*entry, start, end);
}

std::string getSubsequence3(const FastaIndex &index, int &fd,
                            const std::string &sequenceName,
                            const int &start, const int &end) {
    const FastaEntry *entry = findEntry(index, sequenceName);
    return entry == nullptr
           ? std::string()
           : readInterval(*entry, fd, start, end);
}

char getCharByPos(const FastaIndex &index,
                  const std::string &sequenceName, const int &position) {
    const std::string base =
            getSubsequence2(index, sequenceName, position, position);
    return base.empty() ? '\0' : base[0];
}

std::string getSubsequence2(const FastaIndex &index,
                            const std::string &sequenceName) {
    const FastaEntry *entry = findEntry(index, sequenceName);
    if (entry == nullptr) {
        return std::string();
    }
    return readIntervalOpeningFile(
            *entry, 1, static_cast<int>(std::get<1>(*entry)));
}

std::string getSubsequence3(const FastaIndex &index, int &fd,
                            const std::string &sequenceName) {
    const FastaEntry *entry = findEntry(index, sequenceName);
    if (entry == nullptr) {
        return std::string();
    }
    return readInterval(*entry, fd, 1,
                        static_cast<int>(std::get<1>(*entry)));
}

std::string getSubsequence3(const FastaIndex &index, int &fd,
                            const std::string &sequenceName,
                            const int &start, const int &end,
                            const STRAND &strand) {
    std::string sequence =
            getSubsequence3(index, fd, sequenceName, start, end);
    if (strand == NEGATIVE) {
        sequence = getReverseComplementary(sequence);
    }
    return sequence;
}

std::string getSubsequence2(const FastaIndex &index,
                            const std::string &sequenceName,
                            const int &start, const int &end,
                            const STRAND &strand) {
    std::string sequence =
            getSubsequence2(index, sequenceName, start, end);
    if (strand == NEGATIVE) {
        sequence = getReverseComplementary(sequence);
    }
    return sequence;
}

size_t getSequenceSizeFromPath2(const FastaEntry &entry) {
    return static_cast<size_t>(std::get<1>(entry));
}

std::string getSubsequence(const std::string &sequence,
                           const int &requestedStart,
                           const int &requestedEnd) {
    int start = requestedStart;
    int end = requestedEnd;
    if (start > end) {
        std::swap(start, end);
    }
    start = std::max(1, start);
    if (start > static_cast<int>(sequence.size())) {
        start = static_cast<int>(sequence.size());
        end = start;
    } else if (end > static_cast<int>(sequence.size())) {
        end = static_cast<int>(sequence.size());
    }
    return sequence.substr(static_cast<std::size_t>(start - 1),
                           static_cast<std::size_t>(end - start + 1));
}

std::string getSubsequence(const std::string &sequence,
                           const int &start, const int &end,
                           const STRAND &strand) {
    std::string result = getSubsequence(sequence, start, end);
    if (strand == NEGATIVE) {
        result = getReverseComplementary(result);
    }
    return result;
}
