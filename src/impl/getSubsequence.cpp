//
// Created by baoxing on 10/10/17.
//

#include "getSubsequence.h"


const char replaceCharArray[] = {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};
std::set<char> set_replaceChar(replaceCharArray, replaceCharArray + 10);

char char_replace(char c) {
    if (c == 'A' || c == 'C' || c == 'T' || c == 'G') {
        return c;
    }

    if (c >= 'a' && c <= 'z')
        c = c - 32; // toupper, 32 = 'a' - 'A'

    if (set_replaceChar.find(c) != set_replaceChar.end()) {
        return 'N';
    }

    if (c == 'U') {
        return 'T';
    }

    return c;
}

std::string getSubsequence2(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName, const int &_start, const int &_end) {
    if (map.find(seqName) != map.end()) {
        std::string path;
        long size_all;
        long offset;
        int line_bases;

        std::tie(path, size_all, offset, line_bases) = map[seqName];
        int fd = open(path.c_str(), O_RDONLY);

        int start = _start;
        int end = _end;

        if (start > end) {
            int temp = start;
            start = end;
            end = temp;
        }

        if (start < 1) {
            start = 1;
        }

        if (start > size_all) {
            start = size_all;
            end = start;
        } else if (end > size_all) {
            end = size_all;
        }

        int size = end - start;

        // 带换行的start和end
        int count_n_start = 0;
        if(start <= line_bases) {
            count_n_start = 0;
        }
        else {
            count_n_start = (start - 1) / line_bases;
        }

        int start2 = start + count_n_start;
        int count_n_end = 0;
        if(end <= line_bases) {
            count_n_end = 0;
        }
        else {
            count_n_end = (end - 1) / line_bases;
        }

        int end2 = end + count_n_end;
        int size_real = end2 - start2 + 1;
        char *buf = new char[size_real + 2];

        ssize_t ret = pread(fd, buf, size_real + 1, offset + start2 - 1);
        if(ret == -1) {
            std::cerr << "pread error! " << std::endl;
        }

        buf[size_real + 1] = '\0';
        std::string ret_str = std::string(buf);
        ret_str.erase(std::remove(ret_str.begin(), ret_str.end(), '\n'), ret_str.end());

        size_t s2 = ret_str.size();
        if(s2 - size > 1) {
            ret_str = ret_str.substr(0, size + 1);
        }

        transform(ret_str.begin(), ret_str.end(), ret_str.begin(), char_replace);

        close(fd);
        delete [] buf;

        return ret_str;
    }
    else {
        std::cerr << "could not find " << seqName << std::endl;
    }

    return "";
}

std::string getSubsequence3(std::map<std::string, std::tuple<std::string, long, long, int> > &map, int &fd, const std::string &seqName, const int &_start, const int &_end) {
    if (map.find(seqName) != map.end()) {
        std::string path;
        long size_all;
        long offset;
        int line_bases;

        std::tie(path, size_all, offset, line_bases) = map[seqName];
        int start = _start;
        int end = _end;

        if (start > end) {
            int temp = start;
            start = end;
            end = temp;
        }

        if (start < 1) {
            start = 1;
        }

        if (start > size_all) {
            start = size_all;
            end = start;
        } else if (end > size_all) {
            end = size_all;
        }

        int size = end - start;

        // 带换行的start和end
        int count_n_start = 0;
        if(start <= line_bases) {
            count_n_start = 0;
        }
        else {
            count_n_start = (start - 1) / line_bases;
        }

        int start2 = start + count_n_start;

        int count_n_end = 0;
        if(end <= line_bases) {
            count_n_end = 0;
        }
        else {
            count_n_end = (end - 1) / line_bases;
        }

        int end2 = end + count_n_end;
        int size_real = end2 - start2 + 1;

        char *buf = new char[size_real + 2];

        ssize_t ret = pread(fd, buf, size_real + 1, offset + start2 - 1);
        if(ret == -1) {
            std::cerr << "pread error! " << std::endl;
        }

        buf[size_real + 1] = '\0';

        std::string ret_str = std::string(buf);
        ret_str.erase(std::remove(ret_str.begin(), ret_str.end(), '\n'), ret_str.end());
        size_t s2 = ret_str.size();
        if(s2 - size > 1) {
            ret_str = ret_str.substr(0, size + 1);
        }

        transform(ret_str.begin(), ret_str.end(), ret_str.begin(), char_replace);
        delete [] buf;

        return ret_str;
    }
    else {
        std::cerr << "could not find " << seqName << std::endl;
    }

    return "";
}

char getCharByPos(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName, const int &_pos) {
    if (map.find(seqName) != map.end()) {
        std::string path;
        long size_all;
        long offset;
        int line_bases;

        std::tie(path, size_all, offset, line_bases) = map[seqName];
        int fd = open(path.c_str(), O_RDONLY);

        int start = _pos;

        if (start < 1) {
            start = 1;
        }

        if (start > size_all) {
            start = size_all;
        }

        // 带换行的start和end
        int count_n_start = 0;
        if(start <= line_bases) {
            count_n_start = 0;
        }
        else {
            count_n_start = (start - 1) / line_bases;
        }

        int start2 = start + count_n_start;

        int size_real = 1;
        char *buf = new char[size_real + 2];

        ssize_t ret = pread(fd, buf, size_real + 1, offset + start2 - 1);
        if(ret == -1) {
            std::cerr << "pread error! " << std::endl;
        }

        buf[size_real + 1] = '\0';

        std::string ret_str = std::string(buf);

        char ret_c = char_replace(ret_str[0]);

        close(fd);
        delete [] buf;

        return ret_c;
    }
    else {
        std::cerr << "could not find " << seqName << std::endl;
    }

    return '\0';
}

std::string getSubsequence2(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName) {
    if (map.find(seqName) != map.end()) {
        std::string path;
        long size_all;
        long offset;
        int line_bases;

        std::tie(path, size_all, offset, line_bases) = map[seqName];
        int fd = open(path.c_str(), O_RDONLY);

        int start = 0;
        int end = size_all;

        int size = end - start;

        // 带换行的start和end
        int count_n_start = 0;
        if(start <= line_bases) {
            count_n_start = 0;
        }
        else {
            count_n_start = (start - 1) / line_bases;
        }

        int start2 = start + count_n_start;

        int count_n_end = 0;
        if(end <= line_bases) {
            count_n_end = 0;
        }
        else {
            count_n_end = (end - 1) / line_bases;
        }

        int end2 = end + count_n_end;
        int size_real = end2 - start2 + 1;
        char *buf = new char[size_real + 2];

        ssize_t ret = pread(fd, buf, size_real + 1, offset + start2 - 1);
        if(ret == -1) {
            std::cerr << "pread error! " << std::endl;
        }

        buf[size_real + 1] = '\0';

        std::string ret_str = std::string(buf);
        ret_str.erase(std::remove(ret_str.begin(), ret_str.end(), '\n'), ret_str.end());

        size_t s2 = ret_str.size();
        if(s2 - size > 1) {
            ret_str = ret_str.substr(0, size + 1);
        }

        transform(ret_str.begin(), ret_str.end(), ret_str.begin(), char_replace);

        close(fd);
        delete [] buf;

        return ret_str;
    }
    else {
        std::cerr << "could not find " << seqName << std::endl;
    }

    return "";
}

std::string getSubsequence3(std::map<std::string, std::tuple<std::string, long, long, int> > &map, int &fd, const std::string &seqName) {
    if (map.find(seqName) != map.end()) {
        std::string path;
        long size_all;
        long offset;
        int line_bases;

        std::tie(path, size_all, offset, line_bases) = map[seqName];

        int start = 0;
        int end = size_all;

        int size = end - start;

        // 带换行的start和end
        int count_n_start = 0;
        if(start <= line_bases) {
            count_n_start = 0;
        }
        else {
            count_n_start = (start - 1) / line_bases;
        }

        int start2 = start + count_n_start;

        int count_n_end = 0;
        if(end <= line_bases) {
            count_n_end = 0;
        }
        else {
            count_n_end = (end - 1) / line_bases;
        }

        int end2 = end + count_n_end;
        int size_real = end2 - start2 + 1;
        char *buf = new char[size_real + 2];

        ssize_t ret = pread(fd, buf, size_real + 1, offset + start2 - 1);
        if(ret == -1) {
            std::cerr << "pread error! " << std::endl;
        }

        buf[size_real + 1] = '\0';
        std::string ret_str = std::string(buf);
        ret_str.erase(std::remove(ret_str.begin(), ret_str.end(), '\n'), ret_str.end());

        size_t s2 = ret_str.size();
        if(s2 - size > 1) {
            ret_str = ret_str.substr(0, size + 1);
        }

        transform(ret_str.begin(), ret_str.end(), ret_str.begin(), char_replace);
        delete [] buf;

        return ret_str;
    }
    else {
        std::cerr << "could not find " << seqName << std::endl;
    }

    return "";
}

std::string getSubsequence3(std::map<std::string, std::tuple<std::string, long, long, int> > &map, int &fd, const std::string &seqName, const int &_start, const int &_end, const STRAND &strand) {
//    std::cout << "getSubsequence2, seqName " << seqName << ", _start " << _start << ", _end " << _end << ", strand " << strand << std::endl;
    if (map.find(seqName) != map.end()) {
        if (strand == POSITIVE) {
            return getSubsequence3(map, fd, seqName, _start, _end);
        } else {
            std::string seq = getSubsequence3(map, fd, seqName, _start, _end);
            return getReverseComplementary(seq);
        }
    }

    return "";
}

std::string getSubsequence2(std::map<std::string, std::tuple<std::string, long, long, int> > &map, const std::string &seqName, const int &_start, const int &_end, const STRAND &strand) {
//    std::cout << "getSubsequence2, seqName " << seqName << ", _start " << _start << ", _end " << _end << ", strand " << strand << std::endl;
    if (map.find(seqName) != map.end()) {
        if (strand == POSITIVE) {
            return getSubsequence2(map, seqName, _start, _end);
        } else {
            std::string seq = getSubsequence2(map, seqName, _start, _end);
            return getReverseComplementary(seq);
        }
    }

    return "";
}

size_t getSequenceSizeFromPath2(std::tuple<std::string, long, long, int> &t) {
    long size_all;
    std::tie(std::ignore, size_all, std::ignore, std::ignore) = t;

    return size_all;
}

std::string getSubsequence(const std::string &sequence, const int &_start, const int &_end) {
    int start = _start;
    int end = _end;
    if (start > end) {
        int temp = start;
        start = end;
        end = temp;
    }
    if (start < 1) {
        start = 1;
    }
    if (start > sequence.size()) {
        start = sequence.size();
        end = start;
    } else if (end > sequence.size()) {
        end = sequence.size();
    }
    return sequence.substr(start - 1, end - start + 1);
}

std::string getSubsequence(const std::string &sequence, const int &_start, const int &_end, const STRAND &strand) {
    if (strand == POSITIVE) {
        return getSubsequence(sequence, _start, _end);
    } else {
        std::string seq = getSubsequence(sequence, _start, _end);
        return getReverseComplementary(seq);
    }
}
