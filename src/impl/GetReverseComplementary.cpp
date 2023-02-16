//
// Created by baoxing on 10/10/17.
//

#include "GetReverseComplementary.h"

std::string getReverseComplementary(const std::string &seq) {
    std::stringstream ss;

    for (int i = seq.length() - 1; i >= 0; i--) {
        char c = seq[i];
        switch (c) {
            case 'A' : c = 'T'; break;
            case 'T' : c = 'A'; break;
            case 'U' : c = 'A'; break;
            case 'C' : c = 'G'; break;
            case 'G' : c = 'C'; break;
            case 'R' : c = 'Y'; break;
            case 'Y' : c = 'R'; break;
            case 'K' : c = 'M'; break;
            case 'M' : c = 'K'; break;
            case 'B' : c = 'V'; break;
            case 'V' : c = 'B'; break;
            case 'D' : c = 'H'; break;
            case 'H':  c = 'D'; break;
        }

        ss << c;
    }

    return ss.str();
}
