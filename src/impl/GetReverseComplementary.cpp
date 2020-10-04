//
// Created by baoxing on 10/10/17.
//

#include "GetReverseComplementary.h"


std::string getReverseComplementary(const std::string& sequence) {
    std::string reserveString;
    reserveString.reserve(sequence.size());
    std::stringstream reversecomplementary(reserveString);
    for (int i = sequence.length() - 1; i >= 0; i--) {
        char c = sequence[i];
        if ('A' == c) {
            c = 'T';
        } else if ('T' == c) {
            c = 'A';
        } else if ('U' == c) {
            c = 'A';
        } else if ('C' == c) {
            c = 'G';
        } else if ('G' == c) {
            c = 'C';
        } else if ('R' == c) {
            c = 'Y';
        } else if ('Y' == c) {
            c = 'R';
        } else if ('K' == c) {
            c = 'M';
        } else if ('M' == c) {
            c = 'K';
        } else if ('B' == c) {
            c = 'V';
        } else if ('V' == c) {
            c = 'B';
        } else if ('D' == c) {
            c = 'H';
        } else if ('H' == c) {
            c = 'D';
        }
        reversecomplementary<< c;
    }
    return reversecomplementary.str();
}