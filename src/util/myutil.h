// =====================================================================================
// 
//       Filename:  myutil.h
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 04:26:31 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//
// =====================================================================================

#ifndef _MYUTIL_H
#define _MYUTIL_H

#include <algorithm>
#include "../model/model.h"
#include <string>
#include <vector>
#include <map>
#include "nucleotideCodeSubstitutionMatrix.h"

void split(const std::string &s, char &delim, std::vector<std::string> &elems);

std::string songStrReplaceAll(std::string &str, const std::string &pattern, const std::string &pattern2);

int32_t min(const int32_t &a, const int32_t &b);

void splitCIGAR(std::string cigarString, std::vector<std::string> &cigarElems);

// trim from start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

#endif
