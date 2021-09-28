// =====================================================================================
// 
//       Filename:  myutil.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 04:53:23 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//
// =====================================================================================
/*************************************************************************


************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include "myutil.h"
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>

#include <stdio.h>
#include <errno.h>
#include <dirent.h>

#ifdef __unix
#include <linux/limits.h>
#endif

void split(const std::string &s, char& delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (item.length() > 0) {
            elems.push_back(item);
        }
    }
}


std::string songStrReplaceAll( std::string& str, const std::string& pattern, const std::string& pattern2  ){
    std::regex vowel_re(pattern);
    str=std::regex_replace(str, vowel_re, pattern2);
    return str;
}


int32_t min (  const int32_t & a, const int32_t & b ){
    if ( a < b ){
        return 1;
    }
    return b;
}


void splitCIGAR( std::string & cigarString, std::vector<std::string> & cigarElems) {
    std::regex reg("([0-9]+[MIDNSHPX=])");
    std::smatch match;
    while( regex_search(cigarString, match, reg) ){
        cigarElems.push_back( match[1] );
        cigarString = match.suffix().str();
    }
}

