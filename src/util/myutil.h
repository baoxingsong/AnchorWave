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
/*************************************************************************




 ************************************************************************/

#ifndef _MYUTIL_H
#define _MYUTIL_H
#include <algorithm>
#include "../model/model.h"
#include <string>
#include <vector>
#include <map>
#include "nucleotideCodeSubstitutionMatrix.h"

void split(const std::string &s, char& delim,std::vector<std::string> &elems);
std::string songStrReplaceAll(std::string& str, const std::string& pattern, const std::string& pattern2);

#endif
