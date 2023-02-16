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
#pragma once

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

void split(const std::string &s, char &delim, std::vector<std::string> &elems);

void splitCIGAR(std::string cigarString, std::vector<std::string> &cigarElems);
