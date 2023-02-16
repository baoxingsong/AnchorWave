//
// Created by baoxing on 10/10/17.
//

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <unistd.h>

void readFastaFile(const std::string &filePath, std::map<std::string, std::tuple<std::string, long, long, int> > &map);
