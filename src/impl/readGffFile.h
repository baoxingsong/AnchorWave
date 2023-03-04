//
// Created by baoxing on 10/10/17.
//

#pragma once

#include "../model/Transcript.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <regex>

void get_map_from_gff(const std::string &filePath, std::map<std::string, std::string> &transcript_to_gene_map);

void readGffFile(const std::string &filePath, std::map<std::string, std::vector<Transcript> > &transcriptHashSet, const std::string &type, const int &minExon);

void readGffFile1(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon);

void readGffFile_exon(const std::string &filePath, std::map <std::string, std::vector<Transcript>> &transcriptHashSet, const std::string &cdsParentRegex, const int &minExon);