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
//std::vector<std::string> split(const std::string &s, char delim);


////input a string, the string could be trasformed to uppercas
//void songToUpCase( std::string& str );
///*  input a string, the transformed would also be returned */
//void songToLowCase( std::string& str);
std::string songStrRemoveBlank(std::string& str);
std::string songStrReplaceAll(std::string& str, const std::string& pattern, const std::string& pattern2);


bool if_file_exists (const std::string& name);
long GetFileSize(const std::string& filename);

bool createdFolder(const std::string & path);
bool ifFileExists (const std::string& list_file);
bool ifAFolder(const std::string & folder);
void getMafftResultListFromFolder( std::string& folder, std::vector<std::string>& fileNameList);
void getListFromTextFile(const std::string& list_file, std::vector<std::string>& list);
void getSetFromTextFile(const std::string& list_file, std::set<std::string>& set);
bool caseInsensitiveStringCompare(const std::string& str1, const std::string& str2);
//for MSA
void getSdiList(const std::string& list_file, std::map<std::string, std::string>& sdiLists);

void overlapGenomeAndGene(std::map <std::string, std::vector<Gene>> & geneMap, std::map <std::string, Fasta> & sequences);

#endif

