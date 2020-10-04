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

//
//std::vector<std::string> split(const std::string &s, char delim) {
//    std::vector<std::string> elems;
//    split(s, delim, elems);
//    return elems;
//}


//void songToUpCase( std::string& str  ){
//    transform(str.begin(), str.end(), str.begin(),::toupper);
//}
////*  input a string, the transformed would also be returned */
//void songToLowCase( std::string& str   ){
//    transform(str.begin(), str.end(), str.begin(), ::tolower);
//}
/* this implemention is very slow, maybe should find a better way to redo it */
std::string songStrRemoveBlank( std::string& str  ){
    std::string pattern="\\s";
    std::string pattern2="";
    return songStrReplaceAll( str, pattern, pattern2 );
}

std::string songStrReplaceAll( std::string& str, const std::string& pattern, const std::string& pattern2  ){
    std::regex vowel_re(pattern);
    str=std::regex_replace(str, vowel_re, pattern2);
    return str;
}




// orf checking related function end

//translation related begin
// this function is private here


//std::string nA2AA(std::string& seq){
//    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
//    std::string pattern="-";
//    std::string pattern2="";
//    seq=songStrReplaceAll(seq, pattern, pattern2);
//    std::stringstream ssSb;
//    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
//    return ssSb.str();
//}

//std::string nA2AANoDeletedIndel(std::string& seq){
//    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
//    std::stringstream ssSb;
//    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
//    return ssSb.str();
//}



//std::string nA2AANoDeletedIndel(std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
//    std::stringstream ssSb;
//    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
//    return ssSb.str();
//}
//translation related end

bool if_file_exists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

long GetFileSize( const std::string& filename) {
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


bool createdFolder(const std::string & path){
    struct stat info;
    if( stat( &path[0], &info ) != 0 ){
        std::string command = "mkdir " + path;
        system(&command[0]);
        return true;
    }else if( info.st_mode & S_IFDIR ){
        return true;
    } else {
        std::cout << "could not access temp folder: " +  path << std::endl;
        exit(2);
    }
    return false;
}

bool ifFileExists (const std::string& list_file){
    std::ifstream infile(list_file);
    return infile.good();
}

bool ifAFolder(const std::string & folder){
    struct stat info;
    const char * path = folder.c_str();
    if( stat( path, &info ) != 0 ){
        std::cerr << folder << " does not exists" <<std::endl;
        exit(1);
    }else if( info.st_mode & S_IFDIR ){
        return true;
    } else {
        std::cerr << folder << " is not a directory" <<std::endl;
        exit(2);
    }
    return false;
}

void getMafftResultListFromFolder( std::string& folder, std::vector<std::string>& fileNameList){
    int len = folder.length();
    if (folder[len-1] != '/'){
        folder = folder + "/";
    }
    if( ifAFolder(folder)){
        const char * path = folder.c_str();
        DIR * dir = opendir(path);
        struct dirent * file;
        while ((file = readdir(dir)) != NULL) {
            std::string fileName = file->d_name;
            std::string pattern = "mafft";
            if( fileName.find(pattern) != std::string::npos ){
                fileNameList.push_back(std::string(fileName));
//                std::cout << folder << " listing regular file " << file->d_name << std::endl;
            }
        }
    }
}

void getListFromTextFile(const std::string& list_file, std::vector<std::string>& list){
    std::ifstream infile(list_file);
    if( ! infile.good()){
        std::cerr << "error in opening vector list file " << list_file << std::endl;
        exit (1);
    }
    std::string line="";
    std::regex reg("^(\\S+)\\s*");
    while (std::getline(infile, line)) {
        if (line.compare(0, 1, "#") == 0) {
            continue;
        }
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty() ) {
        } else {
            std::string key = match[1];
            list.push_back(key);
        }
    }
}
void getSetFromTextFile(const std::string& list_file, std::set<std::string>& set){
    std::ifstream infile(list_file);
    if( ! infile.good()){
        std::cerr << "error in opening set list file " << list_file << std::endl;
        exit (1);
    }
    std::string line="";
    std::regex reg("^(\\S+)\\s*");
    while (std::getline(infile, line)) {
        if (line.compare(0, 1, "#") == 0) {
            continue;
        }
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty() ) {
        } else {
            std::string key = match[1];
            set.insert(key);
        }
    }
}


bool caseInsensitiveStringCompare(const std::string& str1, const std::string& str2) {
    if (str1.size() != str2.size()) {
        return false;
    }
    for (std::string::const_iterator c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2) {
        if (tolower(*c1) != tolower(*c2)) {
            return false;
        }
    }
    return true;
}

//for MSA
void getSdiList(const std::string& list_file, std::map<std::string, std::string>& sdiLists ){
    std::ifstream infile(list_file);
    if( ! infile.good()){
        std::cerr << "error in opening variants list file " << list_file << std::endl;
        exit (1);
    }
    std::string line="";
    std::regex reg("^(\\S+)\\s+([\\s\\S]+)");
    while (std::getline(infile, line)) {
        if (line.compare(0, 1, "#") == 0) {
            continue;
        }
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty() ) {
        } else {
            std::string key = match[1];
            std::string value = match[2];
            sdiLists[key] = value;
        }
    }
}



void overlapGenomeAndGene(std::map <std::string, std::vector<Gene>> & geneMap, std::map <std::string, Fasta> & sequences){

    // remove those chrs not shared by the fasta file and gff file begin
    std::set <std::string> refChrToRemove;
    for( std::map<std::string, std::vector<Gene> >::iterator it0=geneMap.begin(); it0!=geneMap.end(); ++it0 ){
        if ( sequences.find(it0->first) == sequences.end() ){
            refChrToRemove.insert(it0->first);
        }
    }
    for( std::map<std::string, Fasta >::iterator it0=sequences.begin(); it0!=sequences.end(); ++it0 ){
        if ( geneMap.find(it0->first) == geneMap.end() ){
            refChrToRemove.insert(it0->first);
        }
    }
    for( std::string chr : refChrToRemove ){
        if( geneMap.find(chr) != geneMap.end() ){
            geneMap.erase(chr);
        }
        if( sequences.find(chr) != sequences.end() ){
            sequences.erase(chr);
        }
    }
    // remove those chrs not shared by the fasta file and gff file end
}