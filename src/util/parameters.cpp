/*
 * =====================================================================================
 *
 *       Filename:  parameters.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/03/2017 23:33:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#include "parameters.h"
#include "myutil.h"
#include <fstream>
#include <regex>
#include <stdio.h>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>


std::map<std::string, std::string> initialize_paramters(std::string parameters_file, std::string exepath){
    std::map<std::string, std::string> parameters;
    parameters["exepath"] = exepath;
    std::ifstream infile(parameters_file);
    if( ! infile.good()){
        std::cerr << "error in opening parameter file " << parameters_file << std::endl;
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
            parameters[key] = value;
        }
    }
    return parameters;
}

std::string get_parameters(std::string parameter, std::map<std::string, std::string>& parameters){
    if( parameters.find(parameter) == parameters.end() ){
        std::cerr << "could not find parameter: " << parameter << " in parameters file" << std::endl;
        exit(1);
    }
    return parameters[parameter];
}

std::string createdTempFolder(std::map<std::string, std::string>& parameters){
    std::string tempFolder = get_parameters("tempFolder", parameters);
    struct stat info;
    bool ifTempFolderCreated = false;
    while( !ifTempFolderCreated ){
        if( stat( &tempFolder[0], &info ) != 0 ){
            std::string command = "mkdir " + tempFolder;
            system(&command[0]);
        } else if( info.st_mode & S_IFDIR ){
            ifTempFolderCreated=true;
        } else {
            std::cout << "could not access temp folder: " +  tempFolder << std::endl;
            exit(2);
        }
    }
    return tempFolder;
}

std::string createdMsaPreFolder(std::map<std::string, std::string>& parameters){
    std::string tempFolder = get_parameters("MsaPreFolder", parameters);
    struct stat info;
    bool ifTempFolderCreated = false;
    while( !ifTempFolderCreated ){
        if( stat( &tempFolder[0], &info ) != 0 ){
            std::string command = "mkdir " + tempFolder;
            system(&command[0]);
        } else if( info.st_mode & S_IFDIR ){
            ifTempFolderCreated=true;
        } else {
            std::cout << "could not access MsaPreFolder folder: " +  tempFolder << std::endl;
            exit(2);
        }
    }
    return tempFolder;
}

std::string generateUUID(){
    return sole::uuid0().str();
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    if (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result = buffer.data();
    }
    return result;
}

std::string generateUUID(std::string& prefix){
    return prefix + "." + sole::uuid0().str();
}

std::string getexepath(char** argv) {
    std::string exeFile = std::string(argv[0]);
    std::regex reg("\\/");
    std::smatch match;
    regex_search(exeFile, match, reg);
    if( match.empty() ){
        std::string command = "which " + exeFile;
        exeFile = exec(command.c_str());
        std::regex e ("\\s+$");
        std::string temp="";
        std::regex_replace(std::back_inserter(temp), exeFile.begin(), exeFile.end(), e, "");
        exeFile = temp;
        exeFile.pop_back();
    }
    std::regex e ("\\/[\\w._ ]+$");
    std::string result="";
    std::regex_replace(std::back_inserter(result), exeFile.begin(), exeFile.end(), e, "");
    return result;
}

std::string getexepath(std::map<std::string, std::string>& parameters){
    return parameters["exepath"];
}
