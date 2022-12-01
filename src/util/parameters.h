/*
 * =====================================================================================
 *
 *       Filename:  parameters.h
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

#ifndef VARIANTSANNOTATIONSOFTWARE_PARAMETERS_H
#define VARIANTSANNOTATIONSOFTWARE_PARAMETERS_H

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <map>

std::map<std::string, std::string> initialize_paramters(std::string parameters_file, std::string exepath);

std::string get_parameters(std::string parameter, std::map<std::string, std::string> &parameters);

std::string getexepath(char **argv);

#endif