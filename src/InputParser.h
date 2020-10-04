/*
 * =====================================================================================
 *
 *       Filename:  InputParser.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 01:10:23
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#ifndef _INPUTPARSER_H
#define _INPUTPARSER_H
#include <string>
#include <vector>
#include <algorithm>
class InputParser{
    private:
        std::vector <std::string> tokens;
    public:
        InputParser (int &argc, char **argv);
        const std::string getCmdOption( std::string &option);
        std::string getCmdOption( const char* o);
        bool cmdOptionExists( std::string &option);
        bool cmdOptionExists( const char* o);
};

void usage( );
#endif
