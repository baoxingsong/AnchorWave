/*
 * =====================================================================================
 *
 *       Filename:  FileWriteCatch.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  14/11/2017 15:47:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *        Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************

The class if designed to append content to huge mount of file


 ************************************************************************/

#ifndef ANNOTATIONLIFTOVER_FILEWRITECATCH_H
#define ANNOTATIONLIFTOVER_FILEWRITECATCH_H

#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <thread>
#include <unistd.h>
#include <atomic>
#include <mutex>

class FileWriteCatch {
    private:
        std::map<std::string, std::stringstream> fileName_content_map;
    public:
        void writeOut(int const & maxThread);
        void addContent(const std::string& fileName, const std::string & content);
};

#endif //ANNOTATIONLIFTOVER_FILEWRITECATCH_H
