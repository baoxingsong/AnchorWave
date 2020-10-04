/*
 * =====================================================================================
 *
 *       Filename:  FileWriteCatch.cc
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




 ************************************************************************/


#include "FileWriteCatch.h"

void outputContent(std::string const & filePath, std::stringstream const & content, std::atomic_int & number_of_runing_threads){
    std::ofstream ofile;
    ofile.open( filePath,  std::ofstream::app | std::ofstream::out);
    ofile << content.rdbuf();
    ofile.close();
    number_of_runing_threads--;
}


void FileWriteCatch::writeOut( int const & maxThread){

    std::atomic_int number_of_runing_threads(0);

    for( std::map<std::string, std::stringstream>::iterator fileName_content_mapi=fileName_content_map.begin();
         fileName_content_mapi!=fileName_content_map.end(); ++fileName_content_mapi ){

        bool isThisThreadUnrun = true;
        while (isThisThreadUnrun) {
            if (number_of_runing_threads < (maxThread))  {
                std::thread t(outputContent, std::ref(fileName_content_mapi->first), std::ref(fileName_content_mapi->second), std::ref(number_of_runing_threads));
                ++number_of_runing_threads;
                t.detach();
                isThisThreadUnrun = false;
                break;
            } else {
                usleep(100);
            }
        }
    }
    while (number_of_runing_threads > 0) {// wait for all the thread
        usleep(100);
    }
    for( std::map<std::string, std::stringstream>::iterator fileName_content_mapi=fileName_content_map.begin();
        fileName_content_mapi!=fileName_content_map.end(); ++fileName_content_mapi ){
//
//        std::ofstream ofile;
//        ofile.open( fileName_content_mapi->first,  std::ofstream::app | std::ofstream::out);
//        ofile << fileName_content_mapi->second.rdbuf() << std::endl;
//        ofile.close();

        fileName_content_mapi->second.str(std::string());
    }

}
void FileWriteCatch::addContent(const std::string& fileName, const std::string & content){
    if (this->fileName_content_map.find(fileName) != this->fileName_content_map.end() ){

    }else{
        fileName_content_map[fileName] = std::stringstream();
    }
    fileName_content_map[fileName] << content;
}

