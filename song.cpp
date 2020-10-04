/*
 * =====================================================================================
 *
 *       Filename:  song.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/23/2017 21:51:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include "src/controlLayer.h"
#include "./googletest/googletest/include/gtest/gtest.h"


using namespace std;


int main(int argc, char** argv){
//testing::InitGoogleTest(&argc, argv);
//RUN_ALL_TESTS();
//return 0;
    if( argc<=1 ){
        usage();
        return 1;
    }
    std::string program = argv[1];
    if( program.compare("-h") == 0 || program.compare("--help") == 0 ){
        usage();
        exit(1);
    }
    InputParser inputParser (argc, argv);
    string parameterFile;
    std::string exepath = getexepath(argv);
    if( inputParser.cmdOptionExists("-parameter")){
        parameterFile = inputParser.getCmdOption("-parameter");
    }else{
        parameterFile = exepath + "/configure";
    }
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, exepath);
    if( program.compare("gff2seq") == 0 ) {
        return getSequences(--argc, ++argv, parameters);
    } else if ( program.compare("genoAli") == 0 ) {
        return DenoveAssemblyVariantCalling(--argc, ++argv, parameters);
    } else if ( program.compare("proali") == 0 ) {
        return genomeAlignment(--argc, ++argv, parameters);
    } else{
        usage();
    }
    return 0;
}
