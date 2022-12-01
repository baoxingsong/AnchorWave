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
#include<string.h>
#include<cstring>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "./minimap2/minimap.h"
#include "./minimap2/kseq.h"

int main(int argc, char **argv) {
//
//    testing::InitGoogleTest(&argc, argv);
//    RUN_ALL_TESTS();
//    return 0;

    if (argc <= 1) {
        usage();
        return 1;
    }
    std::string program = argv[1];
    if (program.compare("-h") == 0 || program.compare("--help") == 0) {
        usage();
        exit(1);
    }
    InputParser inputParser(argc, argv);
    std::string parameterFile;
    std::string exepath = getexepath(argv);
    if (inputParser.cmdOptionExists("-parameter")) {
        parameterFile = inputParser.getCmdOption("-parameter");
    } else {
        parameterFile = exepath + "/configure";
    }
    std::map<std::string, std::string> parameters = initialize_paramters(parameterFile, exepath);

#ifdef __AVX512BW__
    std::cerr << "AVX512 is enabled" << std::endl;
#elif __AVX2__
    std::cerr << "AVX2 is enabled" << std::endl;
#elif __SSE4_1__
    std::cerr << "SSE4.1 is enabled" << std::endl;
#elif __SSE2__
    std::cerr << "SSE2 is enabled" << std::endl;
#else
    std::cerr << "The code has not been tested on you hardware platform." << std::endl;
    std::cerr << "If you find anything abnormal, please contact us." << std::endl;
#endif

    if (program.compare("gff2seq") == 0) {
        return gff2seq(--argc, ++argv, parameters);
    } else if (program.compare("proali") == 0) {
        return proportationalAlignment(--argc, ++argv, parameters);
    } else if (program.compare("genoAli") == 0) {
        return genomeAlignment(--argc, ++argv, parameters);
    } else if (program.compare("triAnc") == 0) {
        return tripleAncestral(--argc, ++argv, parameters);
    } else if (program.compare("maf2vcf") == 0) {
        return maf2vcf(--argc, ++argv, parameters);
    } else if (program.compare("sam2maf") == 0) {
        return sam2maf(--argc, ++argv, parameters);
    } else if (program.compare("sam2vcf") == 0) {
        return sam2vcf(--argc, ++argv, parameters);
    } else if (program.compare("evaluateTEAlignment") == 0) {
        return evaluateTEAlignment(--argc, ++argv, parameters);
    } else if (program.compare("sdiToMaf") == 0) {
        return sdiToMaf(--argc, ++argv, parameters);
    } else if (program.compare("ali") == 0) {
        return ali(--argc, ++argv, parameters);
    } else {
        usage();
    }

    return 0;
}
