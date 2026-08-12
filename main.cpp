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
#include "src/cpu_arch.h"
#include "src/controlLayer.h"
#include "src/trio/TrioGraphCommand.h"

int main(int argc, char **argv) {
    if (argc <= 1) {
        usage();
        return 1;
    }
    std::string program = argv[1];
    if (program.compare("-h") == 0 || program.compare("--help") == 0) {
        usage();
        return 0;
    }
    if (program.compare("-v") == 0 || program.compare("--version") == 0) {
        std::cout << PROGRAMNAME << " v" << VERSION << std::endl;
        return 0;
    }

#ifdef __SSE2NEON__
    std::cerr << "ARM NEON is enabled through SSE2NEON" << std::endl;
#elif defined(__AVX512BW__)
    std::cerr << "AVX512 is enabled" << std::endl;
#elif defined(__AVX2__)
    std::cerr << "AVX2 is enabled" << std::endl;
#elif defined(__SSE4_1__)
    std::cerr << "SSE4.1 is enabled" << std::endl;
#elif defined(__SSE2__)
    std::cerr << "SSE2 is enabled" << std::endl;
#else
    std::cerr << "The code has not been tested on you hardware platform." << std::endl;
    std::cerr << "If you find anything abnormal, please contact us." << std::endl;
#endif

    if (program.compare("gff2seq") == 0) {
        return gff2seq(--argc, ++argv);
    }
    else if (program.compare("proali") == 0) {
        return proportionalAlignment(--argc, ++argv);
    }
    else if (program.compare("genoAli") == 0) {
        return genomeAlignment(--argc, ++argv);
    }
    else if (program.compare("ali") == 0) {
        return ali(--argc, ++argv);
    }
    else if (program.compare("geno") == 0) {
        return geno(--argc, ++argv);
    }
    else if (program.compare("pro") == 0) {
        return pro(--argc, ++argv);
    }
    else if (program.compare("trioGraphAli") == 0) {
        return trioGraphAli(--argc, ++argv);
    }

    usage();

    return 0;
}
