cmake_minimum_required(VERSION 3.0)
project(proali)

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 14)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
#else()
#    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()


SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib")
set(CMAKE_INSTALL_RPATH "/usr/local/lib/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_FLAGS "-Wall -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")
set(CMAKE_C_FLAGS "-Wall -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread")


#for google unit test
add_subdirectory(./googletest/googletest)
include_directories(googletest/googletest/include googletest/googletest)


add_executable(proali
        googletest
        googletest/googletest
        #       ksw2_extz2_sse.c
#       ksw2.h
#       src/myImportandFunction/alignNeedlemanForTranscript_simd.cpp
        src/model/GenomeBasicFeature.cpp
        src/model/GenomeBasicFeature.h
        src/controlLayer.cpp
        src/controlLayer.h
        src/InputParser.cpp
        src/InputParser.h
        src/model/model.h
        src/util/myutil.cpp
        src/util/myutil.h
        src/util/nucleotideCodeSubstitutionMatrix.cpp
        src/util/nucleotideCodeSubstitutionMatrix.h
        src/util/parameters.cpp
        src/util/parameters.h
        src/sole.h
        song.cpp
        src/model/Fasta.cpp
        src/model/Fasta.h
        src/model/variant.cpp
        src/model/variant.h
        src/model/Transcript.cpp
        src/model/Transcript.h
        src/model/Gene.cpp
        src/model/Gene.h
        src/model/Matrix.cpp
        src/model/Matrix.h
        src/model/MsaWindow.cpp
        src/model/MsaWindow.h
        src/model/MsaSingleRecord.cpp
        src/model/MsaSingleRecord.h
        src/model/MsaFileRecord.cpp
        src/model/MsaFileRecord.h
        src/model/TwoSeqOfMsaResult.cpp
        src/model/TwoSeqOfMsaResult.h
        src/model/FirstLastList.cpp
        src/model/FirstLastList.h
        src/model/STRAND.h
        src/model/REGION.h
        src/impl/impl.h
        src/impl/GetReverseComplementary.cpp
        src/impl/GetReverseComplementary.h
        src/impl/readFastaFile.cpp
        src/impl/readFastaFile.h
        src/impl/checkOrfState.cpp
        src/impl/checkOrfState.h
        src/util/util.h
        src/impl/nA2AA.cpp
        src/impl/nA2AA.h
        src/impl/getSubsequence.cpp
        src/impl/getSubsequence.h
        src/impl/readGffFile.cpp
        src/impl/readGffFile.h
        src/impl/coordinateLiftOver.cpp
        src/impl/coordinateLiftOver.h
        src/impl/getSequencesFromGff.cpp
        src/impl/getSequencesFromGff.h
        src/service/service.h
        src/myImportandFunction/myImportantFunction.h
        src/impl/TranscriptUpdateInformation.cpp
        src/impl/TranscriptUpdateInformation.h
        src/impl/FileWriteCatch.cpp
        src/impl/FileWriteCatch.h
        src/tests/util/nucleotideCodeSubstitutionMatrix_test.cpp
        src/model/Range.cpp src/model/Range.h src/model/AlignmentMatch.cpp
        src/model/AlignmentMatch.h
        src/service/TransferGffWithNucmerResult.cpp
        src/service/TransferGffWithNucmerResult.h
        src/myImportandFunction/alignSlidingWindow.cpp
        src/myImportandFunction/alignSlidingWindow.h
        src/impl/readGffFileWithEverything.cpp
        src/impl/readGffFileWithEverything.h
        src/tests/impl/readGffFileWithEverything.cpp
        src/impl/CheckAndUpdateTranscriptsEnds.cpp
        src/impl/CheckAndUpdateTranscriptsEnds.h
        src/tests/impl/readFastaFile_test.cpp
        src/impl/deNovoGenomeVariantCalling.cpp
        src/impl/deNovoGenomeVariantCalling.h
        src/model/OrthologPair.cpp
        src/model/OrthologPair.h
        src/impl/geneSyntenic.cpp
        src/impl/geneSyntenic.h

        WFA/system/profiler_counter.c
        WFA/system/profiler_counter.h
        WFA/system/profiler_timer.c
        WFA/system/profiler_timer.h
        WFA/benchmark/benchmark_edit.h
        WFA/benchmark/benchmark_gap_affine.c
        WFA/benchmark/benchmark_gap_affine.h
        WFA/benchmark/benchmark_gap_lineal.c
        WFA/benchmark/benchmark_gap_lineal.h
        WFA/benchmark/benchmark_utils.c
        WFA/benchmark/benchmark_utils.h
        WFA/utils/commons.c
        WFA/utils/commons.h
        WFA/utils/dna_text.c
        WFA/utils/dna_text.h
        WFA/utils/string_padded.c
        WFA/utils/string_padded.h
        WFA/utils/vector.c
        WFA/utils/vector.h
        WFA/system/mm_allocator.c
        WFA/system/mm_allocator.h
        WFA/gap_affine/affine_penalties.c
        WFA/gap_affine/affine_penalties.h
        WFA/gap_affine/affine_table.c
        WFA/gap_affine/affine_table.h
        WFA/gap_affine/affine_wavefront_align.c
        WFA/gap_affine/affine_wavefront_align.h
        WFA/gap_affine/affine_wavefront_backtrace.c
        WFA/gap_affine/affine_wavefront_backtrace.h
        WFA/gap_affine/affine_wavefront.c
        WFA/gap_affine/affine_wavefront_display.c
        WFA/gap_affine/affine_wavefront_display.h
        WFA/gap_affine/affine_wavefront_extend.c
        WFA/gap_affine/affine_wavefront_extend.h
        WFA/gap_affine/affine_wavefront.h
        WFA/gap_affine/affine_wavefront_penalties.c
        WFA/gap_affine/affine_wavefront_penalties.h
        WFA/gap_affine/affine_wavefront_reduction.c
        WFA/gap_affine/affine_wavefront_reduction.h
        WFA/gap_affine/affine_wavefront_utils.c
        WFA/gap_affine/affine_wavefront_utils.h
        WFA/gap_affine/swg.c
        WFA/gap_affine/swg.h
        WFA/gap_affine/wavefront_stats.c
        WFA/gap_affine/wavefront_stats.h
        WFA/gap_lineal/lineal_penalties.h
        WFA/gap_lineal/nw.c
        WFA/gap_lineal/nw.h
        WFA/edit/edit_cigar.c
        WFA/edit/edit_cigar.h
        WFA/edit/edit_dp.c
        WFA/edit/edit_dp.h
        WFA/edit/edit_table.c
        WFA/edit/edit_table.h
        WFA/benchmark/benchmark_edit.c
        WFA/benchmark/benchmark_edit.h
        WFA/benchmark/benchmark_gap_affine.c
        WFA/benchmark/benchmark_gap_affine.h
        WFA/benchmark/benchmark_gap_lineal.c
        WFA/benchmark/benchmark_gap_lineal.h
        WFA/benchmark/benchmark_utils.c
        WFA/benchmark/benchmark_utils.h
        )

target_link_libraries(proali gtest gtest_main)

