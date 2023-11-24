cmake_minimum_required(VERSION 3.0)
project(anchorwave)

#include(CheckCCompilerFlag)
#include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)
set(CMAKE_CXX_STANDARD 14)
#
#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    set(MACOSX TRUE)
##else()
##    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
#endif()

find_package(ZLIB)

set(CMAKE_CXX_FLAGS " -lm -lz -Wall -g3 -O3 -funroll-all-loops -I ./ -std=gnu++14 -pthread -lrt") # the -mavx* is for GNU compiler and xCORE* are for intel compiler
set(CMAKE_C_FLAGS " -lm -lz -Wall -g3 -O3 -funroll-all-loops -I ./ -pthread -lrt")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512bw -mavx512cd -mavx512vl -mavx512f -mavx512dq")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512bw -mavx512cd -mavx512vl -mavx512f -mavx512dq")

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m512i a, b, c;
          const int64_t src[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
          int64_t dst[8];
          a =  _mm512_loadu_si512( (__m512i*)src );
          b =  _mm512_loadu_si512( (__m512i*)src );
          c = _mm512_add_epi64( a, b );
          _mm512_storeu_si512( (__m512i*)dst, c );
          for( int i = 0; i < 8; i++ ){
            if( ( src[i] + src[i] ) != dst[i] ){
              return -1;
            }
          }
          return 0;
        }"
        COMPILER_CAN_HANDLE_AVX512)


if (COMPILER_CAN_HANDLE_AVX512)
    message(STATUS "Support AVX512: True")
else()
    message(STATUS "Support AVX512: FALSE")
    string(REPLACE "-mavx512bw -mavx512cd -mavx512vl -mavx512f -mavx512dq" "" CUSTOM_FLAGS_CXX ${CMAKE_CXX_FLAGS})
    string(REPLACE "-mavx512bw -mavx512cd -mavx512vl -mavx512f -mavx512dq" "" CUSTOM_FLAGS_C ${CMAKE_C_FLAGS})
    set(CMAKE_C_FLAGS "${CUSTOM_FLAGS_C}")
    set(CMAKE_CXX_FLAGS "${CUSTOM_FLAGS_CXX}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m256i a, b, c;
          const int32_t src[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
          int32_t dst[8];
          a =  _mm256_loadu_si256( (__m256i*)src );
          b =  _mm256_loadu_si256( (__m256i*)src );
          c = _mm256_add_epi32( a, b );
          _mm256_storeu_si256( (__m256i*)dst, c );
          for( int i = 0; i < 8; i++ ){
            if( ( src[i] + src[i] ) != dst[i] ){
              return -1;
            }
          }
          return 0;
        }"
        COMPILER_CAN_HANDLE_AVX2)

if (COMPILER_CAN_HANDLE_AVX2)
    message(STATUS "Support AVX2: True")
else()
    message(STATUS "Support AVX2: FALSE")
    string(REPLACE "-mavx2" "" CUSTOM_FLAGS_CXX ${CMAKE_CXX_FLAGS})
    string(REPLACE "-mavx2" "" CUSTOM_FLAGS_C ${CMAKE_C_FLAGS})
    set(CMAKE_C_FLAGS "${CUSTOM_FLAGS_C}")
    set(CMAKE_CXX_FLAGS "${CUSTOM_FLAGS_CXX}")
endif()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m128i a, b, c, d;
          const int8_t src1[16] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
          const int8_t src2[16] = {  16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
          int8_t dstmin[16];
          int8_t dstmax[16];
          a =  _mm_loadu_si128 ( (__m128i*)src1 );
          b =  _mm_loadu_si128 ( (__m128i*)src2 );
          c = _mm_max_epi8( a, b );
          d = _mm_min_epi8( a, b );
          _mm_storeu_si128( (__m128i*)dstmax, c );
          _mm_storeu_si128( (__m128i*)dstmin, d );
          for( int i = 0; i < 16; i++ ){
            if( ( dstmin[i] + dstmax[i] ) != 17 ){
              return -1;
            }
          }
          return 0;
        }"
        COMPILER_CAN_HANDLE_SSE4.1)

if (COMPILER_CAN_HANDLE_SSE4.1)
    message(STATUS "Support SSE4.1: True")
else()
    message(STATUS "Support SSE4.1: FALSE")
    string(REPLACE "-msse4.1" "" CUSTOM_FLAGS_CXX ${CMAKE_CXX_FLAGS})
    string(REPLACE "-msse4.1" "" CUSTOM_FLAGS_C ${CMAKE_C_FLAGS})
    set(CMAKE_C_FLAGS "${CUSTOM_FLAGS_C}")
    set(CMAKE_CXX_FLAGS "${CUSTOM_FLAGS_CXX}")
endif()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")


check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m128i a, b, c;
          const int16_t src[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
          int16_t dst[8];
          a =  _mm_loadu_si128 ( (__m128i*)src );
          b =  _mm_loadu_si128 ( (__m128i*)src );
          c = _mm_add_epi16( a, b );
          _mm_storeu_si128( (__m128i*)dst, c );
          for( int i = 0; i < 8; i++ ){
            if( ( src[i] + src[i] ) != dst[i] ){
              return -1;
            }
          }
          return 0;
        }"
        COMPILER_CAN_HANDLE_SSE2)

if (COMPILER_CAN_HANDLE_SSE2)
    message(STATUS "Support SSE2: True")
else()
    message(STATUS "Support SSE2: FALSE")
    string(REPLACE "-msse2" "" CUSTOM_FLAGS_CXX ${CMAKE_CXX_FLAGS})
    string(REPLACE "-msse2" "" CUSTOM_FLAGS_C ${CMAKE_C_FLAGS})
    set(CMAKE_C_FLAGS "${CUSTOM_FLAGS_C}")
    set(CMAKE_CXX_FLAGS "${CUSTOM_FLAGS_CXX}")
endif()

EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
if(${ARCHITECTURE} STREQUAL "arm64")
    add_definitions(-D__SSE2__ -D__SSE2NEON__)
    include_directories(src minimap2/sse2neon)
    set(CMAKE_CXX_FLAGS " -Wall -g3 -O3 -I ./ -std=gnu++14 -pthread") # the -mavx* is for GNU compiler and xCORE* are for intel compiler
    set(CMAKE_C_FLAGS " -Wall -g3 -O3 -I ./ -pthread")
elseif (${ARCHITECTURE} STREQUAL "x86_64" AND CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_CXX_FLAGS " -lm -lz -Wall -g3 -O3 -funroll-all-loops -msse4.1 -msse4.2 -I ./ -std=gnu++14 -pthread") # the -mavx* is for GNU compiler and xCORE* are for intel compiler
    set(CMAKE_C_FLAGS " -lm -lz -Wall -g3 -O3 -funroll-all-loops -msse4.1 -msse4.2 -I ./ -pthread")
endif()

message(STATUS "${CMAKE_C_FLAGS}")

#for google unit test
add_subdirectory(./googletest/googletest)
include_directories(googletest/googletest/include googletest/googletest)

add_executable(anchorwave
        googletest
        googletest/googletest
        sse2neon.h
        main.cpp

        src/controlLayer.cpp
        src/controlLayer.h
        src/InputParser.cpp
        src/InputParser.h
        src/version.h

        src/impl/deNovoGenomeVariantCallingMultipleThreads.cpp
        src/impl/deNovoGenomeVariantCalling.h
        src/impl/geneSyntenic.cpp
        src/impl/geneSyntenic.h
        src/impl/GetReverseComplementary.cpp
        src/impl/GetReverseComplementary.h
        src/impl/getSequencesFromGff.cpp
        src/impl/getSequencesFromGff.h
        src/impl/getSubsequence.cpp
        src/impl/getSubsequence.h
        src/impl/readFastaFile.cpp
        src/impl/readFastaFile.h
        src/impl/readGffFile.cpp
        src/impl/readGffFile.h
        src/impl/TranscriptUpdateInformation.cpp
        src/impl/TranscriptUpdateInformation.h

        src/model/AlignmentMatch.cpp
        src/model/AlignmentMatch.h
        src/model/GenomeBasicFeature.cpp
        src/model/GenomeBasicFeature.h
        src/model/Score.cpp
        src/model/Score.h
        src/model/STRAND.cpp
        src/model/STRAND.h
        src/model/Transcript.cpp
        src/model/Transcript.h

        src/myImportandFunction/alignSlidingWindow.cpp
        src/myImportandFunction/alignSlidingWindow.h

        src/service/TransferGffWithNucmerResult.cpp
        src/service/TransferGffWithNucmerResult.h

        src/util/myutil.cpp
        src/util/myutil.h


        WFA2-lib/alignment/affine2p_penalties.c
        WFA2-lib/alignment/affine_penalties.c
        WFA2-lib/alignment/cigar.c
        WFA2-lib/alignment/score_matrix.c
        WFA2-lib/system/mm_allocator.c
        WFA2-lib/system/mm_stack.c
        WFA2-lib/system/profiler_counter.c
        WFA2-lib/system/profiler_timer.c
        WFA2-lib/utils/bitmap.c
        WFA2-lib/utils/commons.c
        WFA2-lib/utils/dna_text.c
        WFA2-lib/utils/heatmap.c
        WFA2-lib/utils/sequence_buffer.c
        WFA2-lib/utils/string_padded.c
        WFA2-lib/utils/vector.c
        WFA2-lib/wavefront/wavefront_align.c
        WFA2-lib/wavefront/wavefront_aligner.c
        WFA2-lib/wavefront/wavefront_attributes.c
        WFA2-lib/wavefront/wavefront_backtrace_buffer.c
        WFA2-lib/wavefront/wavefront_backtrace.c
        WFA2-lib/wavefront/wavefront_backtrace_offload.c
        WFA2-lib/wavefront/wavefront_bialign.c
        WFA2-lib/wavefront/wavefront.c
        WFA2-lib/wavefront/wavefront_components.c
        WFA2-lib/wavefront/wavefront_compute_affine2p.c
        WFA2-lib/wavefront/wavefront_compute_affine.c
        WFA2-lib/wavefront/wavefront_compute.c
        WFA2-lib/wavefront/wavefront_compute_edit.c
        WFA2-lib/wavefront/wavefront_compute_linear.c
        WFA2-lib/wavefront/wavefront_debug.c
        WFA2-lib/wavefront/wavefront_display.c
        WFA2-lib/wavefront/wavefront_extend.c
        WFA2-lib/wavefront/wavefront_heuristic.c
        WFA2-lib/wavefront/wavefront_pcigar.c
        WFA2-lib/wavefront/wavefront_penalties.c
        WFA2-lib/wavefront/wavefront_plot.c
        WFA2-lib/wavefront/wavefront_slab.c
        WFA2-lib/wavefront/wavefront_unialign.c
        WFA2-lib/wavefront/wavefront_unialign.h
        WFA2-lib/wavefront/wavefront_bialigner.c
        WFA2-lib/wavefront/wavefront_bialigner.h
        WFA2-lib/alignment/affine2p_penalties.h
        WFA2-lib/alignment/affine_penalties.h
        WFA2-lib/alignment/cigar.h
        WFA2-lib/alignment/linear_penalties.h
        WFA2-lib/alignment/score_matrix.h
        WFA2-lib/scripts/fasta2seq.sh
        WFA2-lib/scripts/seq2fasta.sh
        WFA2-lib/scripts/wfa.alg.cmp.score.sh
        WFA2-lib/system/mm_allocator.h
        WFA2-lib/system/mm_stack.h
        WFA2-lib/system/profiler_counter.h
        WFA2-lib/system/profiler_timer.h
        WFA2-lib/utils/bitmap.h
        WFA2-lib/utils/commons.h
        WFA2-lib/utils/dna_text.h
        WFA2-lib/utils/heatmap.h
        WFA2-lib/utils/sequence_buffer.h
        WFA2-lib/utils/string_padded.h
        WFA2-lib/utils/vector.h
        WFA2-lib/wavefront/wavefront_aligner.h
        WFA2-lib/wavefront/wavefront_align.h
        WFA2-lib/wavefront/wavefront_attributes.h
        WFA2-lib/wavefront/wavefront_backtrace_buffer.h
        WFA2-lib/wavefront/wavefront_backtrace.h
        WFA2-lib/wavefront/wavefront_backtrace_offload.h
        WFA2-lib/wavefront/wavefront_bialign.h
        WFA2-lib/wavefront/wavefront_components.h
        WFA2-lib/wavefront/wavefront_compute_affine2p.h
        WFA2-lib/wavefront/wavefront_compute_affine.h
        WFA2-lib/wavefront/wavefront_compute_edit.h
        WFA2-lib/wavefront/wavefront_compute.h
        WFA2-lib/wavefront/wavefront_compute_linear.h
        WFA2-lib/wavefront/wavefront_debug.h
        WFA2-lib/wavefront/wavefront_display.h
        WFA2-lib/wavefront/wavefront_extend.h
        WFA2-lib/wavefront/wavefront.h
        WFA2-lib/wavefront/wavefront_heuristic.h
        WFA2-lib/wavefront/wavefront_offset.h
        WFA2-lib/wavefront/wavefront_pcigar.h
        WFA2-lib/wavefront/wavefront_penalties.h
        WFA2-lib/wavefront/wavefront_plot.h
        WFA2-lib/wavefront/wavefront_slab.h
        WFA2-lib/bindings/cpp/WFAligner.cpp
        WFA2-lib/bindings/cpp/WFAligner.hpp


        minimap2/align.c
        minimap2/bseq.c
        minimap2/chain.c
        minimap2/esterr.c
        minimap2/format.c
        minimap2/hit.c
        minimap2/index.c
        minimap2/kalloc.c
        minimap2/ksw2_dispatch.c
        minimap2/ksw2_extd2_sse.c
        minimap2/ksw2_extd2_avx.c
        minimap2/ksw2_exts2_sse.c
        minimap2/ksw2_extz2_sse.c
        minimap2/ksw2_ll_sse.c
        minimap2/kthread.c
        minimap2/map.c
        minimap2/misc.c
        minimap2/options.c
        minimap2/pe.c
        minimap2/sdust.c
        minimap2/sketch.c
        minimap2/splitidx.c
        minimap2/sse2neon/emmintrin.h
        minimap2/bseq.h
        minimap2/kalloc.h
        minimap2/kdq.h
        minimap2/ketopt.h
        minimap2/khash.h
        minimap2/kseq.h
        minimap2/ksort.h
        minimap2/ksw2.h
        minimap2/kthread.h
        minimap2/kvec.h
        minimap2/minimap.h
        minimap2/mmpriv.h
        minimap2/sdust.h
        )

target_link_libraries(anchorwave gtest gtest_main ZLIB::ZLIB)

install(TARGETS anchorwave DESTINATION bin)
