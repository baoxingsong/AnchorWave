cmake_minimum_required(VERSION 3.10)
project(anchorwave VERSION 1.3.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# Inner-wavefront parallelism is dormant at runtime unless an attempt sets
# max_num_threads > 1. AnchorWave keeps the default at one and only grants
# 2/4/8 threads to selected tail tasks, avoiding ordinary oversubscription.
option(ANCHORWAVE_ENABLE_WFA_PARALLEL
       "Compile OpenMP support for per-attempt WFA tail parallelism" ON)
if(ANCHORWAVE_ENABLE_WFA_PARALLEL)
    find_package(OpenMP QUIET)
endif()

# The source code selects the x86 intrinsic headers or the SSE2NEON adapter from
# compiler-provided architecture macros.  For a native x86 build, let the
# compiler enable exactly the SIMD extensions supported by the build CPU.  The
# explicit targets below are for reproducible SIMD regression builds; a cross
# build deliberately keeps the toolchain's target settings unchanged.
set(ANCHORWAVE_SIMD_TARGET "native" CACHE STRING
    "x86 SIMD target: native, avx2, sse4, sse2, or toolchain")
set_property(CACHE ANCHORWAVE_SIMD_TARGET PROPERTY STRINGS
             native avx2 sse4 sse2 toolchain)
set(ANCHORWAVE_SIMD_COMPILE_OPTIONS)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$" AND
   CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" AND
   CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(ANCHORWAVE_SIMD_TARGET STREQUAL "native")
        if(NOT CMAKE_CROSSCOMPILING)
            list(APPEND ANCHORWAVE_SIMD_COMPILE_OPTIONS -march=native)
            message(STATUS "AnchorWave CPU architecture: ${CMAKE_SYSTEM_PROCESSOR} (native SIMD auto-detection)")
        else()
            message(STATUS "AnchorWave CPU architecture: ${CMAKE_SYSTEM_PROCESSOR} (cross-toolchain SIMD flags)")
        endif()
    elseif(ANCHORWAVE_SIMD_TARGET STREQUAL "avx2")
        list(APPEND ANCHORWAVE_SIMD_COMPILE_OPTIONS -march=haswell)
        message(STATUS "AnchorWave x86 SIMD target: AVX2")
    elseif(ANCHORWAVE_SIMD_TARGET STREQUAL "sse4")
        # Some supported GCC toolchains do not accept the Penryn alias even
        # though they support SSE4.1.  Express the same ISA floor explicitly.
        list(APPEND ANCHORWAVE_SIMD_COMPILE_OPTIONS -march=core2 -msse4.1)
        message(STATUS "AnchorWave x86 SIMD target: SSE4.1")
    elseif(ANCHORWAVE_SIMD_TARGET STREQUAL "sse2")
        list(APPEND ANCHORWAVE_SIMD_COMPILE_OPTIONS -march=x86-64)
        message(STATUS "AnchorWave x86 SIMD target: SSE2")
    elseif(NOT ANCHORWAVE_SIMD_TARGET STREQUAL "toolchain")
        message(FATAL_ERROR "Unsupported ANCHORWAVE_SIMD_TARGET=${ANCHORWAVE_SIMD_TARGET}")
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64|arm|ARM)")
    message(STATUS "AnchorWave CPU architecture: ${CMAKE_SYSTEM_PROCESSOR} (NEON via SSE2NEON)")
else()
    message(STATUS "AnchorWave CPU architecture: ${CMAKE_SYSTEM_PROCESSOR} (toolchain defaults)")
endif()


set(ANCHORWAVE_TRIO_SOURCES
        src/trio/TrioGraphCommand.cpp
        src/trio/Phylogeny.cpp
        src/trio/io/TaxonManifestReader.cpp
        src/trio/io/PairwiseManifestReader.cpp
        src/trio/io/CopyConstraintReader.cpp
        src/trio/io/BlockProjectionReader.cpp
        src/trio/io/FastaEvidenceValidator.cpp
        src/trio/io/MafEvidenceReader.cpp
        src/trio/io/GraphWriters.cpp
        src/trio/model/TrioTypes.cpp
        src/trio/model/StableId.cpp
        src/trio/model/AlignmentEvidence.cpp
        src/trio/model/AlignmentSiteGraph.cpp
        src/trio/impl/CopyRelationshipResolver.cpp
        src/trio/alignment/TwoPieceAffine.cpp
        src/trio/impl/LocalizedPoaRepair.cpp
        src/trio/impl/GraphRepairEngine.cpp
        src/trio/evolution/AncestralInference.cpp
        src/trio/evolution/AncestralAdjacency.cpp)

add_library(anchorwave_trio_core STATIC ${ANCHORWAVE_TRIO_SOURCES})
target_include_directories(anchorwave_trio_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(anchorwave_trio_core PRIVATE -Wall -Wextra -pedantic)
target_link_libraries(anchorwave_trio_core PUBLIC ZLIB::ZLIB)

# Keep the vendored WFA2 sources in a dedicated target.  This list mirrors
# smarco/WFA2-lib main@2fc4a1afac0f624e7020ee5aadb7692b38157eaa and makes it
# possible to unit-test AnchorWave's BiWFA adapter without linking the complete
# genoAli/proali executable.
set(ANCHORWAVE_WFA2_SOURCES
        WFA2-lib/wavefront/wavefront_align.c
        WFA2-lib/wavefront/wavefront_aligner.c
        WFA2-lib/wavefront/wavefront_attributes.c
        WFA2-lib/wavefront/wavefront_backtrace_buffer.c
        WFA2-lib/wavefront/wavefront_backtrace.c
        WFA2-lib/wavefront/wavefront_backtrace_offload.c
        WFA2-lib/wavefront/wavefront_bialign.c
        WFA2-lib/wavefront/wavefront_bialigner.c
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
        WFA2-lib/wavefront/wavefront_extend_kernels.c
        WFA2-lib/wavefront/wavefront_extend_kernels_avx.c
        WFA2-lib/wavefront/wavefront_heuristic.c
        WFA2-lib/wavefront/wavefront_pcigar.c
        WFA2-lib/wavefront/wavefront_penalties.c
        WFA2-lib/wavefront/wavefront_plot.c
        WFA2-lib/wavefront/wavefront_sequences.c
        WFA2-lib/wavefront/wavefront_slab.c
        WFA2-lib/wavefront/wavefront_termination.c
        WFA2-lib/wavefront/wavefront_unialign.c
        WFA2-lib/system/mm_stack.c
        WFA2-lib/system/mm_allocator.c
        WFA2-lib/system/profiler_counter.c
        WFA2-lib/system/profiler_timer.c
        WFA2-lib/utils/bitmap.c
        WFA2-lib/utils/dna_text.c
        WFA2-lib/utils/sequence_buffer.c
        WFA2-lib/utils/vector.c
        WFA2-lib/utils/commons.c
        WFA2-lib/utils/heatmap.c
        WFA2-lib/alignment/affine2p_penalties.c
        WFA2-lib/alignment/affine_penalties.c
        WFA2-lib/alignment/cigar.c
        WFA2-lib/alignment/cigar_utils.c
        WFA2-lib/alignment/score_matrix.c)

add_library(anchorwave_wfa2_core STATIC
        ${ANCHORWAVE_WFA2_SOURCES}
        src/myImportandFunction/AlignmentAlgorithmSelector.cpp
        src/myImportandFunction/AlignmentResourceScheduler.cpp
        src/myImportandFunction/WfaAlignment.cpp)
target_include_directories(anchorwave_wfa2_core PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/WFA2-lib")
target_compile_options(anchorwave_wfa2_core PRIVATE
        -Wall
        -O3
        -funroll-loops
        ${ANCHORWAVE_SIMD_COMPILE_OPTIONS})
if(ANCHORWAVE_ENABLE_WFA_PARALLEL AND OpenMP_C_FOUND)
    target_compile_definitions(anchorwave_wfa2_core PRIVATE WFA_PARALLEL)
    target_link_libraries(anchorwave_wfa2_core PUBLIC OpenMP::OpenMP_C)
    message(STATUS "AnchorWave inner-WFA OpenMP support: enabled (runtime default: 1 thread)")
elseif(ANCHORWAVE_ENABLE_WFA_PARALLEL)
    message(STATUS "AnchorWave inner-WFA OpenMP support: unavailable; max_num_threads remains advisory")
endif()
if(UNIX)
    target_link_libraries(anchorwave_wfa2_core PUBLIC m Threads::Threads)
endif()

add_library(anchorwave_anchor_scheduler STATIC
        src/service/AnchorTaskExecutor.cpp
        src/service/AlignmentGapPlanner.cpp)
target_include_directories(anchorwave_anchor_scheduler PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(anchorwave_anchor_scheduler PRIVATE
        -Wall -Wextra -pedantic)
target_link_libraries(anchorwave_anchor_scheduler PUBLIC Threads::Threads)

add_library(anchorwave_compressed_input STATIC
        src/io/CompressedInput.cpp)
target_include_directories(anchorwave_compressed_input PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(anchorwave_compressed_input PRIVATE
        -Wall -Wextra -pedantic)
target_link_libraries(anchorwave_compressed_input PUBLIC ZLIB::ZLIB)

option(ANCHORWAVE_ENABLE_SANITIZERS "Enable ASan and UBSan for core targets" OFF)
if(ANCHORWAVE_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(anchorwave_trio_core PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_libraries(anchorwave_trio_core PUBLIC -fsanitize=address,undefined)
    target_compile_options(anchorwave_wfa2_core PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_libraries(anchorwave_wfa2_core PUBLIC -fsanitize=address,undefined)
    target_compile_options(anchorwave_anchor_scheduler PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_libraries(anchorwave_anchor_scheduler PUBLIC
            -fsanitize=address,undefined)
    target_compile_options(anchorwave_compressed_input PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_libraries(anchorwave_compressed_input PUBLIC
            -fsanitize=address,undefined)
endif()

add_executable(anchorwave
        main.cpp

        src/controlLayer.cpp
        src/controlLayer.h
        src/InputParser.cpp
        src/InputParser.h
        src/version.h

        src/impl/deNovoGenomeVariantCallingMultipleThreads.cpp
        src/impl/deNovoGenomeVariantCalling.h
        src/impl/AlignmentBlockBuffer.cpp
        src/impl/AlignmentBlockBuffer.h
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
        src/service/ReadSamUtils.h

        src/util/myutil.cpp
        src/util/myutil.h


        minimap2/align.c
        minimap2/bseq.c
        minimap2/chain.c
        minimap2/esterr.c
        minimap2/format.c
        minimap2/hit.c
        minimap2/index.c
        minimap2/kalloc.c
        minimap2/ksw2_progress.c
        minimap2/ksw2_dispatch.c
        minimap2/ksw2_extd2_sse.c
        minimap2/ksw2_extd2_avx.c
        minimap2/ksw2_singletrack.cpp
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

target_include_directories(anchorwave PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(anchorwave PRIVATE
        -Wall
        -O3
        -funroll-loops
        ${ANCHORWAVE_SIMD_COMPILE_OPTIONS})

target_link_libraries(anchorwave PRIVATE anchorwave_trio_core anchorwave_wfa2_core
        anchorwave_anchor_scheduler anchorwave_compressed_input
        ZLIB::ZLIB Threads::Threads
        ${CMAKE_DL_LIBS})
if(UNIX AND NOT APPLE)
    target_link_libraries(anchorwave PRIVATE rt m)
elseif(UNIX)
    target_link_libraries(anchorwave PRIVATE m)
endif()

install(TARGETS anchorwave DESTINATION bin)
