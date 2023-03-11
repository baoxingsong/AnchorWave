//
// Created by song on 8/5/18.
//

#pragma once

#ifdef __SSE2NEON__
#include "../../sse2neon.h"
#else
#include <immintrin.h>
#endif // __SSE2NEON__

#include "../impl/getSequencesFromGff.h"
#include "../model//Score.h"
#include "../util/myutil.h"

#include "../../WFA2-lib/bindings/cpp/WFAligner.hpp"
#include "../../minimap2/ksw2.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <string>

int64_t alignSlidingWindow(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                           const int64_t &slidingWindowSize, const int32_t &matchingScore,
                           const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2);

int64_t alignSlidingWindowNW(std::string &dna_q, std::string &dna_d, std::string &_alignment_q, std::string &_alignment_d,
                             const int64_t &slidingWindowSize, const int32_t &matchingScore,
                             const int32_t &mismatchingPenalty, const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2);
