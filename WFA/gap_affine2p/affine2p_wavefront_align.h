/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WaveFront-Alignment for pairwise gap-affine 2-pieces alignment
 */

#ifndef AFFINE2P_WAVEFRONT_ALIGN_H_
#define AFFINE2P_WAVEFRONT_ALIGN_H_

#include "affine2p_wavefront.h"

/*
 * Computation using Wavefronts
 */
#ifdef __cplusplus
extern "C" {
#endif
void affine2p_wavefronts_align(
    affine2p_wavefronts_t* const affine2p_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
int64_t affine2p_wavefronts_align_anchorwavenumberthreshold_version(
        affine2p_wavefronts_t* const affine2p_wavefronts,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length, const int64_t anchorwavenumberthreshold);
#ifdef __cplusplus
}
#endif

#endif /* AFFINE2P_WAVEFRONT_ALIGN_H_ */
