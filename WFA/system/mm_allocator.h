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
 * DESCRIPTION: Simple managed-memory allocator that reduces the overhead
 *   of using malloc/calloc/free functions by allocating slabs of memory
 *   and dispatching memory segments in order.
 */

#ifndef MM_ALLOCATOR_H_
#define MM_ALLOCATOR_H_

#include "../utils/vector.h"

/*
 * Debug
 */
//#define MM_ALLOCATOR_LOG

/*
 * MM-Allocator
 */
typedef struct {
  // Metadata
  uint64_t request_ticker;      // Request ticker
  // Memory segments
  uint64_t segment_size;        // Memory segment size (bytes)
  vector_t* segments;           // Memory segments (mm_allocator_segment_t*)
  vector_t* segments_free;      // Completely free segments (mm_allocator_segment_t*)
  uint64_t current_segment_idx; // Current segment being used (serving memory)
  // Malloc Memory
  vector_t* malloc_requests;    // Malloc requests (void*)
} mm_allocator_t;

/*
 * Setup
 */
#ifdef __cplusplus
extern "C" {
#endif
mm_allocator_t* mm_allocator_new( uint64_t segment_size);
#ifdef __cplusplus
}
#endif


void mm_allocator_clear(
    mm_allocator_t* const mm_allocator);

#ifdef __cplusplus
extern "C" {
#endif
void mm_allocator_delete(
    mm_allocator_t* const mm_allocator);
#ifdef __cplusplus
}
#endif
/*
 * Allocate
 */
void* mm_allocator_allocate(
    mm_allocator_t* const mm_allocator,
    uint64_t num_bytes,
    const bool zero_mem
#ifdef MM_ALLOCATOR_LOG
    ,const char* func_name,
    uint64_t line_no
#endif
    );

#ifdef MM_ALLOCATOR_LOG
#define mm_allocator_alloc(mm_allocator,type) \
  ((type*)mm_allocator_allocate(mm_allocator,sizeof(type),false,__func__,(uint64_t)__LINE__))
#define mm_allocator_malloc(mm_allocator,num_bytes) \
  (mm_allocator_allocate(mm_allocator,(num_bytes),false,__func__,(uint64_t)__LINE__))
#define mm_allocator_calloc(mm_allocator,num_elements,type,clear_mem) \
  ((type*)mm_allocator_allocate(mm_allocator,(num_elements)*sizeof(type),clear_mem,__func__,(uint64_t)__LINE__))
#else
#define mm_allocator_alloc(mm_allocator,type) \
  ((type*)mm_allocator_allocate(mm_allocator,sizeof(type),false))
#define mm_allocator_malloc(mm_allocator,num_bytes) \
  (mm_allocator_allocate(mm_allocator,(num_bytes),false))
#define mm_allocator_calloc(mm_allocator,num_elements,type,clear_mem) \
  ((type*)mm_allocator_allocate(mm_allocator,(num_elements)*sizeof(type),clear_mem))
#endif

#define mm_allocator_uint64(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint64_t))
#define mm_allocator_uint32(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint32_t))
#define mm_allocator_uint16(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint16_t))
#define mm_allocator_uint8(mm_allocator)  mm_allocator_malloc(mm_allocator,sizeof(uint8_t))

/*
 * Free
 */
void mm_allocator_free(
    mm_allocator_t* const mm_allocator,
    void* const memory);

/*
 * Utils
 */
void mm_allocator_get_occupation(
    mm_allocator_t* const mm_allocator,
    uint64_t* const bytes_used,
    uint64_t* const bytes_free_available,
    uint64_t* const bytes_free_fragmented);

/*
 * Display
 */
void mm_allocator_print(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool display_requests);

#endif /* MM_ALLOCATOR_H_ */
