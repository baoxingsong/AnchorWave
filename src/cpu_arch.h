#ifndef ANCHORWAVE_CPU_ARCH_H
#define ANCHORWAVE_CPU_ARCH_H

/*
 * Use compiler-defined macros so architecture selection also works for cross
 * builds and Apple universal builds.  Minimap2's SIMD sources use the SSE2
 * interface; on ARM, SSE2NEON supplies that interface using NEON instructions.
 */
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64) || \
        defined(__arm__) || defined(_M_ARM)
#define ANCHORWAVE_ARCH_ARM 1
#ifndef __SSE2NEON__
#define __SSE2NEON__ 1
#endif
#ifndef __SSE2__
#define __SSE2__ 1
#endif
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#define ANCHORWAVE_ARCH_X86 1
#else
#error "AnchorWave requires an x86 CPU with SSE2 or an ARM CPU with NEON"
#endif

#endif
