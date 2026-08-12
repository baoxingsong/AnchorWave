#include "AlignmentAlgorithmSelector.h"

#include "WfaAlignment.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <climits>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace anchorwave {
namespace {

constexpr uint64_t kKmerLength = 15;
constexpr uint64_t kCertificationAnchorSpacing = 16384;
constexpr uint64_t kTargetSketchSize = 4096;
constexpr uint64_t kMaximumSketchSize = 8192;
constexpr uint64_t kMaximumOccurrencesPerHash = 4;
constexpr uint64_t kSingletrackWfaBaseBytes = 32ULL * 1024ULL * 1024ULL;
constexpr uint64_t kStandardWfaBaseBytes = 16ULL * 1024ULL * 1024ULL;
constexpr uint64_t kSuccinctWfaBaseBytes = 72ULL * 1024ULL * 1024ULL;
constexpr uint64_t kStandardWfaBytesPerScoreSquared = 32;
// Singletrack retains only M wavefronts for traceback. A conservative 8-byte
// coefficient is used instead of the theoretical dual-affine 32/5 ratio, and
// the fixed-width recycled I/D pool is accounted for separately below.
constexpr uint64_t kSingletrackWfaBytesPerScoreSquared = 8;
constexpr uint64_t kSingletrackOffsetBytes = 4;
constexpr uint64_t kSingletrackDiagonalPadding = 256;
// Conservative admission coefficients for WFA2's progressively more compact
// piggyback/backtrace modes. WFA2's runtime max_memory_abort guard remains the
// final authority if an empirical estimate is low.
constexpr uint64_t kMediumWfaBytesPerScoreSquared = 16;
constexpr uint64_t kLowWfaBytesPerScoreSquared = 8;
// ksw_extd2 allocates a rectangular virtual traceback matrix but writes only
// the SIMD-rounded live prefix of every anti-diagonal row. RSS admission must
// therefore count touched OS pages, not virtual bytes. The platform SIMD width
// below mirrors the compile-time dispatch in alignSlidingWindow.cpp.
#if defined(__AVX512BW__)
constexpr uint64_t kKsw2SimdBytes = 64;
#elif defined(__AVX2__)
constexpr uint64_t kKsw2SimdBytes = 32;
#else
constexpr uint64_t kKsw2SimdBytes = 16;
#endif
// Initial relative-time coefficients. Short/medium seeded benchmarks show a
// substantial succinct-mode setup/backtrace cost. Singletrack is refined
// below with deterministic strata fitted to completed B73/Mo17 records;
// keeping the constants
// analytic preserves monotonic behavior on other user-provided -t/-M limits.
constexpr long double kMediumWfaWorkMultiplier = 20.0L;
constexpr long double kLowWfaWorkMultiplier = 25.0L;
constexpr long double kStandardWfaWorkMultiplier = 1.25L;
// Cold-start calibration.  Runtime calibration is deliberately expressed as
// engine-specific analytic throughput rather than a hardware-specific length
// threshold.  The coefficients are refined below from the B73/Mo17 trace but
// retain monotonic scaling when users change -t and -M.
constexpr long double kWfaLinearBasesPerSecond = 4300000.0L;
constexpr long double kWfaScoreSquaredPerSecond = 460000000.0L;
// Isolated AVX512 B73/Mo17 matrices consistently sustain about 1.20--1.23e9
// cells/s after traceback and output construction.  KSW2-full, banded KSW2 and
// each sliding chunk use the same kernel; assigning sliding a synthetic 5.2e9
// rate hid long low-memory work from LPT and was a major source of tail delay.
constexpr long double kKsw2CellsPerSecond = 1200000000.0L;
constexpr long double kSlidingCellsPerSecond = 1200000000.0L;
constexpr double kWfaSetupSeconds = 0.00025;
constexpr double kKsw2SetupSeconds = 0.00015;
constexpr double kSlidingSetupSeconds = 0.00040;
// In the first complete 100-Mb B73/Mo17 trace, 8/31 rescue intervals were
// certified by the sparse exact-anchor chain before monolithic traceback.
// Among the remaining intervals the smallest planned band never uniquely
// certified, so execution starts at 2x and retains geometric widening.
constexpr long double kCertifiedChainSuccessProbability = 0.25L;

// B73/Mo17 Singletrack strata.  The old one-size-fits-all estimate was far
// above the median for easy intervals but missed a small difficult tail.  A
// score-stratified P50 makes fast high-WFA work visible to the scheduler;
// P90 uses both the memory-tail score and a slower measured throughput.
constexpr long double kSingletrackP50ScoreZeroMultiplier = 0.15L;
constexpr long double kSingletrackP50SmallScoreMultiplier = 0.25L;
constexpr long double kSingletrackP50MediumScoreMultiplier = 0.50L;
constexpr long double kSingletrackP50LargeScoreMultiplier = 0.55L;
constexpr long double kSingletrackP50UncertaintyMultiplier = 0.15L;
constexpr long double kSingletrackP50MaximumMultiplier = 0.85L;
constexpr long double kSingletrackP90EstimatedScoreMultiplier = 1.25L;
constexpr long double kSingletrackP90ConservativeScoreMultiplier = 0.65L;
constexpr long double kSingletrackP90ScoreSquaredPerSecond = 300000000.0L;


std::atomic<double> exactAlignmentMaximumEstimatedMinutes(0.0);

std::atomic<uint64_t> evaluatedIntervals(0);
std::atomic<uint64_t> exactTierIntervals(0);
std::atomic<uint64_t> bandedOnlyIntervals(0);
std::atomic<uint64_t> slidingOnlyIntervals(0);
std::atomic<uint64_t> singletrackWfaMemoryRejects(0);
std::atomic<uint64_t> singletrackWfaWorkWarnings(0);
std::atomic<uint64_t> singletrackWfaTimeRejects(0);
std::atomic<uint64_t> standardWfaMemoryRejects(0);
std::atomic<uint64_t> standardWfaWorkWarnings(0);
std::atomic<uint64_t> standardWfaTimeRejects(0);
std::atomic<uint64_t> mediumWfaMemoryRejects(0);
std::atomic<uint64_t> mediumWfaWorkWarnings(0);
std::atomic<uint64_t> mediumWfaTimeRejects(0);
std::atomic<uint64_t> lowWfaMemoryRejects(0);
std::atomic<uint64_t> lowWfaWorkWarnings(0);
std::atomic<uint64_t> lowWfaTimeRejects(0);
std::atomic<uint64_t> scoreCertifiedKsw2MemoryRejects(0);
std::atomic<uint64_t> scoreCertifiedKsw2TimeRejects(0);
std::atomic<uint64_t> fullKsw2MemoryRejects(0);
std::atomic<uint64_t> fullKsw2TimeRejects(0);
std::atomic<uint64_t> bandedKsw2MemoryRejects(0);
std::atomic<uint64_t> exactRuntimeMemoryFailures(0);
std::atomic<uint64_t> exactRuntimeOtherFailures(0);
std::atomic<uint64_t> bandedFallbackExecutions(0);
std::atomic<uint64_t> slidingFallbackExecutions(0);
std::atomic<uint64_t> nextTraceIntervalId(1);
std::mutex traceMutex;
std::ofstream traceStream;
std::atomic<bool> traceEnabled(false);
uint64_t traceLinesSinceFlush = 0;

struct SketchEntry {
    uint64_t hash = 0;
    uint32_t position = 0;
};

struct SequenceSummary {
    std::vector<SketchEntry> sketch;
    double ambiguousFraction = 0.0;
    double lowComplexityFraction = 0.0;
};

struct Anchor {
    uint32_t query = 0;
    uint32_t reference = 0;
};

uint64_t saturatingAdd(uint64_t first, uint64_t second) {
    if (second > std::numeric_limits<uint64_t>::max() - first) {
        return std::numeric_limits<uint64_t>::max();
    }
    return first + second;
}

uint64_t saturatingMultiply(uint64_t first, uint64_t second) {
    if (first == 0 || second == 0) {
        return 0;
    }
    if (first > std::numeric_limits<uint64_t>::max() / second) {
        return std::numeric_limits<uint64_t>::max();
    }
    return first * second;
}

uint64_t hostPageSizeBytes() {
#if defined(_WIN32)
    SYSTEM_INFO information;
    GetSystemInfo(&information);
    return std::max<uint64_t>(4096, information.dwPageSize);
#else
    const long pageSize = sysconf(_SC_PAGESIZE);
    return pageSize > 0 ? static_cast<uint64_t>(pageSize) : 4096;
#endif
}

uint64_t divideRoundUp(uint64_t value, uint64_t divisor) {
    return value == 0 ? 0 : 1 + (value - 1) / divisor;
}

uint64_t pageRoundedBytes(uint64_t bytes, uint64_t pageSize) {
    return saturatingMultiply(divideRoundUp(bytes, pageSize), pageSize);
}

int64_t floorDivideByTwo(int64_t value) {
    return value >= 0 ? value / 2 : -static_cast<int64_t>(
            (static_cast<uint64_t>(-(value + 1)) + 2) / 2);
}

uint64_t ksw2TouchedTracePages(uint64_t queryLength,
                               uint64_t referenceLength,
                               uint64_t bandWidth,
                               uint64_t pageSize) {
    if (queryLength == 0 || referenceLength == 0) {
        return 0;
    }
    const uint64_t diagonalCount = queryLength + referenceLength - 1;
    const uint64_t shorterLength = std::min(queryLength, referenceLength);
    const uint64_t effectiveBand = std::min(
            bandWidth, std::max(queryLength, referenceLength));
    const uint64_t retainedBases = std::min(
            shorterLength, saturatingAdd(effectiveBand, 1));
    const uint64_t rowStrideBytes = saturatingMultiply(
            saturatingAdd(divideRoundUp(retainedBases, kKsw2SimdBytes), 1),
            kKsw2SimdBytes);
    uint64_t pageCount = 0;
    uint64_t lastPage = 0;
    bool havePage = false;
    for (uint64_t diagonal = 0; diagonal < diagonalCount; ++diagonal) {
        int64_t start = std::max<int64_t>(
                0, static_cast<int64_t>(diagonal) -
                           static_cast<int64_t>(queryLength) + 1);
        int64_t end = std::min<int64_t>(
                static_cast<int64_t>(referenceLength) - 1,
                static_cast<int64_t>(diagonal));
        const int64_t leftBand = floorDivideByTwo(
                static_cast<int64_t>(diagonal) -
                static_cast<int64_t>(effectiveBand) + 1);
        const int64_t rightBand = static_cast<int64_t>(
                (diagonal + effectiveBand) / 2);
        start = std::max(start, leftBand);
        end = std::min(end, rightBand);
        if (start > end) {
            continue;
        }
        const uint64_t firstVector = static_cast<uint64_t>(start) /
                                     kKsw2SimdBytes;
        const uint64_t lastVector = static_cast<uint64_t>(end) /
                                    kKsw2SimdBytes;
        const uint64_t touchedBytes = saturatingMultiply(
                lastVector - firstVector + 1, kKsw2SimdBytes);
        const uint64_t rowOffset = saturatingMultiply(
                diagonal, rowStrideBytes);
        if (rowOffset == std::numeric_limits<uint64_t>::max() ||
            touchedBytes == std::numeric_limits<uint64_t>::max()) {
            return std::numeric_limits<uint64_t>::max() / pageSize;
        }
        const uint64_t firstPage = rowOffset / pageSize;
        const uint64_t lastTouchedByte = saturatingAdd(
                rowOffset, touchedBytes - 1);
        const uint64_t finalPage = lastTouchedByte / pageSize;
        if (!havePage || firstPage > lastPage) {
            pageCount = saturatingAdd(
                    pageCount, finalPage - firstPage + 1);
        } else if (finalPage > lastPage) {
            pageCount = saturatingAdd(pageCount, finalPage - lastPage);
        }
        lastPage = std::max(lastPage, finalPage);
        havePage = true;
    }
    // Large malloc-backed traceback regions are page aligned on supported
    // platforms. One additional page covers allocator/header displacement
    // without applying an empirical multiplicative discount.
    return saturatingAdd(pageCount, havePage ? 1 : 0);
}

uint64_t ksw2TracebackMemoryBytes(uint64_t queryLength,
                                  uint64_t referenceLength,
                                  uint64_t bandWidth) {
    if (queryLength == 0 || referenceLength == 0) {
        return 0;
    }
    if (queryLength > static_cast<uint64_t>(INT_MAX) ||
        referenceLength > static_cast<uint64_t>(INT_MAX) ||
        queryLength > std::numeric_limits<uint64_t>::max() -
                              referenceLength + 1) {
        return std::numeric_limits<uint64_t>::max();
    }
    const uint64_t pageSize = hostPageSizeBytes();
    const uint64_t targetBlocks = divideRoundUp(
            referenceLength, kKsw2SimdBytes);
    const uint64_t queryBlocks = divideRoundUp(
            queryLength, kKsw2SimdBytes);
    const uint64_t diagonalCount = queryLength + referenceLength - 1;
#if defined(__AVX512BW__) || defined(__AVX2__)
    const uint64_t primaryAllocation = saturatingMultiply(
            saturatingAdd(
                    saturatingAdd(saturatingMultiply(targetBlocks, 8),
                                  queryBlocks),
                    64),
            64);
#else
    const uint64_t primaryAllocation = saturatingMultiply(
            saturatingAdd(
                    saturatingAdd(saturatingMultiply(targetBlocks, 8),
                                  queryBlocks),
                    1),
            kKsw2SimdBytes);
#endif
    const uint64_t hAllocation = saturatingMultiply(
            saturatingMultiply(targetBlocks, kKsw2SimdBytes), 4);
    const uint64_t offsetsAllocation = saturatingMultiply(diagonalCount, 8);
    const uint64_t cigarAllocation = saturatingMultiply(diagonalCount, 4);
    const uint64_t encodedInputs = saturatingAdd(queryLength,
                                                 referenceLength);
    const uint64_t touchedTrace = saturatingMultiply(
            ksw2TouchedTracePages(queryLength, referenceLength, bandWidth,
                                  pageSize),
            pageSize);
    uint64_t resident = touchedTrace;
    for (const uint64_t allocation : {primaryAllocation, hAllocation,
                                      offsetsAllocation, cigarAllocation,
                                      encodedInputs}) {
        resident = saturatingAdd(
                resident, pageRoundedBytes(allocation, pageSize));
    }
    return resident;
}

uint64_t ksw2ScoreOnlyMemoryBytes(uint64_t queryLength,
                                  uint64_t referenceLength) {
    // extd2 score-only keeps SIMD working vectors and the encoded sequences
    // but omits p/off/off_end. AVX512 has the largest linear allocation among
    // the compiled kernels. Ninety-six bytes per input base plus a fixed
    // allocator/dispatch allowance is a portable conservative envelope.
    constexpr uint64_t kScoreOnlyBaseBytes = 32ULL * 1024ULL * 1024ULL;
    constexpr uint64_t kScoreOnlyBytesPerInputBase = 96;
    return saturatingAdd(
            kScoreOnlyBaseBytes,
            saturatingMultiply(
                    saturatingAdd(queryLength, referenceLength),
                    kScoreOnlyBytesPerInputBase));
}

uint64_t ksw2CertifiedReservationBytes(
        uint64_t scoreOnlyBytes,
        uint64_t queryLength,
        uint64_t referenceLength,
        uint64_t initialBand,
        uint64_t maximumBand) {
    if (initialBand > maximumBand) {
        return std::numeric_limits<uint64_t>::max();
    }
    // The unbanded score-only pass and every banded traceback attempt execute
    // sequentially. ksw_extd2_sse frees its score/trace workspace before the
    // next attempt starts, and failed attempt strings are destroyed at the end
    // of that iteration. Reserving the sum of every geometric expansion would
    // therefore price cumulative work as if it were simultaneous resident
    // memory and can incorrectly disable the rescue path. Traceback memory is
    // monotone in the band, so the largest attempt defines its peak; caller
    // transient accounting separately covers input/output strings.
    return std::max(
            scoreOnlyBytes,
            ksw2TracebackMemoryBytes(
                    queryLength, referenceLength, maximumBand));
}

uint64_t ksw2MaximumBandForBudget(uint64_t queryLength,
                                  uint64_t referenceLength,
                                  uint64_t scoreOnlyBytes,
                                  uint64_t initialBand,
                                  uint64_t budgetBytes) {
    const uint64_t maximumLength = std::max(queryLength, referenceLength);
    if (initialBand > maximumLength ||
        ksw2CertifiedReservationBytes(
                scoreOnlyBytes, queryLength, referenceLength,
                initialBand, initialBand) > budgetBytes) {
        return 0;
    }
    // Search the largest geometric-expansion cap whose peak attempt RSS fits
    // the raw algorithm budget. Attempts execute sequentially, so runtime work
    // is cumulative but memory is not. This still avoids a blanket empirical
    // discount and uses the architecture/page-aware resident model directly.
    uint64_t low = initialBand;
    uint64_t high = maximumLength;
    while (low < high) {
        const uint64_t middle = low + (high - low + 1) / 2;
        if (ksw2CertifiedReservationBytes(
                    scoreOnlyBytes, queryLength, referenceLength,
                    initialBand, middle) <= budgetBytes) {
            low = middle;
        } else {
            high = middle - 1;
        }
    }
    return low;
}

uint64_t ksw2MaximumSquareWindowForBudget(uint64_t requestedWindow,
                                           uint64_t budgetBytes) {
    if (requestedWindow == 0 || budgetBytes == 0) {
        return 0;
    }

    // Every interval in one genoAli invocation normally has the same -w.
    // Cache this relatively expensive, page-aware inversion instead of
    // rescanning up to O(w log w) antidiagonals for every gap profile.
    static std::mutex cacheMutex;
    static std::unordered_map<
            uint64_t, std::unordered_map<uint64_t, uint64_t>> cache;
    {
        std::lock_guard<std::mutex> lock(cacheMutex);
        const auto requested = cache.find(requestedWindow);
        if (requested != cache.end()) {
            const auto budget = requested->second.find(budgetBytes);
            if (budget != requested->second.end()) {
                return budget->second;
            }
        }
    }

    uint64_t low = 0;
    uint64_t high = requestedWindow;
    while (low < high) {
        const uint64_t middle = low + (high - low + 1) / 2;
        if (ksw2TracebackMemoryBytes(middle, middle, middle) <=
                budgetBytes) {
            low = middle;
        } else {
            high = middle - 1;
        }
    }
    {
        std::lock_guard<std::mutex> lock(cacheMutex);
        cache[requestedWindow].emplace(budgetBytes, low);
    }
    return low;
}

long double ksw2GeometricTraceWork(uint64_t maximumLength,
                                   uint64_t initialBand,
                                   uint64_t finalBand,
                                   long double scoringMultiplier) {
    if (initialBand > finalBand) {
        return 0.0L;
    }
    long double work = 0.0L;
    uint64_t band = initialBand;
    while (true) {
        work += 2.0L * static_cast<long double>(maximumLength) *
                static_cast<long double>(band) * scoringMultiplier;
        if (band == finalBand) {
            return work;
        }
        const uint64_t doubled = band >
                std::numeric_limits<uint64_t>::max() / 2
                ? finalBand : band * 2;
        const uint64_t nextBand = std::min(
                finalBand, std::max(band + 1, doubled));
        if (nextBand <= band) {
            return std::numeric_limits<long double>::infinity();
        }
        band = nextBand;
    }
}

uint64_t positivePenalty(int32_t penalty) {
    if (penalty == std::numeric_limits<int32_t>::min()) {
        return static_cast<uint64_t>(std::numeric_limits<int32_t>::max()) + 1;
    }
    return static_cast<uint64_t>(penalty < 0 ? -penalty : penalty);
}

uint64_t mixHash(uint64_t value) {
    value ^= value >> 30;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31;
    return value;
}

uint64_t sequenceFingerprint(const std::string &sequence) {
    // FNV-1a followed by the selector's avalanche mix. This is a provenance
    // guard, not a biological sketch: every byte and the length participate.
    uint64_t hash = 1469598103934665603ULL;
    for (const unsigned char value : sequence) {
        hash ^= static_cast<uint64_t>(value);
        hash *= 1099511628211ULL;
    }
    return mixHash(hash ^ mixHash(static_cast<uint64_t>(sequence.size())));
}

int encodedBase(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

SequenceSummary summarizeComposition(const std::string &sequence) {
    SequenceSummary summary;
    if (sequence.empty()) {
        return summary;
    }

    uint64_t ambiguous = 0;
    uint64_t lowComplexity = 0;
    char previous = '\0';
    uint64_t runLength = 0;
    for (char base : sequence) {
        const int code = encodedBase(base);
        if (code < 0) {
            ++ambiguous;
            previous = '\0';
            runLength = 0;
            continue;
        }
        const char normalized = static_cast<char>(std::toupper(
                static_cast<unsigned char>(base)));
        if (normalized == previous) {
            ++runLength;
        } else {
            previous = normalized;
            runLength = 1;
        }
        if (runLength >= 4) {
            ++lowComplexity;
        }
    }
    summary.ambiguousFraction = static_cast<double>(ambiguous) /
                                static_cast<double>(sequence.size());
    summary.lowComplexityFraction = static_cast<double>(lowComplexity) /
                                    static_cast<double>(sequence.size());
    return summary;
}

SequenceSummary summarizeSequence(const std::string &sequence,
                                  uint64_t samplingModulo) {
    SequenceSummary summary = summarizeComposition(sequence);

    if (sequence.size() < kKmerLength) {
        return summary;
    }

    const uint64_t mask = (1ULL << (2 * kKmerLength)) - 1;
    uint64_t forward = 0;
    uint64_t valid = 0;
    std::unordered_map<uint64_t, uint8_t> retainedOccurrences;
    retainedOccurrences.reserve(kTargetSketchSize * 2);
    summary.sketch.reserve(kTargetSketchSize);

    for (uint64_t position = 0; position < sequence.size(); ++position) {
        const int code = encodedBase(sequence[static_cast<std::size_t>(position)]);
        if (code < 0) {
            forward = 0;
            valid = 0;
            continue;
        }
        forward = ((forward << 2) | static_cast<uint64_t>(code)) & mask;
        ++valid;
        if (valid < kKmerLength) {
            continue;
        }
        // Gap sequences have already been oriented to the same strand by the
        // anchor chain. Keep the sketch strand-aware; canonical k-mers would
        // otherwise make, for example, A^n and T^n look identical merely
        // because they are reverse complements.
        const uint64_t hash = mixHash(forward);
        if (hash % samplingModulo != 0) {
            continue;
        }
        uint8_t &occurrences = retainedOccurrences[hash];
        if (occurrences >= kMaximumOccurrencesPerHash) {
            continue;
        }
        ++occurrences;
        const uint64_t start = position + 1 - kKmerLength;
        if (start <= std::numeric_limits<uint32_t>::max()) {
            summary.sketch.push_back(
                    SketchEntry{hash, static_cast<uint32_t>(start)});
        }
        if (summary.sketch.size() >= kMaximumSketchSize) {
            break;
        }
    }
    return summary;
}

using PositionMap = std::unordered_map<uint64_t, std::vector<uint32_t>>;

PositionMap indexSketch(const std::vector<SketchEntry> &sketch) {
    PositionMap index;
    index.reserve(sketch.size());
    for (const SketchEntry &entry : sketch) {
        index[entry.hash].push_back(entry.position);
    }
    return index;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    const std::size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() % 2 != 0) {
        return upper;
    }
    std::nth_element(values.begin(), values.begin() + middle - 1,
                     values.begin() + middle);
    return (upper + values[middle - 1]) / 2.0;
}

double quantile(std::vector<double> values, double probability) {
    if (values.empty()) {
        return 0.0;
    }
    probability = std::max(0.0, std::min(1.0, probability));
    const std::size_t index = static_cast<std::size_t>(std::floor(
            probability * static_cast<double>(values.size() - 1)));
    std::nth_element(values.begin(), values.begin() + index, values.end());
    return values[index];
}

std::vector<Anchor> longestIncreasingChain(std::vector<Anchor> anchors) {
    if (anchors.empty()) {
        return {};
    }
    std::sort(anchors.begin(), anchors.end(),
              [](const Anchor &first, const Anchor &second) {
                  if (first.query != second.query) {
                      return first.query < second.query;
                  }
                  return first.reference > second.reference;
              });

    std::vector<uint32_t> tails;
    std::vector<std::size_t> tailIndices;
    std::vector<std::size_t> predecessors(
            anchors.size(), std::numeric_limits<std::size_t>::max());
    tails.reserve(anchors.size());
    tailIndices.reserve(anchors.size());
    for (std::size_t i = 0; i < anchors.size(); ++i) {
        const auto found = std::lower_bound(
                tails.begin(), tails.end(), anchors[i].reference);
        const std::size_t length = static_cast<std::size_t>(
                found - tails.begin());
        if (found == tails.end()) {
            tails.push_back(anchors[i].reference);
            tailIndices.push_back(i);
        } else {
            *found = anchors[i].reference;
            tailIndices[length] = i;
        }
        if (length > 0) {
            predecessors[i] = tailIndices[length - 1];
        }
    }

    std::vector<Anchor> chain;
    chain.reserve(tails.size());
    std::size_t index = tailIndices.back();
    while (index != std::numeric_limits<std::size_t>::max()) {
        chain.push_back(anchors[index]);
        index = predecessors[index];
    }
    std::reverse(chain.begin(), chain.end());
    return chain;
}

uint64_t gapCost(uint64_t length,
                 uint64_t open1,
                 uint64_t extend1,
                 uint64_t open2,
                 uint64_t extend2) {
    if (length == 0) {
        return 0;
    }
    const uint64_t first = saturatingAdd(open1,
                                         saturatingMultiply(extend1, length));
    const uint64_t second = saturatingAdd(open2,
                                          saturatingMultiply(extend2, length));
    return std::min(first, second);
}

uint64_t roundedSaturating(long double value) {
    if (!(value > 0.0L)) {
        return 0;
    }
    if (value >= static_cast<long double>(
            std::numeric_limits<uint64_t>::max())) {
        return std::numeric_limits<uint64_t>::max();
    }
    return static_cast<uint64_t>(std::ceil(value));
}

AlgorithmCostEstimate estimate(AlignmentCandidate candidate,
                               uint64_t memoryBytes,
                               long double workUnits,
                               bool memoryFeasible,
                               bool workFeasible,
                               double estimatedMinutes = 0.0,
                               double estimatedMinutesP90 = 0.0,
                               double successProbability = 1.0,
                               bool timeFeasible = true,
                               bool withinConfiguredTimeLimit = true) {
    AlgorithmCostEstimate result;
    result.candidate = candidate;
    result.memoryBytes = memoryBytes;
    result.workUnits = workUnits;
    result.estimatedMinutes = estimatedMinutes;
    result.estimatedMinutesP90 = std::max(estimatedMinutes,
                                         estimatedMinutesP90);
    result.successProbability = std::max(
            0.0, std::min(1.0, successProbability));
    result.memoryFeasible = memoryFeasible;
    result.workFeasible = workFeasible;
    result.withinConfiguredTimeLimit = withinConfiguredTimeLimit;
    result.timeFeasible = timeFeasible;
    result.feasible = memoryFeasible && timeFeasible;
    return result;
}

double finiteMinutes(long double workUnits, long double unitsPerSecond,
                     double setupSeconds) {
    const long double minutes =
            (workUnits / unitsPerSecond + setupSeconds) / 60.0L;
    if (!(minutes > 0.0L)) {
        return 0.0;
    }
    if (minutes >= static_cast<long double>(
                           std::numeric_limits<double>::max())) {
        return std::numeric_limits<double>::max();
    }
    return static_cast<double>(minutes);
}

double wfaMinutes(uint64_t maximumLength, long double score,
                  long double multiplier) {
    const long double seconds = multiplier *
            (static_cast<long double>(maximumLength) /
                     kWfaLinearBasesPerSecond +
             score * score / kWfaScoreSquaredPerSecond) +
            kWfaSetupSeconds;
    if (!(seconds > 0.0L)) {
        return 0.0;
    }
    const long double minutes = seconds / 60.0L;
    return minutes >= static_cast<long double>(
                              std::numeric_limits<double>::max())
           ? std::numeric_limits<double>::max()
           : static_cast<double>(minutes);
}

double wfaMinutesAtScoreThroughput(uint64_t maximumLength,
                                   long double score,
                                   long double scoreSquaredPerSecond) {
    const long double seconds =
            static_cast<long double>(maximumLength) /
                    kWfaLinearBasesPerSecond +
            score * score / scoreSquaredPerSecond + kWfaSetupSeconds;
    if (!(seconds > 0.0L)) {
        return 0.0;
    }
    const long double minutes = seconds / 60.0L;
    return minutes >= static_cast<long double>(
                              std::numeric_limits<double>::max())
           ? std::numeric_limits<double>::max()
           : static_cast<double>(minutes);
}

long double singletrackP50Multiplier(const AlignmentProfile &profile) {
    long double multiplier = kSingletrackP50LargeScoreMultiplier;
    if (profile.estimatedScore == 0) {
        multiplier = kSingletrackP50ScoreZeroMultiplier;
    } else if (profile.estimatedScore < 1000) {
        multiplier = kSingletrackP50SmallScoreMultiplier;
    } else if (profile.estimatedScore < 10000) {
        multiplier = kSingletrackP50MediumScoreMultiplier;
    }
    multiplier += kSingletrackP50UncertaintyMultiplier *
                  static_cast<long double>(profile.uncertainty);
    return std::min(kSingletrackP50MaximumMultiplier, multiplier);
}

double p90Inflation(const AlignmentProfile &profile,
                    AlignmentCandidate candidate) {
    double inflation = 1.12 + 0.55 * profile.uncertainty;
    if (candidate == AlignmentCandidate::MediumWfa ||
               candidate == AlignmentCandidate::LowWfa) {
        inflation += 0.08;
    } else if (candidate == AlignmentCandidate::Ksw2ScoreCertified ||
               candidate == AlignmentCandidate::Ksw2Full ||
               candidate == AlignmentCandidate::Ksw2Banded) {
        inflation = 1.08 + 0.20 * profile.uncertainty;
    }
    return inflation;
}

}  // namespace

uint64_t ksw2TracebackResidentMemoryBytes(
        uint64_t queryLength,
        uint64_t referenceLength,
        uint64_t bandWidth) {
    return ksw2TracebackMemoryBytes(
            queryLength, referenceLength, bandWidth);
}

const char *alignmentCandidateName(AlignmentCandidate candidate) {
    switch (candidate) {
        case AlignmentCandidate::SingletrackWfa:
            return "WAVEFRONT_SINGLETRACK";
        case AlignmentCandidate::StandardWfa: return "WAVEFRONT";
        case AlignmentCandidate::MediumWfa: return "WAVEFRONT_MEDIUM";
        case AlignmentCandidate::LowWfa: return "WAVEFRONT_LOW";
        case AlignmentCandidate::Ksw2ScoreCertified:
            return "KSW2_SCORE_CERTIFIED";
        case AlignmentCandidate::Ksw2Full: return "MINIMAP2";
        case AlignmentCandidate::Ksw2Banded: return "BANDED_MINIMAP2";
        case AlignmentCandidate::SlidingWindow: return "SLIDING_WINDOW";
    }
    return "UNKNOWN";
}

std::string alignmentMethodBedLabel(const std::string &internalMethod) {
    if (internalMethod == "WAVEFRONT_SINGLETRACK" ||
        internalMethod == "WAVEFRONT" ||
        internalMethod == "WAVEFRONT_MEDIUM" ||
        internalMethod == "WAVEFRONT_LOW") {
        return "WAVEFRONT";
    }
    if (internalMethod == "KSW2_SCORE_CERTIFIED") {
        return "MINIMAP2";
    }
    return internalMethod;
}

double alignmentRiskAdjustedMinutes(
        const AlgorithmCostEstimate &estimate) {
    const double p50 = std::max(0.0, estimate.estimatedMinutes);
    const double p90 = std::max(p50, estimate.estimatedMinutesP90);
    return p50 + 0.35 * (p90 - p50) +
           (1.0 - std::max(0.0, std::min(
                          1.0, estimate.successProbability))) * p90;
}

double alignmentSelectionPriorityMinutes(
        const AlignmentSelectionPlan &plan) {
    const std::vector<AlignmentCandidate> *tier = nullptr;
    if (!plan.exactCandidates.empty()) {
        tier = &plan.exactCandidates;
    } else if (!plan.approximateCandidates.empty()) {
        tier = &plan.approximateCandidates;
    } else if (!plan.bandedCandidates.empty()) {
        tier = &plan.bandedCandidates;
    } else {
        tier = &plan.lastResortCandidates;
    }
    double priority = std::numeric_limits<double>::infinity();
    for (const AlignmentCandidate candidate : *tier) {
        const auto found = std::find_if(
                plan.estimates.begin(), plan.estimates.end(),
                [candidate](const AlgorithmCostEstimate &estimate) {
                    return estimate.candidate == candidate;
                });
        if (found != plan.estimates.end()) {
            priority = std::min(
                    priority, alignmentRiskAdjustedMinutes(*found));
        }
    }
    return std::isfinite(priority) ? std::max(0.0, priority) : 0.0;
}

uint64_t alignmentSelectionPriorityMemoryBytes(
        const AlignmentSelectionPlan &plan) {
    const std::vector<AlignmentCandidate> *tier = nullptr;
    if (!plan.exactCandidates.empty()) {
        tier = &plan.exactCandidates;
    } else if (!plan.approximateCandidates.empty()) {
        tier = &plan.approximateCandidates;
    } else if (!plan.bandedCandidates.empty()) {
        tier = &plan.bandedCandidates;
    } else {
        tier = &plan.lastResortCandidates;
    }
    const AlgorithmCostEstimate *best = nullptr;
    for (const AlignmentCandidate candidate : *tier) {
        const auto found = std::find_if(
                plan.estimates.begin(), plan.estimates.end(),
                [candidate](const AlgorithmCostEstimate &estimate) {
                    return estimate.candidate == candidate;
                });
        if (found == plan.estimates.end()) {
            continue;
        }
        if (best == nullptr ||
            alignmentRiskAdjustedMinutes(*found) <
                    alignmentRiskAdjustedMinutes(*best)) {
            best = &*found;
        }
    }
    return best != nullptr ? best->memoryBytes : 0;
}

bool exactCandidateWithinTimeLimit(
        double estimatedMinutesP50,
        double estimatedMinutesP90,
        double maximumEstimatedMinutes) {
    if (!std::isfinite(maximumEstimatedMinutes) ||
        maximumEstimatedMinutes < 0.0) {
        throw std::invalid_argument(
                "exact-alignment maximum estimated minutes must be finite and >= 0");
    }
    if (maximumEstimatedMinutes == 0.0) {
        return true;
    }
    const double p50 = std::max(0.0, estimatedMinutesP50);
    const double p90 = std::max(p50, estimatedMinutesP90);
    return p50 <= maximumEstimatedMinutes &&
           p90 <= maximumEstimatedMinutes;
}

bool highWfaFastLaneEligible(
        AlignmentCandidate candidate,
        const AlgorithmCostEstimate &estimate,
        double predictedWaitMinutes,
        double fastestExactRiskAdjustedMinutes) {
    if (candidate != AlignmentCandidate::SingletrackWfa &&
        candidate != AlignmentCandidate::StandardWfa) {
        return false;
    }
    if (predictedWaitMinutes > 0.0 ||
        !std::isfinite(predictedWaitMinutes) ||
        !std::isfinite(fastestExactRiskAdjustedMinutes)) {
        return false;
    }
    const double highRuntime = alignmentRiskAdjustedMinutes(estimate);
    return highRuntime <= fastestExactRiskAdjustedMinutes;
}

bool fastExactDominatesSlowExact(
        AlignmentCandidate fastCandidate,
        const AlgorithmCostEstimate &fastEstimate,
        double fastWaitMinutes,
        AlignmentCandidate slowCandidate,
        const AlgorithmCostEstimate &slowEstimate,
        double fastMemoryShadowMinutes,
        double slowMemoryShadowMinutes,
        double maximumWaitMinutes) {
    const bool fast = fastCandidate == AlignmentCandidate::SingletrackWfa ||
                      fastCandidate == AlignmentCandidate::StandardWfa ||
                      fastCandidate ==
                              AlignmentCandidate::Ksw2ScoreCertified ||
                      fastCandidate == AlignmentCandidate::Ksw2Full;
    const bool slow = slowCandidate == AlignmentCandidate::MediumWfa ||
                      slowCandidate == AlignmentCandidate::LowWfa;
    if (!fast || !slow || !std::isfinite(fastWaitMinutes) ||
        !std::isfinite(maximumWaitMinutes) || fastWaitMinutes < 0.0 ||
        maximumWaitMinutes < 0.0 ||
        !std::isfinite(fastMemoryShadowMinutes) ||
        !std::isfinite(slowMemoryShadowMinutes)) {
        return false;
    }
    const double fastP50 = std::max(0.0, fastEstimate.estimatedMinutes);
    const double fastP90 = std::max(
            fastP50, fastEstimate.estimatedMinutesP90);
    const double slowP50 = std::max(0.0, slowEstimate.estimatedMinutes);
    const double fastShadow = std::max(0.0, fastMemoryShadowMinutes);
    const double slowShadow = std::max(0.0, slowMemoryShadowMinutes);
    if (!std::isfinite(fastP90) || !std::isfinite(slowP50)) {
        return false;
    }

    // One second covers a park/requeue/admission transition.  For longer
    // exact tasks use ten percent of fast P90 so a statistically insignificant
    // difference cannot drain the process for a large reservation.
    const double guardMinutes = std::max(1.0 / 60.0,
                                         0.10 * fastP90);
    const double slowCompletion = slowP50 + slowShadow;
    const double fastWithoutWait = fastP90 + fastShadow + guardMinutes;
    const double positiveSlack = std::max(
            0.0, slowCompletion - fastWithoutWait);
    const double holdBudget = std::min(maximumWaitMinutes, positiveSlack);
    const double fastCompletion = fastWaitMinutes + fastWithoutWait;
    return fastWaitMinutes <= holdBudget &&
           fastCompletion < slowCompletion;
}

void configureExactAlignmentMaximumEstimatedMinutes(double minutes) {
    if (!std::isfinite(minutes) || minutes < 0.0) {
        throw std::invalid_argument(
                "exact-alignment maximum estimated minutes must be finite and >= 0");
    }
    exactAlignmentMaximumEstimatedMinutes.store(minutes,
                                                 std::memory_order_relaxed);
}

double configuredExactAlignmentMaximumEstimatedMinutes() {
    return exactAlignmentMaximumEstimatedMinutes.load(
            std::memory_order_relaxed);
}

AlignmentSelectionPlan makeAlignmentSelectionPlan(
        const std::string &query,
        const std::string &reference,
        int64_t windowWidth,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        double exactAlignmentMaximumMinutes,
        uint64_t elasticHighWfaMemoryBudgetBytes,
        uint64_t elasticFullKsw2MemoryBudgetBytes,
        bool memorySchedulingEnabled) {
    if (!std::isfinite(exactAlignmentMaximumMinutes) ||
        exactAlignmentMaximumMinutes < 0.0) {
        throw std::invalid_argument(
                "exact-alignment maximum estimated minutes must be finite and >= 0");
    }
    AlignmentSelectionPlan plan;
    plan.provenance.queryFingerprint = sequenceFingerprint(query);
    plan.provenance.referenceFingerprint = sequenceFingerprint(reference);
    plan.provenance.windowWidth = windowWidth;
    plan.provenance.mismatchingPenalty = mismatchingPenalty;
    plan.provenance.openGapPenalty1 = openGapPenalty1;
    plan.provenance.extendGapPenalty1 = extendGapPenalty1;
    plan.provenance.openGapPenalty2 = openGapPenalty2;
    plan.provenance.extendGapPenalty2 = extendGapPenalty2;
    plan.provenance.exactAlignmentMaximumEstimatedMinutes =
            exactAlignmentMaximumMinutes;
    plan.provenance.elasticHighWfaMemoryBudgetBytes =
            elasticHighWfaMemoryBudgetBytes;
    plan.provenance.elasticFullKsw2MemoryBudgetBytes =
            elasticFullKsw2MemoryBudgetBytes;
    plan.provenance.memorySchedulingEnabled = memorySchedulingEnabled;
    plan.provenance.valid = true;
    AlignmentProfile &profile = plan.profile;
    profile.queryLength = query.size();
    profile.referenceLength = reference.size();
    profile.identical = query == reference;

    const uint64_t maximumLength = std::max(profile.queryLength,
                                            profile.referenceLength);
    const uint64_t minimumLength = std::min(profile.queryLength,
                                            profile.referenceLength);
    const uint64_t availableKmers = maximumLength >= kKmerLength
                                    ? maximumLength - kKmerLength + 1
                                    : 0;
    const uint64_t samplingModulo = std::max<uint64_t>(
            1, availableKmers / kTargetSketchSize);
    // query == reference already performs one linear comparison.  Avoid two
    // more full scans plus hash-table construction for the common exact case.
    SequenceSummary querySummary;
    SequenceSummary referenceSummary;
    if (!profile.identical) {
        querySummary = summarizeSequence(query, samplingModulo);
        referenceSummary = summarizeSequence(reference, samplingModulo);
    } else {
        // The equality fast path skips the comparatively expensive sketches,
        // but composition still matters: identical mixed A/C/G/T/N input
        // makes KSW2 use its generic scoring kernel and must be priced as such.
        querySummary = summarizeComposition(query);
        referenceSummary = querySummary;
    }
    profile.sampledQueryKmers = querySummary.sketch.size();
    profile.sampledReferenceKmers = referenceSummary.sketch.size();
    profile.ambiguousBaseFraction =
            (querySummary.ambiguousFraction * profile.queryLength +
             referenceSummary.ambiguousFraction * profile.referenceLength) /
            static_cast<double>(std::max<uint64_t>(1,
                    profile.queryLength + profile.referenceLength));
    profile.lowComplexityFraction =
            (querySummary.lowComplexityFraction * profile.queryLength +
             referenceSummary.lowComplexityFraction * profile.referenceLength) /
            static_cast<double>(std::max<uint64_t>(1,
                    profile.queryLength + profile.referenceLength));

    const PositionMap queryIndex = indexSketch(querySummary.sketch);
    const PositionMap referenceIndex = indexSketch(referenceSummary.sketch);
    uint64_t sharedHashes = 0;
    uint64_t uniqueHashes = 0;
    std::vector<Anchor> anchors;
    std::vector<Anchor> uniqueAnchors;
    anchors.reserve(std::min(querySummary.sketch.size(),
                             referenceSummary.sketch.size()) * 2);
    for (const auto &item : queryIndex) {
        const auto found = referenceIndex.find(item.first);
        if (found == referenceIndex.end()) {
            continue;
        }
        ++sharedHashes;
        if (item.second.size() == 1 && found->second.size() == 1) {
            ++uniqueHashes;
            uniqueAnchors.push_back(
                    Anchor{item.second.front(), found->second.front()});
        }
        for (uint32_t queryPosition : item.second) {
            for (uint32_t referencePosition : found->second) {
                anchors.push_back(Anchor{queryPosition, referencePosition});
            }
        }
    }
    const uint64_t unionHashes = queryIndex.size() + referenceIndex.size() -
                                 sharedHashes;
    profile.sketchJaccard = unionHashes == 0
                            ? (profile.identical ? 1.0 : 0.0)
                            : static_cast<double>(sharedHashes) /
                              static_cast<double>(unionHashes);
    profile.uniqueAnchorFraction = sharedHashes == 0
                                   ? 0.0
                                   : static_cast<double>(uniqueHashes) /
                                     static_cast<double>(sharedHashes);

    const std::vector<Anchor> chain = longestIncreasingChain(std::move(anchors));
    const std::vector<Anchor> uniqueChain = longestIncreasingChain(
            std::move(uniqueAnchors));
    uint64_t previousQueryEnd = 0;
    uint64_t previousReferenceEnd = 0;
    bool haveCertificationAnchor = false;
    for (const Anchor &anchor : uniqueChain) {
        const uint64_t queryEnd = static_cast<uint64_t>(anchor.query) +
                                  kKmerLength;
        const uint64_t referenceEnd =
                static_cast<uint64_t>(anchor.reference) + kKmerLength;
        if (queryEnd > query.size() || referenceEnd > reference.size()) {
            continue;
        }
        bool exactMatch = true;
        for (uint64_t offset = 0; offset < kKmerLength; ++offset) {
            if (encodedBase(query[anchor.query + offset]) !=
                encodedBase(reference[anchor.reference + offset])) {
                exactMatch = false;
                break;
            }
        }
        if (!exactMatch ||
            (haveCertificationAnchor &&
             (anchor.query < previousQueryEnd ||
              anchor.reference < previousReferenceEnd))) {
            continue;
        }
        if (haveCertificationAnchor &&
            anchor.query - previousQueryEnd < kCertificationAnchorSpacing &&
            anchor.reference - previousReferenceEnd <
                    kCertificationAnchorSpacing) {
            continue;
        }
        plan.certificationAnchors.push_back(AlignmentChainAnchor{
                anchor.query, anchor.reference,
                static_cast<uint16_t>(kKmerLength)});
        previousQueryEnd = queryEnd;
        previousReferenceEnd = referenceEnd;
        haveCertificationAnchor = true;
    }
    profile.chainedAnchors = chain.size();
    const uint64_t sketchDenominator = std::min(querySummary.sketch.size(),
                                                referenceSummary.sketch.size());
    profile.chainedAnchorFraction = sketchDenominator == 0
                                    ? 0.0
                                    : std::min(1.0,
                                      static_cast<double>(chain.size()) /
                                      static_cast<double>(sketchDenominator));

    profile.lengthRatio = maximumLength == 0 ? 1.0 :
            static_cast<double>(minimumLength) /
            static_cast<double>(maximumLength);
    if (!chain.empty()) {
        std::vector<double> diagonals;
        diagonals.reserve(chain.size());
        double previousDiagonal = 0.0;
        bool havePrevious = false;
        for (const Anchor &anchor : chain) {
            const double diagonal = static_cast<double>(anchor.reference) -
                                    static_cast<double>(anchor.query);
            diagonals.push_back(diagonal);
            if (havePrevious) {
                profile.maximumDiagonalJump = std::max(
                        profile.maximumDiagonalJump,
                        std::fabs(diagonal - previousDiagonal));
            }
            previousDiagonal = diagonal;
            havePrevious = true;
        }
        const double diagonalMedian = median(diagonals);
        for (double &diagonal : diagonals) {
            diagonal = std::fabs(diagonal - diagonalMedian);
        }
        profile.diagonalMad = median(diagonals);
        profile.diagonalP90 = quantile(diagonals, 0.90);
        profile.diagonalP99 = quantile(std::move(diagonals), 0.99);
    }

    if (profile.identical) {
        profile.estimatedMismatchRate = 0.0;
        profile.confidence = 1.0;
    } else if (profile.sketchJaccard > 0.0) {
        const double sharedKmerProbability =
                2.0 * profile.sketchJaccard /
                (1.0 + profile.sketchJaccard);
        profile.estimatedMismatchRate = std::max(
                0.0, std::min(1.0, 1.0 - std::pow(
                        sharedKmerProbability,
                        1.0 / static_cast<double>(kKmerLength))));
        profile.estimatedMismatchRate = std::max(
                profile.estimatedMismatchRate,
                (1.0 - profile.chainedAnchorFraction) * 0.10);
    } else {
        profile.estimatedMismatchRate = minimumLength == 0 ? 0.0 : 0.75;
    }

    const double sketchSupport = std::min(
            1.0, static_cast<double>(sketchDenominator) / 256.0);
    const double repeatConfidence = 0.25 +
                                    0.75 * profile.uniqueAnchorFraction;
    const double chainConfidence = std::min(
            1.0, profile.chainedAnchorFraction * 2.0);
    profile.confidence = profile.identical ? 1.0 :
            sketchSupport * repeatConfidence * chainConfidence *
            (1.0 - 0.5 * profile.ambiguousBaseFraction);
    profile.confidence = std::max(0.0, std::min(1.0, profile.confidence));
    const double normalizedDiagonalRisk = maximumLength == 0 ? 0.0 :
            std::min(1.0,
                     (profile.diagonalP99 * 2.0 +
                      profile.maximumDiagonalJump) /
                     static_cast<double>(maximumLength));
    profile.uncertainty = profile.identical ? 0.0 : std::min(
            1.0,
            0.62 * (1.0 - profile.confidence) +
            0.18 * profile.lowComplexityFraction +
            0.20 * profile.ambiguousBaseFraction +
            0.35 * normalizedDiagonalRisk);

    const uint64_t mismatch = positivePenalty(mismatchingPenalty);
    const uint64_t gapOpen1 = positivePenalty(openGapPenalty1);
    const uint64_t gapExtend1 = positivePenalty(extendGapPenalty1);
    const uint64_t gapOpen2 = positivePenalty(openGapPenalty2);
    const uint64_t gapExtend2 = positivePenalty(extendGapPenalty2);
    const uint64_t lengthDifference = maximumLength - minimumLength;
    const uint64_t unavoidableGapScore = gapCost(
            lengthDifference, gapOpen1, gapExtend1,
            gapOpen2, gapExtend2);
    const long double mismatchScore =
            static_cast<long double>(minimumLength) *
            profile.estimatedMismatchRate * mismatch;
    profile.estimatedScore = saturatingAdd(
            unavoidableGapScore, roundedSaturating(mismatchScore));
    // This is a memory-tail score, not a time estimate.  It deliberately
    // combines sketch confidence, repeat content, and chained-diagonal jumps;
    // balanced long indels are otherwise invisible to net length difference.
    const long double structuralBoost = std::min<long double>(
            1.50L, 4.0L * normalizedDiagonalRisk);
    const long double weakSketchBoost = std::max<long double>(
            0.0L, (0.45L - profile.sketchJaccard) * 2.50L);
    const long double uncertaintyMultiplier = std::min<long double>(
            3.50L,
            1.15L + static_cast<long double>(profile.uncertainty) * 0.85L +
            structuralBoost + weakSketchBoost);
    profile.conservativeScore = std::max(
            unavoidableGapScore,
            roundedSaturating(profile.estimatedScore *
                              uncertaintyMultiplier));
    if (profile.identical) {
        profile.estimatedScore = 0;
        profile.conservativeScore = 0;
    }
    profile.estimatedBandWidth = std::max<uint64_t>(
            lengthDifference,
            roundedSaturating(
                    2.0L * profile.diagonalP99 +
                    profile.maximumDiagonalJump + 64.0L));

    const uint64_t workerBudget = wfaMemoryBudgetBytes(windowWidth);
    const uint64_t standardBudget =
            standardWfaTrialMemoryBudgetBytes(workerBudget);
    const uint64_t scoreSquared = saturatingMultiply(
            profile.conservativeScore, profile.conservativeScore);
    const uint64_t singletrackDiagonalSpan = saturatingAdd(
            saturatingAdd(profile.queryLength, profile.referenceLength),
            kSingletrackDiagonalPadding);
    // Each I/D matrix keeps at most 2*gap_extension wavefronts alive. There
    // are two matrices for each of the two affine pieces.
    const uint64_t singletrackRecycledWavefronts = saturatingMultiply(
            4, saturatingAdd(gapExtend1, gapExtend2));
    const uint64_t singletrackIndelPool = saturatingMultiply(
            saturatingMultiply(singletrackDiagonalSpan,
                               singletrackRecycledWavefronts),
            kSingletrackOffsetBytes);
    const uint64_t singletrackMemory = saturatingAdd(
            saturatingAdd(
                    kSingletrackWfaBaseBytes,
                    saturatingMultiply(
                            kSingletrackWfaBytesPerScoreSquared,
                            scoreSquared)),
            singletrackIndelPool);
    const uint64_t standardMemory = saturatingAdd(
            kStandardWfaBaseBytes,
            saturatingMultiply(kStandardWfaBytesPerScoreSquared,
                               scoreSquared));
    const uint64_t mediumWfaMemory = saturatingAdd(
            kSuccinctWfaBaseBytes,
            saturatingMultiply(kMediumWfaBytesPerScoreSquared,
                               scoreSquared));
    const uint64_t lowWfaMemory = saturatingAdd(
            kSuccinctWfaBaseBytes,
            saturatingMultiply(kLowWfaBytesPerScoreSquared,
                               scoreSquared));
    // Time and memory intentionally use different score quantiles.  The old
    // model fed the memory-tail score into time prediction, coupling safe
    // admission to slow and unstable engine ordering.  B73/Mo17 calibration
    // supports 1.5*estimatedScore as P50; P90 widens with uncertainty.
    const long double estimatedScore = static_cast<long double>(
            profile.estimatedScore);
    const long double timeScoreP50 = estimatedScore * 1.50L;
    const long double timeScoreP90 = estimatedScore *
            (1.80L + 0.80L * profile.uncertainty);
    const long double wfaCoreWork =
            static_cast<long double>(maximumLength) * timeScoreP50 +
            timeScoreP50 * timeScoreP50;
    const long double singletrackP50Factor =
            singletrackP50Multiplier(profile);
    const long double standardWfaWork =
            wfaCoreWork * kStandardWfaWorkMultiplier;
    const long double singletrackWfaWork =
            wfaCoreWork * singletrackP50Factor;
    const long double mediumWfaWork =
            wfaCoreWork * kMediumWfaWorkMultiplier;
    const long double lowWfaWork =
            wfaCoreWork * kLowWfaWorkMultiplier;
    const double singletrackWfaEstimatedMinutes = wfaMinutes(
            maximumLength, timeScoreP50,
            singletrackP50Factor);
    const long double singletrackP90Score = std::max(
            kSingletrackP90EstimatedScoreMultiplier * estimatedScore,
            kSingletrackP90ConservativeScoreMultiplier *
                    static_cast<long double>(profile.conservativeScore));
    const double singletrackWfaEstimatedMinutesP90 =
            wfaMinutesAtScoreThroughput(
                    maximumLength, singletrackP90Score,
                    kSingletrackP90ScoreSquaredPerSecond);
    const double standardWfaEstimatedMinutes = wfaMinutes(
            maximumLength, timeScoreP50, kStandardWfaWorkMultiplier);
    const double standardWfaEstimatedMinutesP90 = wfaMinutes(
            maximumLength, timeScoreP90, kStandardWfaWorkMultiplier);
    const double mediumWfaEstimatedMinutes = wfaMinutes(
            maximumLength, timeScoreP50, kMediumWfaWorkMultiplier);
    const double mediumWfaEstimatedMinutesP90 = wfaMinutes(
            maximumLength, timeScoreP90, kMediumWfaWorkMultiplier);
    const double lowWfaEstimatedMinutes = wfaMinutes(
            maximumLength, timeScoreP50, kLowWfaWorkMultiplier);
    const double lowWfaEstimatedMinutesP90 = wfaMinutes(
            maximumLength, timeScoreP90, kLowWfaWorkMultiplier);
    const uint64_t requestedSlidingWidth = windowWidth > 0
            ? static_cast<uint64_t>(windowWidth) : 0;
    const uint64_t memorySafeSlidingWidth =
            ksw2MaximumSquareWindowForBudget(
                    requestedSlidingWidth, workerBudget);
    // Keep the widest sliding window permitted by the original w^2 worker
    // budget. Reducing this width is a biological-quality change, so resource
    // scheduling must never use it as a throughput control.
    const uint64_t slidingWidth = std::min<uint64_t>(
            maximumLength, memorySafeSlidingWidth);
    const uint64_t slidingMemory = ksw2TracebackMemoryBytes(
            std::min(profile.queryLength, slidingWidth),
            std::min(profile.referenceLength, slidingWidth),
            slidingWidth);
    // A sliding pass stops as soon as either sequence is exhausted; any
    // remaining tail is emitted as one gap and does not create another DP
    // matrix.  Model the actual sequence of square chunks instead of the old
    // max(length)*window rectangle, which substantially overestimated highly
    // unbalanced maize gaps while underestimating kernel time via its unrelated
    // throughput coefficient.
    long double slidingWork = 0.0L;
    if (slidingWidth > 0) {
        const uint64_t completeChunks = minimumLength / slidingWidth;
        slidingWork = static_cast<long double>(completeChunks) *
                      slidingWidth * slidingWidth;
        const uint64_t consumed = saturatingMultiply(completeChunks,
                                                     slidingWidth);
        const uint64_t remainingMinimum = minimumLength > consumed
                ? minimumLength - consumed : 0;
        const uint64_t remainingMaximum = maximumLength > consumed
                ? maximumLength - consumed : 0;
        slidingWork += static_cast<long double>(remainingMinimum) *
                       std::min<uint64_t>(slidingWidth, remainingMaximum);
    }
    const bool lengthsFitInt = profile.queryLength <= INT_MAX &&
                               profile.referenceLength <= INT_MAX;
    const auto workWithinSlidingEnvelope = [slidingWork](long double work) {
        return work == 0.0L ||
               (slidingWork > 0.0L && work <= slidingWork);
    };
    // -w is a hard per-alignment algorithm-memory ceiling.  A process-wide
    // -M limit controls concurrency, but must never make an individual WFA or
    // KSW2 attempt larger than w^2.  Keep the legacy elastic parameters as
    // optional lower ceilings for source compatibility; they can no longer
    // raise either high-WFA mode above the -w-derived budget.
    const uint64_t singletrackBudget = elasticHighWfaMemoryBudgetBytes > 0
                                       ? std::min(workerBudget,
                                                  elasticHighWfaMemoryBudgetBytes)
                                       : workerBudget;
    const uint64_t highStandardBudget = elasticHighWfaMemoryBudgetBytes > 0
                                        ? std::min(standardBudget,
                                                   elasticHighWfaMemoryBudgetBytes)
                                        : standardBudget;
    const bool singletrackMemoryFeasible =
            lengthsFitInt && singletrackBudget > 0 &&
            singletrackMemory <= singletrackBudget;
    const bool standardMemoryFeasible = lengthsFitInt &&
                                        highStandardBudget > 0 &&
                                        standardMemory <= highStandardBudget;
    const bool mediumWfaMemoryFeasible =
            lengthsFitInt && workerBudget > 0 &&
            mediumWfaMemory <= workerBudget;
    const bool lowWfaMemoryFeasible =
            lengthsFitInt && workerBudget > 0 &&
            lowWfaMemory <= workerBudget;
    const uint64_t scoreOnlyKswMemory = ksw2ScoreOnlyMemoryBytes(
            profile.queryLength, profile.referenceLength);
    plan.ksw2CertifiedScoreOnlyMemoryBytes = scoreOnlyKswMemory;
    const uint64_t fullKswMemory = ksw2TracebackMemoryBytes(
            profile.queryLength, profile.referenceLength,
            std::max(profile.queryLength, profile.referenceLength));
    const uint64_t fullKswBudget = elasticFullKsw2MemoryBudgetBytes > 0
                                   ? std::min(workerBudget,
                                              elasticFullKsw2MemoryBudgetBytes)
                                   : workerBudget;
    const bool fullKswFeasible = lengthsFitInt && fullKswBudget > 0 &&
                                 fullKswMemory <= fullKswBudget;
    const bool fastExactStructurallyUnavailable =
            !singletrackMemoryFeasible &&
            !standardMemoryFeasible &&
            !fullKswFeasible;

    // Score-certified KSW2 separates exact scoring from traceback. Start from
    // a chain-informed band, then geometrically double inside one guarded
    // reservation.  Like every other sequence-alignment engine, the
    // score-only plus traceback peak is capped by the per-task w^2 budget.
    const long double diagonalSignal = std::min(
            profile.diagonalP99 > 0.0
                    ? static_cast<long double>(profile.diagonalP99)
                    : static_cast<long double>(profile.maximumDiagonalJump),
            profile.maximumDiagonalJump > 0.0
                    ? static_cast<long double>(profile.maximumDiagonalJump)
                    : static_cast<long double>(profile.diagonalP99));
    const long double diagonalFraction = 0.20L +
            0.30L * static_cast<long double>(profile.uncertainty);
    uint64_t certifiedInitialBand = std::max<uint64_t>(
            lengthDifference,
            roundedSaturating(diagonalSignal * diagonalFraction + 64.0L));
    certifiedInitialBand = std::max<uint64_t>(64, certifiedInitialBand);
    certifiedInitialBand = std::min<uint64_t>(certifiedInitialBand,
                                              maximumLength);
    uint64_t certifiedInitialMemory = scoreOnlyKswMemory;
    // This path is constructed only when both high-WFA modes and full-KSW2 are
    // structurally unavailable.  The process scheduler may delay this large
    // exact rescue candidate, but it cannot expand the candidate beyond w^2.
    uint64_t certifiedBudget = 0;
    uint64_t certifiedMaximumBand = 0;
    bool scoreCertifiedKswFeasible = false;
    if (fastExactStructurallyUnavailable) {
        const uint64_t certifiedInitialTraceMemory =
                ksw2TracebackMemoryBytes(
                        profile.queryLength, profile.referenceLength,
                        certifiedInitialBand);
        certifiedInitialMemory = saturatingAdd(
                scoreOnlyKswMemory, certifiedInitialTraceMemory);
        certifiedBudget = fullKswBudget;
        certifiedMaximumBand = ksw2MaximumBandForBudget(
                profile.queryLength, profile.referenceLength,
                scoreOnlyKswMemory, certifiedInitialBand,
                certifiedBudget);
        scoreCertifiedKswFeasible = lengthsFitInt &&
                certifiedBudget > 0 &&
                certifiedInitialBand <= certifiedMaximumBand;
    }
    const uint64_t scoreCertifiedKswMemory =
            scoreCertifiedKswFeasible
            ? ksw2CertifiedReservationBytes(
                      scoreOnlyKswMemory,
                      profile.queryLength,
                      profile.referenceLength,
                      certifiedInitialBand,
                      certifiedMaximumBand)
            : certifiedInitialMemory;
    // For an approximate banded attempt, spend the complete per-task w^2
    // allowance on the widest traceback band that fits.  Endpoint reachability
    // requires at least |q-r|; a narrower band is not a valid candidate.  Reuse
    // the already-computed certified cap when possible, otherwise invert the
    // architecture/page-aware KSW2 resident model.
    const uint64_t minimumReachableBand = std::max<uint64_t>(
            1, lengthDifference);
    uint64_t allowedBand = 0;
    if (lengthsFitInt && windowWidth > 0 && maximumLength > 0 &&
        minimumReachableBand <= maximumLength && workerBudget > 0) {
        if (fullKswMemory <= workerBudget) {
            allowedBand = maximumLength;
        } else if (certifiedMaximumBand >= minimumReachableBand) {
            allowedBand = certifiedMaximumBand;
        } else {
            allowedBand = ksw2MaximumBandForBudget(
                    profile.queryLength, profile.referenceLength, 0,
                    minimumReachableBand, workerBudget);
        }
    }
    const bool bandedKswFeasible = lengthsFitInt && allowedBand > 0 &&
            allowedBand >= lengthDifference &&
            allowedBand <= static_cast<uint64_t>(INT_MAX);
    const uint64_t bandMemory = bandedKswFeasible
                                ? ksw2TracebackMemoryBytes(
                                          profile.queryLength,
                                          profile.referenceLength,
                                          allowedBand)
                                : 0;
    // KSW2's optimized match/mismatch kernel treats the fifth alphabet
    // symbol as a wildcard.  N-containing tasks therefore use the generic
    // scoring kernel so their objective remains identical to WFA.  The first
    // B73/Mo17 calibration measured a ~1.86x cost for that path; round upward
    // to keep its completion-time ranking conservative.
    const long double kswScoringMultiplier =
            profile.ambiguousBaseFraction > 0.0 ? 1.90L : 1.0L;
    const long double fullKswWork = static_cast<long double>(
            profile.queryLength) * profile.referenceLength *
            kswScoringMultiplier;
    const uint64_t certifiedFirstTraceBand = scoreCertifiedKswFeasible
            ? std::min<uint64_t>(
                      certifiedMaximumBand,
                      certifiedInitialBand >
                              std::numeric_limits<uint64_t>::max() / 2
                      ? certifiedMaximumBand
                      : std::max<uint64_t>(
                                certifiedInitialBand,
                                certifiedInitialBand * 2))
            : certifiedInitialBand;
    const uint64_t certifiedP50Band = scoreCertifiedKswFeasible
            ? std::min<uint64_t>(
                      certifiedMaximumBand,
                      std::max<uint64_t>(
                              certifiedFirstTraceBand,
                              profile.estimatedBandWidth))
            : certifiedFirstTraceBand;
    const uint64_t certifiedP90Band = scoreCertifiedKswFeasible
            ? std::min<uint64_t>(
                      certifiedMaximumBand,
                      std::max<uint64_t>(
                              certifiedP50Band,
                              roundedSaturating(
                                      static_cast<long double>(
                                              profile.estimatedBandWidth) *
                                      (1.0L + profile.uncertainty))))
            : certifiedFirstTraceBand;
    const long double certifiedP50GeometricTraceWork =
            ksw2GeometricTraceWork(
                    maximumLength, certifiedFirstTraceBand,
                    certifiedP50Band, kswScoringMultiplier);
    const long double certifiedP90TraceWork = ksw2GeometricTraceWork(
            maximumLength, certifiedFirstTraceBand,
            certifiedP90Band, kswScoringMultiplier);
    // P50 includes the observed chance that the sparse exact-anchor chain
    // certifies before monolithic DP. P90 prices the complete geometric path
    // plus the smallest-band best-effort recovery. The chain itself is
    // deliberately omitted from this analytic term because its independently
    // traced segments are small compared with the monolithic matrices;
    // runtime residual calibration still protects the upper quantile.
    const long double certifiedSmallestTraceWork =
            scoreCertifiedKswFeasible &&
                    certifiedInitialBand < certifiedFirstTraceBand
            ? 2.0L * static_cast<long double>(maximumLength) *
                      static_cast<long double>(certifiedInitialBand) *
                      kswScoringMultiplier
            : 0.0L;
    const long double certifiedP50TraceWork =
            (1.0L - kCertifiedChainSuccessProbability) *
            certifiedP50GeometricTraceWork +
            0.60L * certifiedSmallestTraceWork;
    const long double certifiedP90TraceWorkWithFallback =
            certifiedP90TraceWork + certifiedSmallestTraceWork;
    // The first 100 Mb AVX512 B73/Mo17 run observed 8/20 successful
    // certifications. Coverage of the chain-estimated diagonal envelope helps,
    // but is not a guarantee in repetitive maize sequence. Replace the former
    // unsupported 70--97% prior with a bounded, coverage-aware estimate.
    const double certifiedBandCoverage = profile.estimatedBandWidth == 0
            ? 1.0 : std::min(
                    1.0,
                    static_cast<double>(certifiedMaximumBand) /
                    static_cast<double>(profile.estimatedBandWidth));
    const double certifiedSuccessProbability = std::max(
            0.30, std::min(
                    0.70, 0.32 + 0.28 * certifiedBandCoverage +
                          0.10 * profile.confidence));
    double nextExactP50Minutes = std::numeric_limits<double>::infinity();
    double nextExactP90Minutes = std::numeric_limits<double>::infinity();
    if (mediumWfaMemoryFeasible) {
        nextExactP50Minutes = mediumWfaEstimatedMinutes;
        nextExactP90Minutes = mediumWfaEstimatedMinutesP90;
    }
    if (lowWfaMemoryFeasible &&
        lowWfaEstimatedMinutes < nextExactP50Minutes) {
        nextExactP50Minutes = lowWfaEstimatedMinutes;
        nextExactP90Minutes = lowWfaEstimatedMinutesP90;
    }
    const double fallbackP50Minutes = std::isfinite(nextExactP50Minutes)
            ? nextExactP50Minutes : 0.0;
    const double fallbackP90Minutes = std::isfinite(nextExactP90Minutes)
            ? nextExactP90Minutes : 0.0;
    const long double scoreCertifiedKswWork = fullKswWork +
            certifiedP50TraceWork;
    const long double scoreCertifiedKswP90Work = fullKswWork +
            certifiedP90TraceWorkWithFallback;
    const long double bandKswWork = static_cast<long double>(maximumLength) *
                                    allowedBand * kswScoringMultiplier;
    const double fullKswEstimatedMinutes = finiteMinutes(
            fullKswWork, kKsw2CellsPerSecond, kKsw2SetupSeconds);
    const double scoreCertifiedAttemptMinutes = finiteMinutes(
            scoreCertifiedKswWork, kKsw2CellsPerSecond,
            2.0 * kKsw2SetupSeconds);
    const double scoreCertifiedKswEstimatedMinutes =
            scoreCertifiedAttemptMinutes +
            (1.0 - certifiedSuccessProbability) * fallbackP50Minutes;
    const double bandKswEstimatedMinutes = finiteMinutes(
            bandKswWork, kKsw2CellsPerSecond, kKsw2SetupSeconds);
    const double slidingEstimatedMinutes = finiteMinutes(
            slidingWork * kswScoringMultiplier,
            kSlidingCellsPerSecond, kSlidingSetupSeconds);

    // Predict score degradation rather than assuming that the nominal
    // Banded KSW2 and the sliding decomposition share one fallback-quality
    // tier.  For KSW2 the principal failure signal is a band narrower than
    // the chain-estimated diagonal excursion.  For sliding alignment it is
    // the number of forced chunk boundaries, weighted by the probability
    // that a structural event crosses a boundary.  Both quantities are
    // expressed in objective-score units, so the selector can prefer the
    // result expected to be closer to full DP before using runtime as a
    // tie-breaker.
    const double cheapestExtension = static_cast<double>(std::max<uint64_t>(
            1, std::min(gapExtend1, gapExtend2)));
    const double cheapestOpen = static_cast<double>(std::max<uint64_t>(
            1, std::min(gapOpen1, gapOpen2)));
    const double mismatchUnit = static_cast<double>(std::max<uint64_t>(
            1, mismatch));
    const double expectedPathWidth = static_cast<double>(std::max<uint64_t>(
            lengthDifference, profile.estimatedBandWidth));
    const double tailPathWidth = std::max(
            expectedPathWidth,
            expectedPathWidth * (1.0 + profile.uncertainty) +
                    profile.maximumDiagonalJump *
                    (0.25 + 0.75 * profile.uncertainty));
    const double bandShortfall = std::max(
            0.0, expectedPathWidth - static_cast<double>(allowedBand));
    const double bandTailShortfall = std::max(
            bandShortfall,
            tailPathWidth - static_cast<double>(allowedBand));
    // A chain whose point estimate lies inside the band is not proof that the
    // band contains the best path.  In low-Jaccard/high-uncertainty maize
    // intervals the old model assigned zero loss and confidently selected a
    // catastrophically truncated band.  Charge a small unresolved-path floor
    // equivalent to 64 uncertain mismatches at P50 and 128 at P90.  On 502
    // paired B73/Mo17 fallback intervals this removed every large false-banded
    // decision while retaining the high-confidence banded wins.
    const double unresolvedPathProbability = std::max(
            0.0, std::min(
                    1.0,
                    (1.0 - profile.sketchJaccard) * profile.uncertainty));
    const double bandedUnresolvedLoss = unresolvedPathProbability *
            (cheapestOpen + mismatchUnit * 64.0);
    const double bandedUnresolvedLossP90 = unresolvedPathProbability *
            (cheapestOpen + mismatchUnit * 128.0);
    const double bandedScoreLoss = bandedKswFeasible
            ? std::max(
                      bandShortfall * cheapestExtension *
                              (0.75 + 0.50 * profile.uncertainty),
                      bandedUnresolvedLoss)
            : std::numeric_limits<double>::infinity();
    const double bandedScoreLossP90 = bandedKswFeasible
            ? std::max(
                      bandTailShortfall * cheapestExtension +
                              (bandTailShortfall > 0.0 ? cheapestOpen : 0.0),
                      bandedUnresolvedLossP90)
            : std::numeric_limits<double>::infinity();

    const uint64_t slidingChunks = slidingWidth == 0
            ? std::numeric_limits<uint64_t>::max()
            : divideRoundUp(maximumLength, slidingWidth);
    const uint64_t slidingBoundaries = slidingChunks ==
            std::numeric_limits<uint64_t>::max()
            ? slidingChunks : (slidingChunks > 0 ? slidingChunks - 1 : 0);
    const double diagonalRisk = maximumLength == 0 ? 0.0 : std::min(
            1.0, (profile.diagonalP99 + profile.maximumDiagonalJump) /
                         static_cast<double>(maximumLength));
    const double boundaryRisk = std::min(
            1.0, 0.08 + 0.62 * profile.uncertainty +
                         0.20 * profile.lowComplexityFraction +
                         0.35 * diagonalRisk);
    const double slidingBoundaryPenalty = cheapestOpen +
            mismatchUnit * (8.0 + 24.0 * boundaryRisk);
    const double slidingScoreLoss = slidingBoundaries ==
            std::numeric_limits<uint64_t>::max()
            ? std::numeric_limits<double>::infinity()
            : static_cast<double>(slidingBoundaries) * boundaryRisk *
                      slidingBoundaryPenalty;
    const double slidingScoreLossP90 = slidingBoundaries ==
            std::numeric_limits<uint64_t>::max()
            ? std::numeric_limits<double>::infinity()
            : static_cast<double>(slidingBoundaries) *
                      std::min(1.0, boundaryRisk + 0.35) *
                      (cheapestOpen + mismatchUnit * 64.0);
    const auto p90 = [&profile](AlignmentCandidate candidate,
                                double p50) {
        return p50 * p90Inflation(profile, candidate);
    };
    const double fullKswEstimatedMinutesP90 = p90(
            AlignmentCandidate::Ksw2Full, fullKswEstimatedMinutes);
    const double scoreCertifiedKswEstimatedMinutesP90 = std::max(
            scoreCertifiedKswEstimatedMinutes,
            p90(AlignmentCandidate::Ksw2ScoreCertified,
                finiteMinutes(scoreCertifiedKswP90Work,
                              kKsw2CellsPerSecond,
                              2.0 * kKsw2SetupSeconds) +
                fallbackP90Minutes));
    const auto withinExactTimeLimit = [exactAlignmentMaximumMinutes](
            double p50, double p90Minutes) {
        return exactCandidateWithinTimeLimit(
                p50, p90Minutes, exactAlignmentMaximumMinutes);
    };
    const bool singletrackWfaWithinTimeLimit = withinExactTimeLimit(
            singletrackWfaEstimatedMinutes,
            singletrackWfaEstimatedMinutesP90);
    const bool standardWfaWithinTimeLimit = withinExactTimeLimit(
            standardWfaEstimatedMinutes,
            standardWfaEstimatedMinutesP90);
    const bool mediumWfaWithinTimeLimit = withinExactTimeLimit(
            mediumWfaEstimatedMinutes,
            mediumWfaEstimatedMinutesP90);
    const bool lowWfaWithinTimeLimit = withinExactTimeLimit(
            lowWfaEstimatedMinutes,
            lowWfaEstimatedMinutesP90);
    const bool fullKswWithinTimeLimit = withinExactTimeLimit(
            fullKswEstimatedMinutes, fullKswEstimatedMinutesP90);
    const bool scoreCertifiedKswWithinTimeLimit = withinExactTimeLimit(
            scoreCertifiedKswEstimatedMinutes,
            scoreCertifiedKswEstimatedMinutesP90);
    const double highSuccess = std::max(
            0.75, 0.985 - 0.12 * profile.uncertainty);
    const double succinctSuccess = std::max(
            0.82, 0.992 - 0.08 * profile.uncertainty);
    plan.estimates.push_back(estimate(
            AlignmentCandidate::SingletrackWfa,
            singletrackMemory, singletrackWfaWork,
            singletrackMemoryFeasible,
            workWithinSlidingEnvelope(singletrackWfaWork),
            singletrackWfaEstimatedMinutes,
            singletrackWfaEstimatedMinutesP90, highSuccess,
            singletrackWfaWithinTimeLimit,
            singletrackWfaWithinTimeLimit));
    plan.estimates.push_back(estimate(
            AlignmentCandidate::StandardWfa, standardMemory, standardWfaWork,
            standardMemoryFeasible,
            workWithinSlidingEnvelope(standardWfaWork),
            standardWfaEstimatedMinutes,
            standardWfaEstimatedMinutesP90, highSuccess,
            standardWfaWithinTimeLimit,
            standardWfaWithinTimeLimit));
    plan.estimates.push_back(estimate(
            AlignmentCandidate::MediumWfa, mediumWfaMemory, mediumWfaWork,
            mediumWfaMemoryFeasible,
            workWithinSlidingEnvelope(mediumWfaWork),
            mediumWfaEstimatedMinutes,
            mediumWfaEstimatedMinutesP90, succinctSuccess,
            mediumWfaWithinTimeLimit, mediumWfaWithinTimeLimit));
    plan.estimates.push_back(estimate(
            AlignmentCandidate::LowWfa, lowWfaMemory, lowWfaWork,
            lowWfaMemoryFeasible,
            workWithinSlidingEnvelope(lowWfaWork),
            lowWfaEstimatedMinutes,
            lowWfaEstimatedMinutesP90, succinctSuccess,
            lowWfaWithinTimeLimit, lowWfaWithinTimeLimit));
    plan.estimates.push_back(estimate(
            AlignmentCandidate::Ksw2ScoreCertified,
            scoreCertifiedKswMemory, scoreCertifiedKswWork,
            scoreCertifiedKswFeasible &&
                    fastExactStructurallyUnavailable,
            true,
            scoreCertifiedKswEstimatedMinutes,
            scoreCertifiedKswEstimatedMinutesP90,
            certifiedSuccessProbability,
            scoreCertifiedKswWithinTimeLimit,
            scoreCertifiedKswWithinTimeLimit &&
                    fastExactStructurallyUnavailable));
    plan.estimates.push_back(estimate(
            AlignmentCandidate::Ksw2Full, fullKswMemory, fullKswWork,
            fullKswFeasible, true, fullKswEstimatedMinutes,
            fullKswEstimatedMinutesP90, 0.998,
            fullKswWithinTimeLimit, fullKswWithinTimeLimit));
    plan.estimates.push_back(estimate(
            AlignmentCandidate::Ksw2Banded, bandMemory, bandKswWork,
            bandedKswFeasible, true, bandKswEstimatedMinutes,
            p90(AlignmentCandidate::Ksw2Banded,
                bandKswEstimatedMinutes), 0.995));
    plan.estimates.back().predictedScoreLoss = bandedScoreLoss;
    plan.estimates.back().predictedScoreLossP90 = bandedScoreLossP90;
    plan.estimates.push_back(estimate(
            AlignmentCandidate::SlidingWindow, slidingMemory, slidingWork,
            windowWidth > 0, true, slidingEstimatedMinutes,
            p90(AlignmentCandidate::SlidingWindow,
                slidingEstimatedMinutes), 0.999));
    plan.estimates.back().predictedScoreLoss = slidingScoreLoss;
    plan.estimates.back().predictedScoreLossP90 = slidingScoreLossP90;

    plan.exactCandidates.reserve(6);
    if (singletrackMemoryFeasible && singletrackWfaWithinTimeLimit) {
        plan.exactCandidates.push_back(AlignmentCandidate::SingletrackWfa);
    }
    if (standardMemoryFeasible && standardWfaWithinTimeLimit) {
        plan.exactCandidates.push_back(AlignmentCandidate::StandardWfa);
    }
    if (mediumWfaMemoryFeasible && mediumWfaWithinTimeLimit) {
        plan.exactCandidates.push_back(AlignmentCandidate::MediumWfa);
    }
    if (lowWfaMemoryFeasible && lowWfaWithinTimeLimit) {
        plan.exactCandidates.push_back(AlignmentCandidate::LowWfa);
    }
    if (fastExactStructurallyUnavailable &&
        scoreCertifiedKswFeasible &&
        scoreCertifiedKswWithinTimeLimit) {
        plan.ksw2CertifiedInitialBandWidth = static_cast<int64_t>(
                certifiedInitialBand);
        plan.ksw2CertifiedMaximumBandWidth = static_cast<int64_t>(
                certifiedMaximumBand);
        plan.exactCandidates.push_back(
                AlignmentCandidate::Ksw2ScoreCertified);
    }
    if (fullKswFeasible && fullKswWithinTimeLimit) {
        plan.exactCandidates.push_back(AlignmentCandidate::Ksw2Full);
    }
    if (bandedKswFeasible) {
        plan.ksw2BandWidth = static_cast<int64_t>(allowedBand);
        plan.bandedCandidates.push_back(AlignmentCandidate::Ksw2Banded);
    }
    plan.slidingWindowWidth = static_cast<int64_t>(slidingWidth);
    plan.lastResortCandidates.push_back(AlignmentCandidate::SlidingWindow);
    plan.approximateCandidates = plan.bandedCandidates;
    plan.approximateCandidates.insert(
            plan.approximateCandidates.end(),
            plan.lastResortCandidates.begin(),
            plan.lastResortCandidates.end());
    // Banded KSW2 and sliding-window alignment share the approximate tier.
    // Paired B73/Mo17 measurements across 502 independent fallback intervals
    // showed that a small hand-written boolean classifier systematically
    // reversed their quality order.  Rank them by the selector's predicted
    // score loss instead: the lower risk-adjusted loss is preferred, and
    // runtime breaks only an equal-quality tie.  Production normally executes
    // just this first candidate, so this does not pay for both alignments.
    const auto approximateQuality = [&plan](AlignmentCandidate candidate) {
        const auto found = std::find_if(
                plan.estimates.begin(), plan.estimates.end(),
                [candidate](const AlgorithmCostEstimate &cost) {
                    return cost.candidate == candidate;
                });
        if (found == plan.estimates.end()) {
            return std::make_pair(
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity());
        }
        const double p50 = std::max(0.0, found->predictedScoreLoss);
        const double p90Loss = std::max(p50, found->predictedScoreLossP90);
        const double riskAdjustedLoss = p50 + 0.35 * (p90Loss - p50);
        return std::make_pair(riskAdjustedLoss,
                              alignmentRiskAdjustedMinutes(*found));
    };
    std::stable_sort(
            plan.approximateCandidates.begin(),
            plan.approximateCandidates.end(),
            [&approximateQuality](
                    AlignmentCandidate first,
                    AlignmentCandidate second) {
                return approximateQuality(first) < approximateQuality(second);
            });
    return plan;
}

bool alignmentSelectionPlanMatches(
        const AlignmentSelectionPlan &plan,
        const std::string &query,
        const std::string &reference,
        int64_t windowWidth,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        double exactAlignmentMaximumEstimatedMinutes,
        uint64_t elasticHighWfaMemoryBudgetBytes,
        uint64_t elasticFullKsw2MemoryBudgetBytes,
        bool memorySchedulingEnabled) {
    const AlignmentSelectionProvenance &source = plan.provenance;
    return source.valid &&
           plan.profile.queryLength == query.size() &&
           plan.profile.referenceLength == reference.size() &&
           source.queryFingerprint == sequenceFingerprint(query) &&
           source.referenceFingerprint == sequenceFingerprint(reference) &&
           source.windowWidth == windowWidth &&
           source.mismatchingPenalty == mismatchingPenalty &&
           source.openGapPenalty1 == openGapPenalty1 &&
           source.extendGapPenalty1 == extendGapPenalty1 &&
           source.openGapPenalty2 == openGapPenalty2 &&
           source.extendGapPenalty2 == extendGapPenalty2 &&
           source.exactAlignmentMaximumEstimatedMinutes ==
                   exactAlignmentMaximumEstimatedMinutes &&
           source.elasticHighWfaMemoryBudgetBytes ==
                   elasticHighWfaMemoryBudgetBytes &&
           source.elasticFullKsw2MemoryBudgetBytes ==
                   elasticFullKsw2MemoryBudgetBytes &&
           source.memorySchedulingEnabled == memorySchedulingEnabled;
}

void resetAlignmentSelectionTelemetry() {
    evaluatedIntervals.store(0, std::memory_order_relaxed);
    exactTierIntervals.store(0, std::memory_order_relaxed);
    bandedOnlyIntervals.store(0, std::memory_order_relaxed);
    slidingOnlyIntervals.store(0, std::memory_order_relaxed);
    singletrackWfaMemoryRejects.store(0, std::memory_order_relaxed);
    singletrackWfaWorkWarnings.store(0, std::memory_order_relaxed);
    singletrackWfaTimeRejects.store(0, std::memory_order_relaxed);
    standardWfaMemoryRejects.store(0, std::memory_order_relaxed);
    standardWfaWorkWarnings.store(0, std::memory_order_relaxed);
    standardWfaTimeRejects.store(0, std::memory_order_relaxed);
    mediumWfaMemoryRejects.store(0, std::memory_order_relaxed);
    mediumWfaWorkWarnings.store(0, std::memory_order_relaxed);
    mediumWfaTimeRejects.store(0, std::memory_order_relaxed);
    lowWfaMemoryRejects.store(0, std::memory_order_relaxed);
    lowWfaWorkWarnings.store(0, std::memory_order_relaxed);
    lowWfaTimeRejects.store(0, std::memory_order_relaxed);
    scoreCertifiedKsw2MemoryRejects.store(0,
                                          std::memory_order_relaxed);
    scoreCertifiedKsw2TimeRejects.store(0,
                                        std::memory_order_relaxed);
    fullKsw2MemoryRejects.store(0, std::memory_order_relaxed);
    fullKsw2TimeRejects.store(0, std::memory_order_relaxed);
    bandedKsw2MemoryRejects.store(0, std::memory_order_relaxed);
    exactRuntimeMemoryFailures.store(0, std::memory_order_relaxed);
    exactRuntimeOtherFailures.store(0, std::memory_order_relaxed);
    bandedFallbackExecutions.store(0, std::memory_order_relaxed);
    slidingFallbackExecutions.store(0, std::memory_order_relaxed);
}

void recordAlignmentSelectionPlan(const AlignmentSelectionPlan &plan) {
    evaluatedIntervals.fetch_add(1, std::memory_order_relaxed);
    if (!plan.exactCandidates.empty()) {
        exactTierIntervals.fetch_add(1, std::memory_order_relaxed);
    } else if (!plan.bandedCandidates.empty()) {
        bandedOnlyIntervals.fetch_add(1, std::memory_order_relaxed);
    } else {
        slidingOnlyIntervals.fetch_add(1, std::memory_order_relaxed);
    }
    const AlgorithmCostEstimate *singletrack = nullptr;
    const AlgorithmCostEstimate *standard = nullptr;
    const AlgorithmCostEstimate *medium = nullptr;
    const AlgorithmCostEstimate *low = nullptr;
    const AlgorithmCostEstimate *scoreCertifiedKsw2 = nullptr;
    const AlgorithmCostEstimate *fullKsw2 = nullptr;
    const AlgorithmCostEstimate *bandedKsw2 = nullptr;
    for (const AlgorithmCostEstimate &item : plan.estimates) {
        switch (item.candidate) {
            case AlignmentCandidate::SingletrackWfa:
                singletrack = &item;
                break;
            case AlignmentCandidate::StandardWfa: standard = &item; break;
            case AlignmentCandidate::MediumWfa: medium = &item; break;
            case AlignmentCandidate::LowWfa: low = &item; break;
            case AlignmentCandidate::Ksw2ScoreCertified:
                scoreCertifiedKsw2 = &item;
                break;
            case AlignmentCandidate::Ksw2Full: fullKsw2 = &item; break;
            case AlignmentCandidate::Ksw2Banded: bandedKsw2 = &item; break;
            case AlignmentCandidate::SlidingWindow: break;
        }
    }
    if (singletrack != nullptr) {
        if (!singletrack->memoryFeasible) {
            singletrackWfaMemoryRejects.fetch_add(
                    1, std::memory_order_relaxed);
        }
        if (singletrack->memoryFeasible && !singletrack->workFeasible) {
            singletrackWfaWorkWarnings.fetch_add(
                    1, std::memory_order_relaxed);
        }
        if (singletrack->memoryFeasible && !singletrack->timeFeasible) {
            singletrackWfaTimeRejects.fetch_add(
                    1, std::memory_order_relaxed);
        }
    }
    if (standard != nullptr) {
        if (!standard->memoryFeasible) {
            standardWfaMemoryRejects.fetch_add(1, std::memory_order_relaxed);
        }
        if (standard->memoryFeasible && !standard->workFeasible) {
            standardWfaWorkWarnings.fetch_add(1, std::memory_order_relaxed);
        }
        if (standard->memoryFeasible && !standard->timeFeasible) {
            standardWfaTimeRejects.fetch_add(1, std::memory_order_relaxed);
        }
    }
    if (medium != nullptr) {
        if (!medium->memoryFeasible) {
            mediumWfaMemoryRejects.fetch_add(1, std::memory_order_relaxed);
        }
        if (medium->memoryFeasible && !medium->workFeasible) {
            mediumWfaWorkWarnings.fetch_add(1, std::memory_order_relaxed);
        }
        if (medium->memoryFeasible && !medium->timeFeasible) {
            mediumWfaTimeRejects.fetch_add(1, std::memory_order_relaxed);
        }
    }
    if (low != nullptr) {
        if (!low->memoryFeasible) {
            lowWfaMemoryRejects.fetch_add(1, std::memory_order_relaxed);
        }
        if (low->memoryFeasible && !low->workFeasible) {
            lowWfaWorkWarnings.fetch_add(1, std::memory_order_relaxed);
        }
        if (low->memoryFeasible && !low->timeFeasible) {
            lowWfaTimeRejects.fetch_add(1, std::memory_order_relaxed);
        }
    }
    if (scoreCertifiedKsw2 == nullptr ||
        !scoreCertifiedKsw2->memoryFeasible) {
        scoreCertifiedKsw2MemoryRejects.fetch_add(
                1, std::memory_order_relaxed);
    }
    if (scoreCertifiedKsw2 != nullptr &&
        scoreCertifiedKsw2->memoryFeasible &&
        !scoreCertifiedKsw2->timeFeasible) {
        scoreCertifiedKsw2TimeRejects.fetch_add(
                1, std::memory_order_relaxed);
    }
    if (fullKsw2 == nullptr || !fullKsw2->memoryFeasible) {
        fullKsw2MemoryRejects.fetch_add(1, std::memory_order_relaxed);
    }
    if (fullKsw2 != nullptr && fullKsw2->memoryFeasible &&
        !fullKsw2->timeFeasible) {
        fullKsw2TimeRejects.fetch_add(1, std::memory_order_relaxed);
    }
    if (bandedKsw2 == nullptr || !bandedKsw2->memoryFeasible) {
        bandedKsw2MemoryRejects.fetch_add(1, std::memory_order_relaxed);
    }
}

void recordExactAlignmentRuntimeFailure(bool memoryFailure) {
    (memoryFailure ? exactRuntimeMemoryFailures : exactRuntimeOtherFailures)
            .fetch_add(1, std::memory_order_relaxed);
}

void recordBandedFallbackExecution() {
    bandedFallbackExecutions.fetch_add(1, std::memory_order_relaxed);
}

void recordSlidingFallbackExecution() {
    slidingFallbackExecutions.fetch_add(1, std::memory_order_relaxed);
}

AlignmentSelectionTelemetry alignmentSelectionTelemetrySnapshot() {
    AlignmentSelectionTelemetry snapshot;
    snapshot.evaluatedIntervals = evaluatedIntervals.load(
            std::memory_order_relaxed);
    snapshot.exactTierIntervals = exactTierIntervals.load(
            std::memory_order_relaxed);
    snapshot.bandedOnlyIntervals = bandedOnlyIntervals.load(
            std::memory_order_relaxed);
    snapshot.slidingOnlyIntervals = slidingOnlyIntervals.load(
            std::memory_order_relaxed);
    snapshot.singletrackWfaMemoryRejects =
            singletrackWfaMemoryRejects.load(std::memory_order_relaxed);
    snapshot.singletrackWfaWorkWarnings =
            singletrackWfaWorkWarnings.load(std::memory_order_relaxed);
    snapshot.singletrackWfaTimeRejects =
            singletrackWfaTimeRejects.load(std::memory_order_relaxed);
    snapshot.standardWfaMemoryRejects = standardWfaMemoryRejects.load(
            std::memory_order_relaxed);
    snapshot.standardWfaWorkWarnings = standardWfaWorkWarnings.load(
            std::memory_order_relaxed);
    snapshot.standardWfaTimeRejects = standardWfaTimeRejects.load(
            std::memory_order_relaxed);
    snapshot.mediumWfaMemoryRejects = mediumWfaMemoryRejects.load(
            std::memory_order_relaxed);
    snapshot.mediumWfaWorkWarnings = mediumWfaWorkWarnings.load(
            std::memory_order_relaxed);
    snapshot.mediumWfaTimeRejects = mediumWfaTimeRejects.load(
            std::memory_order_relaxed);
    snapshot.lowWfaMemoryRejects = lowWfaMemoryRejects.load(
            std::memory_order_relaxed);
    snapshot.lowWfaWorkWarnings = lowWfaWorkWarnings.load(
            std::memory_order_relaxed);
    snapshot.lowWfaTimeRejects = lowWfaTimeRejects.load(
            std::memory_order_relaxed);
    snapshot.scoreCertifiedKsw2MemoryRejects =
            scoreCertifiedKsw2MemoryRejects.load(
                    std::memory_order_relaxed);
    snapshot.scoreCertifiedKsw2TimeRejects =
            scoreCertifiedKsw2TimeRejects.load(
                    std::memory_order_relaxed);
    snapshot.fullKsw2MemoryRejects = fullKsw2MemoryRejects.load(
            std::memory_order_relaxed);
    snapshot.fullKsw2TimeRejects = fullKsw2TimeRejects.load(
            std::memory_order_relaxed);
    snapshot.bandedKsw2MemoryRejects = bandedKsw2MemoryRejects.load(
            std::memory_order_relaxed);
    snapshot.exactRuntimeMemoryFailures = exactRuntimeMemoryFailures.load(
            std::memory_order_relaxed);
    snapshot.exactRuntimeOtherFailures = exactRuntimeOtherFailures.load(
            std::memory_order_relaxed);
    snapshot.bandedFallbackExecutions = bandedFallbackExecutions.load(
            std::memory_order_relaxed);
    snapshot.slidingFallbackExecutions = slidingFallbackExecutions.load(
            std::memory_order_relaxed);
    return snapshot;
}

void configureAlignmentTraceFile(const std::string &path) {
    traceEnabled.store(false, std::memory_order_release);
    std::lock_guard<std::mutex> lock(traceMutex);
    if (traceStream.is_open()) {
        traceStream.flush();
        traceStream.close();
    }
    traceLinesSinceFlush = 0;
    nextTraceIntervalId.store(1, std::memory_order_relaxed);
    if (path.empty()) {
        return;
    }
    traceStream.open(path, std::ios::out | std::ios::trunc);
    if (!traceStream) {
        throw std::runtime_error(
                "cannot open alignment trace file: " + path);
    }
    traceStream
            << "model_version\tinterval_id\tattempt\tquery_length"
               "\treference_length\testimated_score\tconservative_score"
               "\tconfidence\tuncertainty\tsketch_jaccard"
               "\tambiguous_base_fraction\tlow_complexity_fraction"
               "\tdiagonal_p90\tdiagonal_p99"
               "\tmaximum_diagonal_jump\tcandidate"
               "\tpredicted_minutes_p50\tpredicted_minutes_p90"
               "\tpredicted_memory_bytes\treserved_memory_bytes"
               "\tpredicted_wait_minutes\tactual_seconds"
               "\tconfigured_exact_limit_minutes"
               "\texact_runtime_spent_seconds"
               "\texact_runtime_remaining_seconds"
               "\tactual_score\tcertified_optimal_score"
               "\tcertified_initial_band\tcertified_maximum_band"
               "\tcertified_final_band\tcertified_traceback_attempts"
               "\tactual_memory_bytes\tprocess_rss_bytes"
               "\tprocess_memory_limit_bytes\tactive_reserved_bytes"
               "\tworker_count\tactive_tasks\tready_tasks\tdeferred_tasks"
               "\tplanned_tasks\tschedulable_tasks\tglobal_future_tasks"
               "\tblocked_ordered_heads\toutstanding_tasks"
               "\toutstanding_estimated_cost\tcritical_estimated_cost"
               "\tglobal_tail_phase\tadmission_tail_phase"
               "\ttail_phase\tresult_method\tstatus\n";
    traceStream.flush();
    traceEnabled.store(true, std::memory_order_release);
}

bool alignmentTraceEnabled() {
    return traceEnabled.load(std::memory_order_acquire);
}

uint64_t nextAlignmentTraceIntervalId() {
    return nextTraceIntervalId.fetch_add(1, std::memory_order_relaxed);
}

void recordAlignmentAttemptTrace(const AlignmentAttemptTrace &record) {
    if (!traceEnabled.load(std::memory_order_acquire)) {
        return;
    }
    std::lock_guard<std::mutex> lock(traceMutex);
    if (!traceStream.is_open()) {
        return;
    }
    std::string status = record.status;
    std::replace(status.begin(), status.end(), '\t', ' ');
    std::replace(status.begin(), status.end(), '\n', ' ');
    const AlignmentProfile &profile = record.profile;
    std::string resultMethod = record.resultMethod;
    std::replace(resultMethod.begin(), resultMethod.end(), '\t', ' ');
    std::replace(resultMethod.begin(), resultMethod.end(), '\n', ' ');
    traceStream << "b73-mo17-v8"
                << '\t' << record.intervalId
                << '\t' << record.attempt
                << '\t' << profile.queryLength
                << '\t' << profile.referenceLength
                << '\t' << profile.estimatedScore
                << '\t' << profile.conservativeScore
                << '\t' << std::setprecision(9) << profile.confidence
                << '\t' << profile.uncertainty
                << '\t' << profile.sketchJaccard
                << '\t' << profile.ambiguousBaseFraction
                << '\t' << profile.lowComplexityFraction
                << '\t' << profile.diagonalP90
                << '\t' << profile.diagonalP99
                << '\t' << profile.maximumDiagonalJump
                << '\t' << alignmentCandidateName(record.candidate)
                << '\t' << record.predictedMinutesP50
                << '\t' << record.predictedMinutesP90
                << '\t' << record.predictedMemoryBytes
                << '\t' << record.reservedMemoryBytes
                << '\t' << record.predictedWaitMinutes
                << '\t' << record.actualSeconds
                << '\t' << record.configuredExactLimitMinutes
                << '\t' << record.exactRuntimeSpentSeconds
                << '\t' << record.exactRuntimeRemainingSeconds
                << '\t' << record.actualScore
                << '\t' << record.certifiedOptimalScore
                << '\t' << record.certifiedInitialBandWidth
                << '\t' << record.certifiedMaximumBandWidth
                << '\t' << record.certifiedFinalBandWidth
                << '\t' << record.certifiedTracebackAttempts
                << '\t' << record.actualMemoryBytes
                << '\t' << record.processResidentBytes
                << '\t' << record.processMemoryLimitBytes
                << '\t' << record.activeReservedBytes
                << '\t' << record.workerCount
                << '\t' << record.activeTasks
                << '\t' << record.readyTasks
                << '\t' << record.deferredTasks
                << '\t' << record.plannedTasks
                << '\t' << record.schedulableTasks
                << '\t' << record.globalFutureTasks
                << '\t' << record.blockedOrderedHeads
                << '\t' << record.outstandingTasks
                << '\t' << record.outstandingEstimatedCost
                << '\t' << record.criticalEstimatedCost
                << '\t' << (record.globalTailPhase ? 1 : 0)
                << '\t' << (record.admissionTailPhase ? 1 : 0)
                << '\t' << (record.tailPhase ? 1 : 0)
                << '\t' << resultMethod
                << '\t' << status << '\n';
    if (++traceLinesSinceFlush >= 32 || status == "started") {
        traceStream.flush();
        traceLinesSinceFlush = 0;
    }
}

}  // namespace anchorwave
