#include "AlignmentGapPlanner.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace anchorwave {

namespace {
constexpr uint64_t kMinimumParallelAlignmentGapCost = 100000000ULL;
}

uint64_t minimumParallelAlignmentGapCost() {
    return kMinimumParallelAlignmentGapCost;
}

uint64_t alignmentGapEstimatedCost(uint64_t referenceLength,
                                   uint64_t queryLength) {
    if (queryLength != 0 &&
        referenceLength >
                std::numeric_limits<uint64_t>::max() / queryLength) {
        return std::numeric_limits<uint64_t>::max();
    }
    return referenceLength * queryLength;
}

uint64_t alignmentGapRuntimePriorityCost(double predictedMinutes) {
    if (!(predictedMinutes > 0.0)) {
        return 1;
    }
    constexpr long double kMicrosecondsPerMinute = 60000000.0L;
    const long double cost = static_cast<long double>(predictedMinutes) *
                             kMicrosecondsPerMinute;
    if (cost >= static_cast<long double>(
                        std::numeric_limits<uint64_t>::max())) {
        return std::numeric_limits<uint64_t>::max();
    }
    return std::max<uint64_t>(1, static_cast<uint64_t>(std::ceil(cost)));
}

bool shouldParallelizeAlignmentGap(uint64_t referenceLength,
                                   uint64_t queryLength,
                                   int workerCount) {
    return workerCount > 1 && referenceLength != 0 && queryLength != 0 &&
           alignmentGapEstimatedCost(referenceLength, queryLength) >=
                   kMinimumParallelAlignmentGapCost;
}

std::vector<AlignmentGapDescriptor> planParallelInterAnchorGaps(
        const std::vector<AlignmentAnchorSpan> &anchors,
        int workerCount) {
    std::vector<AlignmentGapDescriptor> gaps;
    if (workerCount <= 1 || anchors.size() < 2) {
        return gaps;
    }

    gaps.reserve(anchors.size() - 1);
    for (std::size_t index = 1; index < anchors.size(); ++index) {
        const AlignmentAnchorSpan &previous = anchors[index - 1];
        const AlignmentAnchorSpan &current = anchors[index];
        if (previous.reverse != current.reverse ||
            previous.referenceEnd >= current.referenceStart) {
            continue;
        }

        AlignmentGapDescriptor gap;
        gap.anchorIndex = index;
        gap.referenceStart = previous.referenceEnd + 1;
        gap.referenceEnd = current.referenceStart - 1;
        gap.reverse = current.reverse;

        if (!current.reverse) {
            if (previous.queryEnd >= current.queryStart) {
                continue;
            }
            gap.queryStart = previous.queryEnd + 1;
            gap.queryEnd = current.queryStart - 1;
        } else {
            // Consecutive reverse-strand anchors must move towards smaller
            // query coordinates as reference coordinates increase.
            if (previous.queryStart <= current.queryEnd) {
                continue;
            }
            gap.queryStart = current.queryEnd + 1;
            gap.queryEnd = previous.queryStart - 1;
        }

        const uint64_t referenceLength =
                gap.referenceEnd - gap.referenceStart + 1;
        const uint64_t queryLength = gap.queryEnd - gap.queryStart + 1;
        gap.geometricCost = alignmentGapEstimatedCost(referenceLength,
                                                       queryLength);
        gap.estimatedCost = gap.geometricCost;
        if (shouldParallelizeAlignmentGap(referenceLength, queryLength,
                                          workerCount)) {
            gaps.push_back(gap);
        }
    }
    return gaps;
}

}  // namespace anchorwave
