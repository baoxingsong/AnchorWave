#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace anchorwave {

struct AlignmentSelectionPlan;

// Lightweight anchor geometry used by the sequence-alignment scheduler.  It is
// intentionally independent of AlignmentMatch so the planner can be tested
// without linking the complete AnchorWave model.
struct AlignmentAnchorSpan {
    uint64_t referenceStart = 0;
    uint64_t referenceEnd = 0;
    uint64_t queryStart = 0;
    uint64_t queryEnd = 0;
    bool reverse = false;
};

struct AlignmentGapDescriptor {
    std::size_t anchorIndex = 0;  // gap immediately before this anchor
    uint64_t referenceStart = 0;
    uint64_t referenceEnd = 0;
    uint64_t queryStart = 0;
    uint64_t queryEnd = 0;
    bool reverse = false;
    uint64_t geometricCost = 0;
    double predictedMinutes = 0.0;
    uint64_t predictedMemoryBytes = 0;
    // Used only by the bounded gap-prefetch heap. Keep it separate from
    // estimatedCost: the latter shares one executor queue with chromosome and
    // block tasks and must remain in the historical geometric-cost units.
    uint64_t schedulingPriorityCost = 0;
    uint64_t estimatedCost = 0;
    std::shared_ptr<const AlignmentSelectionPlan> selectionPlan;
    // Populated lazily by the submitted task's private descriptor copy. Every
    // memory-deferral retry reuses the same strings instead of repeating FASTA
    // I/O and allocation; planned-but-unsubmitted descriptors retain no input
    // cache, and completion releases the strings with the task closure.
    std::shared_ptr<std::string> preparedReferenceSequence;
    std::shared_ptr<std::string> preparedQuerySequence;
};

// The first implementation parallelizes only expensive, two-sided
// inter-anchor regions.  A product of 1e8 corresponds to roughly two 10-kb
// sequences and keeps task overhead small on the B73/Mo17 workload while
// exposing enough work to keep six alignment workers busy.
uint64_t minimumParallelAlignmentGapCost();
uint64_t alignmentGapEstimatedCost(uint64_t referenceLength,
                                   uint64_t queryLength);
uint64_t alignmentGapRuntimePriorityCost(double predictedMinutes);
bool shouldParallelizeAlignmentGap(uint64_t referenceLength,
                                   uint64_t queryLength,
                                   int workerCount);

std::vector<AlignmentGapDescriptor> planParallelInterAnchorGaps(
        const std::vector<AlignmentAnchorSpan> &anchors,
        int workerCount);

}  // namespace anchorwave
