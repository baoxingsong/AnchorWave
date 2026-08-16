#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace anchorwave {

enum class AlignmentCandidate {
    SingletrackWfa,
    StandardWfa,
    MediumWfa,
    LowWfa,
    Ksw2Singletrack,
    Ksw2ScoreCertified,
    Ksw2Full,
    Ksw2Banded,
    SlidingWindow
};

constexpr std::size_t kAlignmentCandidateCount =
        static_cast<std::size_t>(AlignmentCandidate::SlidingWindow) + 1;

struct AlignmentProfile {
    // Incremented whenever the published feature definition changes. Trace
    // and calibration data from different versions must not be pooled.
    uint32_t version = 2;
    uint64_t queryLength = 0;
    uint64_t referenceLength = 0;
    uint64_t sampledQueryKmers = 0;
    uint64_t sampledReferenceKmers = 0;
    uint64_t sampledQueryKmersK9 = 0;
    uint64_t sampledReferenceKmersK9 = 0;
    uint64_t sampledQueryKmersK21 = 0;
    uint64_t sampledReferenceKmersK21 = 0;
    uint64_t chainedAnchors = 0;
    uint64_t estimatedScore = 0;
    uint64_t conservativeScore = 0;
    uint64_t estimatedBandWidth = 0;
    double ambiguousBaseFraction = 0.0;
    // Mean normalized Shannon entropy of deterministic 64-base blocks. It is
    // zero for a homopolymer and one for equal A/C/G/T composition.
    double sequenceEntropy = 0.0;
    double lowComplexityFraction = 0.0;
    double sketchJaccard = 0.0;
    double sketchJaccardK9 = 0.0;
    double sketchJaccardK21 = 0.0;
    double sketchJaccardDispersion = 0.0;
    double chainedAnchorFraction = 0.0;
    double chainQueryCoverage = 0.0;
    double chainReferenceCoverage = 0.0;
    double chainSpanCoverage = 0.0;
    double chainGapP90 = 0.0;
    double uniqueAnchorFraction = 0.0;
    double estimatedMismatchRate = 0.0;
    double diagonalMad = 0.0;
    double diagonalP90 = 0.0;
    double diagonalP99 = 0.0;
    double maximumDiagonalJump = 0.0;
    double lengthRatio = 0.0;
    double uncertainty = 1.0;
    double confidence = 0.0;
    bool identical = false;
};

// Paper-facing decision record for the shared approximate tier. The response
// is modeled directly as D = score(KSW2-banded) - score(sliding-window), so a
// positive value recommends KSW2-banded. Scores use the same two-piece affine
// objective as Tier 1; because the Tier-1 optimum is an upper bound on both
// approximate scores, choosing the larger approximate score also minimizes
// score loss relative to that optimum.
struct ApproximateQualityDecision {
    uint32_t modelVersion = 0;
    AlignmentCandidate selected = AlignmentCandidate::SlidingWindow;
    double predictedBandedMinusSlidingScore = 0.0;
    double residualAdjustedBandedMinusSlidingScore = 0.0;
    double applicabilityMaximumAbsoluteZ = 0.0;
    bool calibrated = false;
};

struct AlgorithmCostEstimate {
    AlignmentCandidate candidate = AlignmentCandidate::SlidingWindow;
    uint64_t memoryBytes = 0;
    long double workUnits = 0.0L;
    // P50 is retained under the historical field name for source
    // compatibility.  The scheduler uses P90 for long-tail protection and
    // memoryBytes as the conservative admission estimate.
    double estimatedMinutes = 0.0;
    double estimatedMinutesP90 = 0.0;
    double successProbability = 1.0;
    // Approximate engines are ranked by their predicted loss relative to the
    // unrestricted dynamic-programming optimum.  A value of zero means that
    // the available evidence predicts the same optimum score.  Exact engines
    // always keep zero loss.  The P90 value makes an uncertain narrow band
    // lose to a more robust sliding decomposition even when their median
    // predictions are similar.
    double predictedScoreLoss = 0.0;
    double predictedScoreLossP90 = 0.0;
    bool memoryFeasible = false;
    bool workFeasible = false;
    // Balanced-mode (-bt > 0) comparison against the configured per-interval
    // time ceiling. Every duration quantile used for admission must stay at
    // or below -bt. With -bt 0 (exact-first), this flag is always true.
    bool withinConfiguredTimeLimit = true;
    bool timeFeasible = true;
    // Memory controls every candidate. Work remains advisory. In balanced
    // mode a positive -bt value may remove any exact-DP candidate whose
    // calibrated duration distribution is clearly beyond the ceiling.
    bool feasible = false;
};

// One state transition in the exact-tier stochastic shortest-path problem.
// attemptMinutes already includes the current resource wait and any
// makespan-oriented memory opportunity cost. Success terminates the interval;
// failure continues with one of the remaining exact candidates.
struct ExactAttemptSystemCost {
    AlignmentCandidate candidate = AlignmentCandidate::SlidingWindow;
    double attemptMinutes = 0.0;
    double successProbability = 1.0;
};

// A sparse, increasing chain of exact sampled k-mer matches.  The
// score-certified KSW2 path may use these matches to build a low-memory
// candidate alignment from independently traced segments.  The candidate is
// accepted only when its score equals the unrestricted score-only DP result,
// so an incorrect/noisy anchor can reduce certification rate but cannot
// reduce alignment quality.
struct AlignmentChainAnchor {
    uint32_t queryPosition = 0;
    uint32_t referencePosition = 0;
    uint16_t length = 0;
};

// Inputs that define a reusable selection plan. Sequence fingerprints protect
// against accidentally reusing a plan for different same-length strings;
// policy and resource fields protect against stale -w/-bt/-M or scoring
// settings. A mismatch is not an error: the caller simply rebuilds the plan.
struct AlignmentSelectionProvenance {
    uint64_t queryFingerprint = 0;
    uint64_t referenceFingerprint = 0;
    int64_t windowWidth = 0;
    int32_t mismatchingPenalty = 0;
    int32_t openGapPenalty1 = 0;
    int32_t extendGapPenalty1 = 0;
    int32_t openGapPenalty2 = 0;
    int32_t extendGapPenalty2 = 0;
    double exactAlignmentMaximumEstimatedMinutes = 0.0;
    uint64_t elasticHighWfaMemoryBudgetBytes = 0;
    uint64_t elasticFullKsw2MemoryBudgetBytes = 0;
    bool memorySchedulingEnabled = false;
    bool valid = false;
};

struct AlignmentSelectionPlan {
    AlignmentSelectionProvenance provenance;
    AlignmentProfile profile;
    std::vector<AlgorithmCostEstimate> estimates;
    // Two strict quality tiers. The executor must exhaust the exact tier
    // before entering the approximate tier. Within the exact tier the runtime
    // scheduler is free to
    // reorder candidates using predicted completion time and current resource
    // pressure; the vector order is only a stable deterministic tie-break.
    // Predicted work cannot remove an exact candidate.  A positive -bt value
    // selects balanced mode and may remove a clearly over-limit exact
    // candidate; -bt 0 selects exact-first mode and disables all time gates.
    std::vector<AlignmentCandidate> exactCandidates;
    // Compatibility views retained for callers and telemetry that distinguish
    // the two implementations. They are not separate quality tiers.
    std::vector<AlignmentCandidate> bandedCandidates;
    std::vector<AlignmentCandidate> lastResortCandidates;
    // Tier two contains both banded KSW2 and the historical sliding-window
    // decomposition. The selector predicts which one will have the larger
    // alignment score and executes that candidate first. A small prediction
    // error rate is accepted to avoid paying for both approximate alignments.
    std::vector<AlignmentCandidate> approximateCandidates;
    ApproximateQualityDecision approximateDecision;
    // Tier-level quality contract: when qlen*rlen fits in the raw -w^2 cell
    // envelope, the interval must remain in Tier 1 even if an implementation's
    // bookkeeping makes its byte estimate larger than w^2. This is not a
    // KSW2-full priority rule: every structurally valid exact engine competes
    // by expected wall time, and process-wide -M remains the hard admission
    // ceiling.
    bool exactCellEnvelopeException = false;
    std::vector<AlignmentChainAnchor> certificationAnchors;
    uint64_t ksw2CertifiedScoreOnlyMemoryBytes = 0;
    int64_t ksw2CertifiedInitialBandWidth = -1;
    int64_t ksw2CertifiedMaximumBandWidth = -1;
    int64_t ksw2BandWidth = -1;
    // Effective sliding-window chunk size after enforcing the same raw w^2
    // per-alignment memory ceiling used by every exact and banded engine.
    // This may be slightly smaller than the requested -w because KSW2 also
    // retains SIMD/page-rounded traceback metadata around the DP cells.
    int64_t slidingWindowWidth = -1;
};

struct AlignmentSelectionTelemetry {
    uint64_t evaluatedIntervals = 0;
    uint64_t exactTierIntervals = 0;
    uint64_t bandedOnlyIntervals = 0;
    uint64_t slidingOnlyIntervals = 0;
    uint64_t singletrackWfaMemoryRejects = 0;
    uint64_t singletrackWfaWorkWarnings = 0;
    uint64_t singletrackWfaTimeRejects = 0;
    uint64_t standardWfaMemoryRejects = 0;
    uint64_t standardWfaWorkWarnings = 0;
    uint64_t standardWfaTimeRejects = 0;
    uint64_t mediumWfaMemoryRejects = 0;
    uint64_t mediumWfaWorkWarnings = 0;
    uint64_t mediumWfaTimeRejects = 0;
    uint64_t lowWfaMemoryRejects = 0;
    uint64_t lowWfaWorkWarnings = 0;
    uint64_t lowWfaTimeRejects = 0;
    uint64_t ksw2SingletrackMemoryRejects = 0;
    uint64_t ksw2SingletrackTimeRejects = 0;
    uint64_t scoreCertifiedKsw2MemoryRejects = 0;
    uint64_t scoreCertifiedKsw2TimeRejects = 0;
    uint64_t fullKsw2MemoryRejects = 0;
    uint64_t fullKsw2TimeRejects = 0;
    uint64_t bandedKsw2MemoryRejects = 0;
    uint64_t exactRuntimeMemoryFailures = 0;
    uint64_t exactRuntimeOtherFailures = 0;
    uint64_t bandedFallbackExecutions = 0;
    uint64_t slidingFallbackExecutions = 0;
};

const char *alignmentCandidateName(AlignmentCandidate candidate);

// Public -b output intentionally hides the exact WFA execution mode so BED
// consumers see one stable method label. Internal selector/trace telemetry
// continues to use alignmentCandidateName() and therefore retains the
// Singletrack/high/medium/low distinction. The low-level standalone BiWFA API
// and benchmark live outside the production selector.
std::string alignmentMethodBedLabel(const std::string &internalMethod);

// Build a bounded, deterministic sequence sketch and predict the resource
// envelope of every available pairwise aligner. Candidates are returned in
// two strict quality tiers. Tier one contains Singletrack high-memory WFA,
// standard high/medium/low WFA, KSW2-Singletrack, full KSW2, and conditional
// score-certified KSW2.
// Score-certified KSW2 first computes the unbanded optimal score without a
// traceback matrix, then accepts a banded traceback only when its score is
// identical to that global optimum. The runtime scheduler
// dynamically ranks feasible tier-one alternatives by predicted start/finish
// time, memory opportunity cost, failure risk, and workload-tail state; the
// declaration order is only a deterministic tie-break. Tier two contains
// banded KSW2 and sliding-window alignment; the selector predicts which will
// produce the larger alignment score and normally executes only that method.
AlignmentSelectionPlan makeAlignmentSelectionPlan(
        const std::string &query,
        const std::string &reference,
        int64_t windowWidth,
        int32_t mismatchingPenalty,
        int32_t openGapPenalty1,
        int32_t extendGapPenalty1,
        int32_t openGapPenalty2,
        int32_t extendGapPenalty2,
        double exactAlignmentMaximumEstimatedMinutes = 0.0,
        // Optional lower raw-prediction ceiling retained for source
        // compatibility. Zero uses the -w-derived ceiling. A nonzero value is
        // clamped to w^2 except for the tier-wide qlen*rlen <= w^2 quality
        // contract; process-wide -M always remains the hard ceiling.
        uint64_t elasticHighWfaMemoryBudgetBytes = 0,
        // Optional lower ceiling for full/score-certified KSW2. It is normally
        // clamped to w^2 except for the tier-wide qlen*rlen <= w^2 quality
        // contract; process-wide -M still controls its admission.
        uint64_t elasticFullKsw2MemoryBudgetBytes = 0,
        bool memorySchedulingEnabled = false);

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
        uint64_t elasticHighWfaMemoryBudgetBytes = 0,
        uint64_t elasticFullKsw2MemoryBudgetBytes = 0,
        bool memorySchedulingEnabled = false);

// Analytic resident-memory model for the KSW2 traceback kernel. Unlike the
// virtual rectangular allocation, this counts the union of SIMD-rounded
// anti-diagonal pages that the kernel actually touches, plus its fully
// touched linear work arrays. Exposed to keep admission tests tied to the
// implementation rather than to an empirical discount factor.
uint64_t ksw2TracebackResidentMemoryBytes(
        uint64_t queryLength,
        uint64_t referenceLength,
        uint64_t bandWidth);

// Resident-memory model for the active SIMD KSW2-Singletrack traceback.
// Singletrack
// stores two byte difference tracks instead of KSW2's packed decision track;
// the score DP remains exact and uses the same -w^2 task ceiling.
uint64_t ksw2SingletrackResidentMemoryBytes(
        uint64_t queryLength,
        uint64_t referenceLength,
        uint64_t bandWidth);

// Shared runtime-ranking helpers. The risk-adjusted value is the cost of one
// attempt; failure-chain cost is added explicitly so it is never embedded in
// both an engine estimate and the dispatcher a second time.
double alignmentRiskAdjustedMinutes(const AlgorithmCostEstimate &estimate);
double alignmentExpectedCompletionMinutes(
        const AlgorithmCostEstimate &estimate,
        double fallbackExactMinutes);
double bestExactExpectedCompletionMinutes(
        const std::vector<AlgorithmCostEstimate> &estimates);
double exactAttemptExpectedSystemMinutes(
        const std::vector<ExactAttemptSystemCost> &attempts,
        std::size_t firstAttempt);
double bestExactExpectedSystemMinutes(
        const std::vector<ExactAttemptSystemCost> &attempts);

// Epsilon-Pareto dominance for candidates in the same exact tier.  A
// candidate is dominated when another needs no more memory or modeled time;
// a small reliability tolerance prevents statistically indistinguishable
// implementation priors from retaining an otherwise strictly worse engine.
bool exactCandidateEpsilonDominates(
        const AlgorithmCostEstimate &preferred,
        const AlgorithmCostEstimate &other,
        double successProbabilityTolerance = 0.0);

AlgorithmCostEstimate calibratedAlignmentCostEstimate(
        const AlgorithmCostEstimate &estimate,
        int activeTasks = 1,
        int workerCount = 1);
void recordAlignmentRuntimeObservation(
        AlignmentCandidate candidate,
        double predictedMinutesP50,
        double actualSeconds);
void recordAlignmentRuntimeObservation(
        AlignmentCandidate candidate,
        double predictedMinutesP50,
        double predictedMinutesP90,
        double actualSeconds,
        int activeTasks,
        int workerCount);
void resetAlignmentRuntimeCalibration();
double alignmentSelectionPriorityMinutes(
        const AlignmentSelectionPlan &plan);
uint64_t alignmentSelectionPriorityMemoryBytes(
        const AlignmentSelectionPlan &plan);

// A cooperative WFA abort is a censored observation: the true memory need is
// above the first budget, but not known exactly. Permit one Singletrack retry
// at the normal w^2 ceiling only when that ceiling is genuinely larger than
// both the attempted budget and the observed abort peak. Zero means no retry.
uint64_t expandedSingletrackRetryBudgetBytes(
        uint64_t attemptedBudgetBytes,
        uint64_t observedPeakBytes,
        uint64_t workerBudgetBytes,
        bool expandedRetryAlreadyUsed);

// Standard/medium/low WFA share the same score-model uncertainty that caused
// a censored Singletrack abort. They are temporarily withheld while the one
// expanded Singletrack retry competes with the independent KSW2 exact family.
bool sharesSingletrackMemoryUnderpredictionRisk(
        AlignmentCandidate candidate);

// After a censored Singletrack memory failure, these exact candidates provide
// an independent DP implementation and are tried before expanding the failed
// WFA attempt to the full worker budget.
bool isIndependentExactAfterSingletrackMemoryFailure(
        AlignmentCandidate candidate);

// A score-only pass followed by full KSW2 repeats the same O(mn) score DP.
// Prefer the one-pass full traceback whenever it is currently admissible.
// An observed Singletrack memory abort is the exception: the independent,
// low-resident score-certified path remains an early rescue candidate.
bool deferScoreCertifiedKsw2BehindFullKsw2(
        bool singletrackMemoryModelInvalidated,
        bool fullKsw2RuntimeAvailable);

// The minimizer-derived WFA score model becomes out-of-domain when ambiguous
// bases remove substantial sketch evidence and the median/tail score models
// disagree. In that regime, try an independent exact KSW2 implementation
// first; WFA remains in Tier 1 as a fallback if KSW2 cannot complete.
bool preferIndependentKsw2ForAmbiguousWfaProfile(
        const AlignmentProfile &profile);

// -bt policy shared by all exact-DP engines. Zero is exact-first. A positive
// value is balanced mode: every modeled admission quantile must be at or
// below the user ceiling. At present that means both P50 and P90; adding a
// higher admission quantile must extend this conjunction rather than bypass
// the ceiling.
bool exactCandidateWithinTimeLimit(
        double estimatedMinutesP50,
        double estimatedMinutesP90,
        double maximumEstimatedMinutes);
bool exactCandidateFastLaneEligible(
        double predictedWaitMinutes,
        double candidateExpectedMinutes,
        double fastestExactExpectedMinutes,
        double memoryOpportunityMinutes);

// Conservative, engine-neutral short-wait comparison. A waiting Tier-1
// candidate may be held ahead of an immediately runnable Tier-1 candidate only
// when its P90 plus wait/pressure still beats the runnable candidate's P50.
bool waitingExactDominatesRunnableExact(
        const AlgorithmCostEstimate &waitingEstimate,
        double waitingMinutes,
        const AlgorithmCostEstimate &runnableEstimate,
        double waitingMemoryShadowMinutes = 0.0,
        double runnableMemoryShadowMinutes = 0.0,
        double maximumWaitMinutes = 2.0);

// Set the process-wide policy used by alignSlidingWindow(). The CLI sets this
// once before worker threads start. Zero selects exact-first and disables all
// exact-DP time limits. A positive value selects balanced mode for every exact
// engine. Temporary memory occupancy never activates the gate.
void configureExactAlignmentMaximumEstimatedMinutes(double minutes);
double configuredExactAlignmentMaximumEstimatedMinutes();

void resetAlignmentSelectionTelemetry();
void recordAlignmentSelectionPlan(const AlignmentSelectionPlan &plan);
void recordExactAlignmentRuntimeFailure(bool memoryFailure);
void recordBandedFallbackExecution();
void recordSlidingFallbackExecution();
AlignmentSelectionTelemetry alignmentSelectionTelemetrySnapshot();

// Optional compact per-attempt TSV trace.  This is intentionally opt-in so a
// production run has no per-interval terminal output.  Reconfiguring with an
// empty path disables tracing; a non-empty path is truncated and receives one
// header followed by thread-safe records.
void configureAlignmentTraceFile(const std::string &path);
bool alignmentTraceEnabled();

struct AlignmentAttemptTrace {
    uint64_t intervalId = 0;
    uint64_t attempt = 0;
    AlignmentProfile profile;
    AlignmentSelectionProvenance provenance;
    ApproximateQualityDecision approximateDecision;
    int64_t approximateKsw2BandWidth = -1;
    int64_t approximateSlidingWindowWidth = -1;
    AlignmentCandidate candidate = AlignmentCandidate::SlidingWindow;
    double predictedMinutesP50 = 0.0;
    double predictedMinutesP90 = 0.0;
    uint64_t predictedMemoryBytes = 0;
    uint64_t reservedMemoryBytes = 0;
    double predictedWaitMinutes = 0.0;
    double actualSeconds = 0.0;
    double configuredExactLimitMinutes = 0.0;
    double exactRuntimeSpentSeconds = 0.0;
    double exactRuntimeRemainingSeconds = 0.0;
    int64_t actualScore = 0;
    int64_t certifiedOptimalScore = 0;
    int64_t certifiedInitialBandWidth = -1;
    int64_t certifiedMaximumBandWidth = -1;
    int64_t certifiedFinalBandWidth = -1;
    uint64_t certifiedTracebackAttempts = 0;
    uint64_t actualMemoryBytes = 0;
    uint64_t processResidentBytes = 0;
    uint64_t processMemoryLimitBytes = 0;
    uint64_t activeReservedBytes = 0;
    int workerCount = 1;
    int activeTasks = 0;
    uint64_t readyTasks = 0;
    uint64_t deferredTasks = 0;
    uint64_t plannedTasks = 0;
    uint64_t schedulableTasks = 0;
    uint64_t globalFutureTasks = 0;
    uint64_t blockedOrderedHeads = 0;
    uint64_t outstandingTasks = 0;
    uint64_t outstandingEstimatedCost = 0;
    uint64_t criticalEstimatedCost = 0;
    bool globalTailPhase = false;
    bool admissionTailPhase = false;
    bool tailPhase = false;
    // Empty for planning/failed attempts.  On completion this is the exact
    // label emitted to the methods BED for the alignment that was retained.
    // It can differ from candidate when a lower-scoring approximate result is
    // discarded in favour of an already-computed retained alignment.
    std::string resultMethod;
    std::string status;
};

uint64_t nextAlignmentTraceIntervalId();
void recordAlignmentAttemptTrace(const AlignmentAttemptTrace &record);

}  // namespace anchorwave
