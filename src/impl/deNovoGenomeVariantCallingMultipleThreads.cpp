//
// Created by Baoxing song on 20.10.18.
//

#include "deNovoGenomeVariantCalling.h"
#include "../myImportandFunction/WfaAlignment.h"
#include "../service/ParallelGapBatchScheduler.h"

std::mutex g_num_mutex;

namespace {

using FastaIndex =
        std::map<std::string, std::tuple<std::string, long, long, int>>;
using FastaEntry = std::tuple<std::string, long, long, int>;

struct ParallelGapResult {
    std::string queryAlignment;
    std::string referenceAlignment;
    std::string method;
    int64_t score = 0;
};

void reportAlignmentSelectionTelemetry() {
    const anchorwave::AlignmentSelectionTelemetry telemetry =
            anchorwave::alignmentSelectionTelemetrySnapshot();
    std::cerr << "AnchorWave alignment selector: evaluated_intervals="
              << telemetry.evaluatedIntervals
              << ", exact_tier_intervals="
              << telemetry.exactTierIntervals
              << ", banded_only_intervals="
              << telemetry.bandedOnlyIntervals
              << ", sliding_only_intervals="
              << telemetry.slidingOnlyIntervals
              << ", singletrack_wfa_memory_rejects="
              << telemetry.singletrackWfaMemoryRejects
              << ", singletrack_wfa_time_rejects="
              << telemetry.singletrackWfaTimeRejects
              << ", wfa_memory_rejects="
              << telemetry.standardWfaMemoryRejects
              << ", wfa_work_warnings="
              << telemetry.standardWfaWorkWarnings
              << ", wfa_time_rejects="
              << telemetry.standardWfaTimeRejects
              << ", wfa_medium_memory_rejects="
              << telemetry.mediumWfaMemoryRejects
              << ", wfa_medium_work_warnings="
              << telemetry.mediumWfaWorkWarnings
              << ", wfa_medium_time_rejects="
              << telemetry.mediumWfaTimeRejects
              << ", wfa_low_memory_rejects="
              << telemetry.lowWfaMemoryRejects
              << ", wfa_low_work_warnings="
              << telemetry.lowWfaWorkWarnings
              << ", wfa_low_time_rejects="
              << telemetry.lowWfaTimeRejects
              << ", score_certified_ksw2_memory_rejects="
              << telemetry.scoreCertifiedKsw2MemoryRejects
              << ", score_certified_ksw2_time_rejects="
              << telemetry.scoreCertifiedKsw2TimeRejects
              << ", full_ksw2_memory_rejects="
              << telemetry.fullKsw2MemoryRejects
              << ", full_ksw2_time_rejects="
              << telemetry.fullKsw2TimeRejects
              << ", banded_ksw2_memory_rejects="
              << telemetry.bandedKsw2MemoryRejects
              << ", exact_runtime_memory_failures="
              << telemetry.exactRuntimeMemoryFailures
              << ", exact_runtime_other_failures="
              << telemetry.exactRuntimeOtherFailures
              << ", banded_fallback_executions="
              << telemetry.bandedFallbackExecutions
              << ", sliding_fallback_executions="
              << telemetry.slidingFallbackExecutions << std::endl;
}

void reportAlignmentResourcePlan(
        const anchorwave::AlignmentResourcePlan &plan) {
    std::cerr << "AnchorWave resource scheduler: requested_threads="
              << plan.requestedMaxThreads
              << ", active_thread_limit=" << plan.effectiveMaxThreads
              << ", per_worker_budget_bytes=" << plan.perWorkerMemoryBytes
              << ", process_memory_limit_bytes="
              << plan.maxProcessMemoryBytes
              << ", baseline_rss_bytes=" << plan.baselineResidentBytes
              << ", in_limit_safety_reserve_bytes="
              << plan.safetyReserveBytes
              << ", task_memory_capacity_bytes="
              << plan.taskMemoryCapacityBytes
              << ", alignment_mode="
              << (anchorwave::configuredExactAlignmentMaximumEstimatedMinutes()
                          == 0.0 ? "exact-first" : "balanced")
              << ", exact_dp_max_estimated_minutes="
              << anchorwave::configuredExactAlignmentMaximumEstimatedMinutes()
              << std::endl;
}

void reportAlignmentMemoryScheduler(
        const anchorwave::AlignmentMemoryScheduler &scheduler) {
    const anchorwave::AlignmentMemorySchedulerStats stats = scheduler.stats();
    std::cerr << "AnchorWave memory scheduler: reservations="
              << stats.reservationCount
              << ", waited_reservations=" << stats.waitedReservationCount
              << ", temporary_reservation_deferrals="
              << stats.temporaryReservationDeferrals
              << ", preferred_reservations="
              << stats.preferredReservationCount
              << ", preferred_waits=" << stats.preferredWaitCount
              << ", impossible_reservations="
              << stats.impossibleReservationCount
              << ", peak_concurrent_reservations="
              << stats.peakConcurrentReservations
              << ", peak_reserved_bytes=" << stats.peakReservedBytes
              << ", peak_projected_process_bytes="
              << stats.peakProjectedProcessBytes
              << ", peak_sampled_rss_bytes="
              << stats.peakObservedResidentBytes << std::endl;
}

std::vector<anchorwave::AlignmentAnchorSpan> makeAlignmentAnchorSpans(
        const std::vector<AlignmentMatch> &matches) {
    std::vector<anchorwave::AlignmentAnchorSpan> spans;
    spans.reserve(matches.size());
    for (const AlignmentMatch &match : matches) {
        anchorwave::AlignmentAnchorSpan span;
        span.referenceStart = match.getRefStartPos();
        span.referenceEnd = match.getRefEndPos();
        span.queryStart = match.getQueryStartPos();
        span.queryEnd = match.getQueryEndPos();
        span.reverse = match.getStrand() == NEGATIVE;
        spans.push_back(span);
    }
    return spans;
}

ParallelGapResult alignParallelGap(
        const FastaEntry &referenceEntry, const FastaEntry &queryEntry,
        const std::string &referenceChromosome,
        const std::string &queryChromosome,
        anchorwave::AlignmentGapDescriptor &gap,
        int32_t windowWidth, int32_t matchingScore,
        int32_t mismatchingPenalty, int32_t openGapPenalty1,
        int32_t extendGapPenalty1, int32_t openGapPenalty2,
        int32_t extendGapPenalty2) {
    // getSubsequence2 currently accepts a mutable map.  Each worker receives
    // private one-entry maps so concurrent lookups cannot mutate shared state.
    FastaIndex referenceIndex;
    referenceIndex.emplace(referenceChromosome, referenceEntry);
    FastaIndex queryIndex;
    queryIndex.emplace(queryChromosome, queryEntry);

    std::shared_ptr<std::string> referenceSequence =
            gap.preparedReferenceSequence;
    std::shared_ptr<std::string> querySequence = gap.preparedQuerySequence;
    if (!referenceSequence) {
        referenceSequence = std::make_shared<std::string>(getSubsequence2(
                referenceIndex, referenceChromosome,
                static_cast<int>(gap.referenceStart),
                static_cast<int>(gap.referenceEnd)));
        gap.preparedReferenceSequence = referenceSequence;
    }
    if (!querySequence) {
        querySequence = std::make_shared<std::string>(getSubsequence2(
                queryIndex, queryChromosome,
                static_cast<int>(gap.queryStart),
                static_cast<int>(gap.queryEnd),
                gap.reverse ? NEGATIVE : POSITIVE));
        gap.preparedQuerySequence = querySequence;
    }

    ParallelGapResult result;
    result.score = alignSlidingWindow(
            *querySequence, *referenceSequence, result.queryAlignment,
            result.referenceAlignment, result.method, windowWidth,
            matchingScore, mismatchingPenalty, openGapPenalty1,
            extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
            gap.selectionPlan);

    if (result.method == "BANDED_MINIMAP2") {
        if (!anchorwave::ungappedSequenceEquals(
                    result.referenceAlignment, *referenceSequence) ||
            !anchorwave::ungappedSequenceEquals(
                    result.queryAlignment, *querySequence)) {
            result.method = "SLIDING_WINDOW";
            result.score = alignSlidingWindowNW(
                    *querySequence, *referenceSequence, result.queryAlignment,
                    result.referenceAlignment, windowWidth, matchingScore,
                    mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                    openGapPenalty2, extendGapPenalty2);
        }
    }
    return result;
}

void prepareGapSchedulingPlans(
        std::vector<anchorwave::AlignmentGapDescriptor> &descriptors,
        const FastaEntry &referenceEntry, const FastaEntry &queryEntry,
        const std::string &referenceChromosome,
        const std::string &queryChromosome,
        int32_t windowWidth, int32_t mismatchingPenalty,
        int32_t openGapPenalty1, int32_t extendGapPenalty1,
        int32_t openGapPenalty2, int32_t extendGapPenalty2) {
    FastaIndex referenceIndex;
    referenceIndex.emplace(referenceChromosome, referenceEntry);
    FastaIndex queryIndex;
    queryIndex.emplace(queryChromosome, queryEntry);
    for (anchorwave::AlignmentGapDescriptor &descriptor : descriptors) {
        std::string referenceSequence = getSubsequence2(
                referenceIndex, referenceChromosome,
                static_cast<int>(descriptor.referenceStart),
                static_cast<int>(descriptor.referenceEnd));
        std::string querySequence = getSubsequence2(
                queryIndex, queryChromosome,
                static_cast<int>(descriptor.queryStart),
                static_cast<int>(descriptor.queryEnd),
                descriptor.reverse ? NEGATIVE : POSITIVE);
        const bool filling = std::all_of(
                referenceSequence.begin(), referenceSequence.end(),
                [](char base) { return base == 'N'; }) ||
                std::all_of(
                        querySequence.begin(), querySequence.end(),
                        [](char base) { return base == 'N'; });
        if (filling) {
            descriptor.predictedMinutes = 0.0;
            descriptor.schedulingPriorityCost =
                    anchorwave::alignmentGapRuntimePriorityCost(0.0);
            continue;
        }
        descriptor.selectionPlan = prepareAlignmentSelectionPlan(
                querySequence, referenceSequence, windowWidth,
                mismatchingPenalty,
                openGapPenalty1, extendGapPenalty1,
                openGapPenalty2, extendGapPenalty2);
        descriptor.predictedMinutes =
                anchorwave::alignmentSelectionPriorityMinutes(
                        *descriptor.selectionPlan);
        descriptor.predictedMemoryBytes =
                anchorwave::alignmentSelectionPriorityMemoryBytes(
                        *descriptor.selectionPlan);
        const uint64_t workerBudget = anchorwave::wfaMemoryBudgetBytes(
                windowWidth);
        const double memoryShare = workerBudget == 0 ? 1.0 : std::min(
                1.0,
                static_cast<double>(descriptor.predictedMemoryBytes) /
                static_cast<double>(workerBudget));
        // LPT still controls the critical path, but expose long, low-memory
        // descriptors early enough to backfill a partially occupied process
        // envelope. B73/Mo17 exact-first runs leave many CPU slots idle while
        // the memory envelope is full; a bounded 150% bonus lets useful work
        // reach the global queue without changing resource admission.
        const double backfillPriorityMinutes = descriptor.predictedMinutes *
                (1.0 + 1.50 * (1.0 - memoryShare));
        descriptor.schedulingPriorityCost =
                anchorwave::alignmentGapRuntimePriorityCost(
                        backfillPriorityMinutes);
    }
}

// Maintains a worker-sized rolling window of alignment tasks/results per
// collinear block. Every result has independent readiness: requesting an early
// gap never waits for a slow unrelated gap in the same initial window. As each
// result is consumed, one later descriptor is submitted, preventing
// chromosome-scale result accumulation outside the -M plan.
class ParallelGapBatchScheduler {
public:
    ParallelGapBatchScheduler(
            const std::vector<AlignmentMatch> &matches,
            const FastaIndex &referenceIndex, const FastaIndex &queryIndex,
            const std::string &referenceChromosome,
            const std::string &queryChromosome,
            int32_t windowWidth, int32_t matchingScore,
            int32_t mismatchingPenalty, int32_t openGapPenalty1,
            int32_t extendGapPenalty1, int32_t openGapPenalty2,
            int32_t extendGapPenalty2,
            anchorwave::AnchorTaskExecutor &executor,
            std::atomic<uint64_t> &submittedTaskCount)
            : referenceEntry_(referenceIndex.at(referenceChromosome)),
              queryEntry_(queryIndex.at(queryChromosome)),
              referenceChromosome_(referenceChromosome),
              queryChromosome_(queryChromosome),
              windowWidth_(windowWidth),
              matchingScore_(matchingScore),
              mismatchingPenalty_(mismatchingPenalty),
              openGapPenalty1_(openGapPenalty1),
              extendGapPenalty1_(extendGapPenalty1),
              openGapPenalty2_(openGapPenalty2),
              extendGapPenalty2_(extendGapPenalty2) {
        std::vector<anchorwave::AlignmentGapDescriptor> descriptors =
                anchorwave::planParallelInterAnchorGaps(
                        makeAlignmentAnchorSpans(matches),
                        executor.workerCount());
        prepareGapSchedulingPlans(
                descriptors, referenceEntry_, queryEntry_,
                referenceChromosome_, queryChromosome_,
                windowWidth_, mismatchingPenalty_,
                openGapPenalty1_, extendGapPenalty1_,
                openGapPenalty2_, extendGapPenalty2_);
        resultWindow_.reset(new GapResultWindow(
                std::move(descriptors),
                // Keep one worker-width runnable. The core may retain at most
                // two worker-widths of submitted-but-unconsumed results (2T),
                // plus one bounded emergency backfill when memory deferrals
                // otherwise leave a CPU idle. A full B73/Mo17 run with a 2T
                // runnable input (4T retained results) accumulated enough
                // ordered output strings to exceed -M despite correct per-
                // attempt reservations; the tighter bound keeps that
                // persistent RSS inside the process memory plan.
                static_cast<std::size_t>(executor.workerCount()), executor,
                [this](anchorwave::AlignmentGapDescriptor &descriptor) {
                    return alignParallelGap(
                            referenceEntry_, queryEntry_,
                            referenceChromosome_, queryChromosome_, descriptor,
                            windowWidth_, matchingScore_, mismatchingPenalty_,
                            openGapPenalty1_, extendGapPenalty1_,
                            openGapPenalty2_, extendGapPenalty2_);
                },
                &submittedTaskCount));
        resultWindow_->start();
    }

    ParallelGapBatchScheduler(const ParallelGapBatchScheduler &) = delete;
    ParallelGapBatchScheduler &operator=(const ParallelGapBatchScheduler &) =
            delete;

    ~ParallelGapBatchScheduler() = default;

    std::shared_ptr<const ParallelGapResult> resultBeforeAnchor(
            std::size_t anchorIndex) {
        return resultWindow_->resultBeforeAnchor(anchorIndex);
    }

private:
    using GapResultWindow =
            anchorwave::ParallelGapBatchSchedulerCore<ParallelGapResult>;

    FastaEntry referenceEntry_;
    FastaEntry queryEntry_;
    std::string referenceChromosome_;
    std::string queryChromosome_;
    int32_t windowWidth_;
    int32_t matchingScore_;
    int32_t mismatchingPenalty_;
    int32_t openGapPenalty1_;
    int32_t extendGapPenalty1_;
    int32_t openGapPenalty2_;
    int32_t extendGapPenalty2_;
    // Keep this member last: it is destroyed first and drains every pending
    // task before the captured sequence/context members above are released.
    std::unique_ptr<GapResultWindow> resultWindow_;
};

}  // namespace

void genomeAlignmentSingleThread(const std::vector<AlignmentMatch> &alignmentMatchs,
                                 const bool outPutMaf, const bool outPutFraged, const bool oMethodBed,
                                 std::ofstream &omaffile, std::ofstream &ofragfile, std::ofstream &oMethodBedfile,
                                 std::map<std::string, std::tuple<std::string, long, long, int> > &map_ref,
                                 std::map<std::string, std::tuple<std::string, long, long, int> > &map_qry,
                                 const int chrWidth, const std::string refFileName, const std::string queryFileName,
                                 const int32_t windowWidth,
                                 const int32_t matchingScore, const int32_t mismatchingPenalty,
                                 const int32_t openGapPenalty1, const int32_t extendGapPenalty1,
                                 const int32_t openGapPenalty2, const int32_t extendGapPenalty2,
                                 anchorwave::AnchorTaskExecutor &alignmentExecutor,
                                 std::atomic<uint64_t> &parallelGapTaskCount) {

    std::string refChr = alignmentMatchs[0].getRefChr();
    std::string queryChr = alignmentMatchs[0].getQueryChr();

    const bool checkResult = true;
    anchorwave::AlignmentBlockBuffer alignmentBuffer(outPutMaf);

    STRAND strand = alignmentMatchs[0].getStrand();

    size_t size_refSequence = getSequenceSizeFromPath2(map_ref[refChr]);
    size_t size_targetSequence = getSequenceSizeFromPath2(map_qry[queryChr]);

    ParallelGapBatchScheduler gapScheduler(
            alignmentMatchs, map_ref, map_qry, refChr, queryChr,
            windowWidth, matchingScore, mismatchingPenalty,
            openGapPenalty1, extendGapPenalty1, openGapPenalty2,
            extendGapPenalty2, alignmentExecutor, parallelGapTaskCount);

    if (POSITIVE == strand) {
        size_t startRef = alignmentMatchs[0].getRefStartPos();
        size_t startQuery;
        size_t endRef;
        size_t endQuery;
        int64_t alignmentScore = 0;
        startQuery = alignmentMatchs[0].getQueryStartPos();
        for (std::size_t anchorIndex = 0;
             anchorIndex < alignmentMatchs.size(); ++anchorIndex) {
            const AlignmentMatch &orthologPair = alignmentMatchs[anchorIndex];
            if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() != startQuery) {
                endQuery = orthologPair.getQueryStartPos() - 1;
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.size(), '-');
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, std::string(), querySeq);
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * querySeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 << "\t" << "FILLING" << "\t" << 0 << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1+ querySeq.size()<< std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryStartPos() == startQuery) {
                endRef = orthologPair.getRefStartPos() - 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.size(), '-');
                std::string _alignment_d = refSeq;
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, refSeq, std::string());
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * refSeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << "FILLING" << "\t" << 0 << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1<< std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryStartPos() == startQuery) {

            } else {
                endRef = orthologPair.getRefStartPos() - 1;
                endQuery = orthologPair.getQueryStartPos() - 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);
                std::string alignMethod;
                {
                    std::string _alignment_q;
                    std::string _alignment_d;
                    std::shared_ptr<const ParallelGapResult> parallelGap =
                            gapScheduler.resultBeforeAnchor(anchorIndex);
                    int64_t thiScore;
                    if (parallelGap) {
                        _alignment_q = parallelGap->queryAlignment;
                        _alignment_d = parallelGap->referenceAlignment;
                        alignMethod = parallelGap->method;
                        thiScore = parallelGap->score;
                    } else {
                        thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    }
                    if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                        if (!anchorwave::ungappedSequenceEquals(
                                    _alignment_d, refSeq) ||
                            !anchorwave::ungappedSequenceEquals(
                                    _alignment_q, querySeq)) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            alignMethod="SLIDING_WINDOW";
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                        }
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_d, refSeq));
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_q, querySeq));
                    }
                    alignmentScore += thiScore;

                    alignmentBuffer.append(
                            _alignment_d, _alignment_q, refSeq, querySeq);

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                    if(oMethodBed){
                        g_num_mutex.lock();
                        oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod)  << "\t" << thiScore << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + querySeq.size() << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            {
                startRef = orthologPair.getRefStartPos();
                startQuery = orthologPair.getQueryStartPos();
                endRef = orthologPair.getRefEndPos();
                endQuery = orthologPair.getQueryEndPos();
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                {
                    std::string _alignment_q;
                    std::string _alignment_d;
                    std::string alignMethod;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                        if (!anchorwave::ungappedSequenceEquals(
                                    _alignment_d, refSeq) ||
                            !anchorwave::ungappedSequenceEquals(
                                    _alignment_q, querySeq)) {
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                            alignMethod="SLIDING_WINDOW";
                        }
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_d, refSeq));
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_q, querySeq));
                    }
                    alignmentScore += thiScore;

                    alignmentBuffer.append(
                            _alignment_d, _alignment_q, refSeq, querySeq);

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << querySeq.size() << "\t+\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                    if(oMethodBed){
                        g_num_mutex.lock();
                        oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + querySeq.size() << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            startRef = orthologPair.getRefEndPos() + 1;
            startQuery = orthologPair.getQueryEndPos() + 1;
        }

        {
            const std::size_t referenceSpan =
                    alignmentMatchs.back().getRefEndPos() -
                    alignmentMatchs.front().getRefStartPos() + 1;
            const std::size_t querySpan =
                    alignmentMatchs.back().getQueryEndPos() -
                    alignmentMatchs.front().getQueryStartPos() + 1;
            assert(alignmentBuffer.referenceLength() == referenceSpan);
            assert(alignmentBuffer.queryLength() == querySpan);
            if (outPutMaf) {
                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getRefStartPos() - 1 << "\t" << std::setw(9) << referenceSpan << "\t+\t" << size_refSequence << "\t";
                alignmentBuffer.writeReference(omaffile);
                omaffile << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getQueryStartPos() - 1 << "\t" << std::setw(9) << querySpan << "\t+\t"
                         << size_targetSequence << "\t";
                alignmentBuffer.writeQuery(omaffile);
                omaffile << std::endl << std::endl;
                g_num_mutex.unlock();
            }
        }
    } else {
        size_t startRef = alignmentMatchs[0].getRefStartPos();
        size_t startQuery=0;
        size_t endRef;
        size_t endQuery = alignmentMatchs[0].getQueryEndPos();
        std::string refChr = alignmentMatchs[0].getRefChr();
        std::string queryChr = alignmentMatchs[0].getQueryChr();

        int64_t alignmentScore = 0;
        for (std::size_t anchorIndex = 0;
             anchorIndex < alignmentMatchs.size(); ++anchorIndex) {
            const AlignmentMatch &orthologPair = alignmentMatchs[anchorIndex];

            if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() != endQuery) {
                startQuery = orthologPair.getQueryEndPos() + 1;
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);
                std::string _alignment_q = querySeq;
                std::string _alignment_d = std::string(querySeq.size(), '-');
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, std::string(), querySeq);
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * querySeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * querySeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    int64_t negative_startQuery = size_targetSequence - querySeq.size() - (startQuery-1);
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1  << "\t" << "FILLING" << "\t" << "0" << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + querySeq.size() << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() != startRef && orthologPair.getQueryEndPos() == endQuery) {
                endRef = orthologPair.getRefStartPos() - 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(refSeq.size(), '-');
                std::string _alignment_d = refSeq;
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, refSeq, std::string());
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * refSeq.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    int64_t negative_startQuery = size_targetSequence -0 - (startQuery-1);
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << 0 << "\t-\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << "FILLING" << "\t" << "0" << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1  << std::endl;
                    g_num_mutex.unlock();
                }
            } else if (orthologPair.getRefStartPos() == startRef && orthologPair.getQueryEndPos() == endQuery) {

            } else {
                endRef = orthologPair.getRefStartPos() - 1;
                startQuery = orthologPair.getQueryEndPos() + 1;
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;
                    std::string alignMethod;
                    std::shared_ptr<const ParallelGapResult> parallelGap =
                            gapScheduler.resultBeforeAnchor(anchorIndex);
                    int64_t thiScore;
                    if (parallelGap) {
                        _alignment_q = parallelGap->queryAlignment;
                        _alignment_d = parallelGap->referenceAlignment;
                        alignMethod = parallelGap->method;
                        thiScore = parallelGap->score;
                    } else {
                        thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, alignMethod, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    }

                    if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                        if (!anchorwave::ungappedSequenceEquals(
                                    _alignment_d, refSeq) ||
                            !anchorwave::ungappedSequenceEquals(
                                    _alignment_q, querySeq)) {
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth,  matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                            alignMethod="SLIDING_WINDOW";
                        }
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_d, refSeq));
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_q, querySeq));
                    }
                    alignmentScore += thiScore;

                    alignmentBuffer.append(
                            _alignment_d, _alignment_q, refSeq, querySeq);

                    if (outPutFraged) {
                        int64_t negative_startQuery = size_targetSequence - querySeq.size() - (startQuery-1);
                        g_num_mutex.lock();
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << size_targetSequence << "\t" << _alignment_q
                                  << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                    if(oMethodBed){
                        g_num_mutex.lock();
                        oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + querySeq.size() << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            {
                startRef = orthologPair.getRefStartPos();
                endQuery = orthologPair.getQueryEndPos();
                endRef = orthologPair.getRefEndPos();
                startQuery = orthologPair.getQueryStartPos();
                std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, strand);
                {
                    std::string _alignment_q;
                    std::string _alignment_d;
                    std::string alignMethod;
//                    std::cout << refSeq << std::endl << querySeq << std::endl << "line 276" << std::endl;
                    int64_t thiScore = alignSlidingWindow(querySeq, refSeq, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

                    if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                        if (!anchorwave::ungappedSequenceEquals(
                                    _alignment_d, refSeq) ||
                            !anchorwave::ungappedSequenceEquals(
                                    _alignment_q, querySeq)) {
//                            std::cout << "align error:" << std::endl << refSeq << std::endl << querySeq << std::endl;
                            thiScore = alignSlidingWindowNW(querySeq, refSeq, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                            alignMethod="SLIDING_WINDOW";
                        }
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_d, refSeq));
                        assert(anchorwave::ungappedSequenceEquals(
                                _alignment_q, querySeq));
                    }
                    alignmentScore += thiScore;

                    alignmentBuffer.append(
                            _alignment_d, _alignment_q, refSeq, querySeq);

                    if (outPutFraged) {
                        g_num_mutex.lock();
                        int64_t negative_startQuery = size_targetSequence - querySeq.size() - (startQuery-1);
                        ofragfile << "a\tscore=" << thiScore << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_refSequence << "\t" << _alignment_d << std::endl
                                  << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << querySeq.size() << "\t-\t" << size_targetSequence << "\t" << _alignment_q << std::endl
                                  << std::endl;
                        g_num_mutex.unlock();
                    }
                    if(oMethodBed){
                        g_num_mutex.lock();
                        oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + querySeq.size() << std::endl;
                        g_num_mutex.unlock();
                    }
                }
            }
            startRef = orthologPair.getRefEndPos() + 1;
            endQuery = orthologPair.getQueryStartPos() - 1;
        }
        {
            const std::size_t referenceSpan =
                    alignmentMatchs.back().getRefEndPos() -
                    alignmentMatchs.front().getRefStartPos() + 1;
            const std::size_t querySpan =
                    alignmentMatchs.front().getQueryEndPos() -
                    alignmentMatchs.back().getQueryStartPos() + 1;
            assert(alignmentBuffer.referenceLength() == referenceSpan);
            assert(alignmentBuffer.queryLength() == querySpan);
            if (outPutMaf) {
                int64_t negative_startQuery = size_targetSequence - querySpan - (alignmentMatchs[alignmentMatchs.size() - 1].getQueryStartPos() - 1);
                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << alignmentMatchs[0].getRefStartPos() - 1 << "\t" << std::setw(9) << referenceSpan << "\t+\t" << size_refSequence << "\t";
                alignmentBuffer.writeReference(omaffile);
                omaffile << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << querySpan
                         << "\t-\t" << size_targetSequence << "\t";
                alignmentBuffer.writeQuery(omaffile);
                omaffile << std::endl << std::endl;
                g_num_mutex.unlock();
            }
        }
    }

}

void genomeAlignment(std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                     const std::string &refFastaFilePath, const std::string &targetFastaFilePath,
                     const int32_t &windowWidth,
                     const std::string &outPutMafFile, const std::string &outPutFragedFile, const std::string &outPutBedFile,
                     const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                     const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1,
                     const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                     const int &maxThread,
                     const uint64_t &maxProcessMemoryBytes) {

    anchorwave::resetAlignmentSelectionTelemetry();

    bool outPutMaf = false;
    bool outPutFraged = false;
    bool oMethodBed = false;

    if (outPutMafFile.size() > 0) {
        outPutMaf = true;
    }
    if (outPutFragedFile.size() > 0) {
        outPutFraged = true;
    }
    if (outPutBedFile.size() > 0) {
        oMethodBed = true;
    }

    std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
    readFastaFile(refFastaFilePath, map_ref);

    std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
    readFastaFile(targetFastaFilePath, map_qry);

    anchorwave::AlignmentResourcePlan resourcePlan;
    try {
        resourcePlan = anchorwave::makeAlignmentResourcePlan(
                maxThread, maxProcessMemoryBytes, windowWidth,
                anchorwave::currentProcessResidentBytes());
    } catch (const std::exception &error) {
        std::cerr << "cannot schedule alignments: " << error.what() << std::endl;
        return;
    }
    anchorwave::AlignmentMemoryScheduler memoryScheduler(resourcePlan);
    anchorwave::ScopedAlignmentMemoryScheduler memorySchedulerScope(
            memoryScheduler);
    anchorwave::AnchorTaskExecutor alignmentExecutor(
            resourcePlan.effectiveMaxThreads);
    std::atomic<uint64_t> parallelGapTaskCount(0);
    reportAlignmentResourcePlan(resourcePlan);

    long unsigned int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refFastaFilePath, delim, elems);
    refFileName = elems.back();
    split(targetFastaFilePath, delim, elems);
    queryFileName = elems.back();

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_ref.begin(); itchr != map_ref.end(); ++itchr) {
        if ((refFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }
    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_qry.begin(); itchr != map_qry.end(); ++itchr) {
        if ((queryFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ofragfile;
    std::ofstream oMethodBedfile;

    if (outPutMaf) {
        omaffile.open(outPutMafFile);
//        omaffile << "##maf version=1" << std::endl;
    }

    if (outPutFraged) {
        ofragfile.open(outPutFragedFile);
//        ofragfile << "##maf version=1" << std::endl;
    }

    if (oMethodBed){
        oMethodBedfile.open(outPutBedFile);
    }

    int32_t size = alignmentMatchsMap.size();

    for (int32_t i = size - 1; i >= 0; --i) {
        std::vector<AlignmentMatch> alignmentMatchs =
                std::move(alignmentMatchsMap[i]);
        if (alignmentMatchs.empty()) {
            continue;
        }
        alignmentExecutor.submit(
                anchorwave::anchorTaskEstimatedCost(alignmentMatchs.size()),
                [alignmentMatchs = std::move(alignmentMatchs), outPutMaf,
                 outPutFraged, oMethodBed,
                 &omaffile, &ofragfile, &oMethodBedfile, &map_ref, &map_qry,
                 chrWidth, refFileName, queryFileName, windowWidth,
                 matchingScore, mismatchingPenalty, openGapPenalty1,
                 extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
                 &alignmentExecutor, &parallelGapTaskCount]() mutable {
                    genomeAlignmentSingleThread(
                            alignmentMatchs, outPutMaf,
                            outPutFraged, oMethodBed, omaffile, ofragfile,
                            oMethodBedfile, map_ref, map_qry, chrWidth,
                            refFileName, queryFileName, windowWidth,
                            matchingScore, mismatchingPenalty,
                            openGapPenalty1, extendGapPenalty1,
                            openGapPenalty2, extendGapPenalty2,
                            alignmentExecutor, parallelGapTaskCount);
                });
    }

    alignmentExecutor.waitForIdle();
    std::cerr << "AnchorWave resource scheduler: peak_active_threads="
              << alignmentExecutor.peakActiveTasks()
              << ", parallel_inter_anchor_tasks="
              << parallelGapTaskCount.load(std::memory_order_relaxed)
              << std::endl;
    reportAlignmentMemoryScheduler(memoryScheduler);
    reportAlignmentSelectionTelemetry();
    if (outPutMaf) {
        omaffile.close();
    }
    if (outPutFraged) {
        ofragfile.close();
    }
}

void genomeAlignmentAndVariantCallingSingleThread(
        std::map<std::string, std::tuple<std::string, long, long, int> > & map_ref,
        std::map<std::string, std::tuple<std::string, long, long, int> > & map_qry,
        const std::vector<AlignmentMatch> & v_am,
        const int & chrWidth, const std::string & refFileName, const std::string & queryFileName,
        const bool & outPutMaf, const bool & outPutFraged, const bool & oMethodBed,
        std::ofstream &omaffile,
        std::ofstream &ofragfile, std::ofstream &oMethodBedfile,
        const int32_t & windowWidth,
        const int32_t & matchingScore, const int32_t & mismatchingPenalty,
        const int32_t & openGapPenalty1, const int32_t & extendGapPenalty1, const int32_t & openGapPenalty2, const int32_t & extendGapPenalty2,
        anchorwave::AnchorTaskExecutor &alignmentExecutor,
        std::atomic<uint64_t> &parallelGapTaskCount) {

    std::string refChr = v_am[0].getRefChr();
    std::string queryChr = v_am[0].getQueryChr();

    size_t startRef = 1;
    size_t startQuery = 1;
    size_t endRef;
    size_t endQuery;
    anchorwave::AlignmentBlockBuffer alignmentBuffer(outPutMaf);

    int64_t alignmentScore = 0;
    STRAND lastStrand = POSITIVE;
    AlignmentMatch lastAlignmentMatch;

    int mafRefStart = 0;
    int mafQueryStart = 0;
    std::string mafStrand = "+";
    bool hasInversion = false;
    bool checkResult = true;

    size_t size_ref_sq = getSequenceSizeFromPath2(map_ref[refChr]);
    size_t size_target_sq = getSequenceSizeFromPath2(map_qry[queryChr]);

    std::string path_ref = std::get<0>(map_ref[refChr]);
    std::string path_qry = std::get<0>(map_qry[queryChr]);

    int fd_ref = open(path_ref.c_str(), O_RDONLY);
    int fd_qry = open(path_qry.c_str(), O_RDONLY);

    ParallelGapBatchScheduler gapScheduler(
            v_am, map_ref, map_qry, refChr, queryChr, windowWidth,
            matchingScore, mismatchingPenalty, openGapPenalty1,
            extendGapPenalty1, openGapPenalty2, extendGapPenalty2,
            alignmentExecutor, parallelGapTaskCount);

    for (std::size_t anchorIndex = 0; anchorIndex < v_am.size();
         ++anchorIndex) {
        const AlignmentMatch &alignmentMatch = v_am[anchorIndex];
        if (alignmentMatch.getStrand() == NEGATIVE) {
            hasInversion = true;
        }

        if (lastStrand == POSITIVE && alignmentMatch.getStrand() == POSITIVE) {
            if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() != startQuery) {
                endQuery = alignmentMatch.getQueryStartPos() - 1;
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

                std::string _alignment_q = seq_qry;
                std::string _alignment_d = std::string(seq_qry.size(), '-');
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, std::string(), seq_qry);
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_qry.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_qry.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }

                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 << "\t" << "FILLING" << "\t" << "0" << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryStartPos() == startQuery) {
                endRef = alignmentMatch.getRefStartPos() - 1;
                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(seq_ref.size(), '-');
                std::string _alignment_d = seq_ref;
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, seq_ref, std::string());
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_ref.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_ref.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << "FILLING" << "\t" << "0" << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1  << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() == startQuery) {

            }
            else {
                endRef = alignmentMatch.getRefStartPos() - 1;
                endQuery = alignmentMatch.getQueryStartPos() - 1;

                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

                std::string _alignment_q;
                std::string _alignment_d;
                std::string alignMethod;
                std::shared_ptr<const ParallelGapResult> parallelGap =
                        gapScheduler.resultBeforeAnchor(anchorIndex);
                int64_t thiScore;
                if (parallelGap) {
                    _alignment_q = parallelGap->queryAlignment;
                    _alignment_d = parallelGap->referenceAlignment;
                    alignMethod = parallelGap->method;
                    thiScore = parallelGap->score;
                } else {
                    thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                }
                if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                    if (!anchorwave::ungappedSequenceEquals(
                                _alignment_d, seq_ref) ||
                        !anchorwave::ungappedSequenceEquals(
                                _alignment_q, seq_qry)) {
//                        std::cout << "seq_ref:" << seq_ref << std::endl;
//                        std::cout << "seq_qry:" << seq_qry << std::endl;
//                        std::cout << "_alignment_d1:" << _alignment_d << std::endl;
//                        std::cout << "_alignment_q1:" << _alignment_q << std::endl;
                        alignMethod = "SLIDING_WINDOW";
                        thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//                        std::cout << "_alignment_d2:" << _alignment_d << std::endl;
//                        std::cout << "_alignment_q2:" << _alignment_q << std::endl;
//                        std::cout  << "line 686\t"  << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << alignMethod << "\t" << thiScore << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;


                    }

                    assert(anchorwave::ungappedSequenceEquals(
                            _alignment_d, seq_ref));
                    assert(anchorwave::ungappedSequenceEquals(
                            _alignment_q, seq_qry));
                }

                alignmentScore += thiScore;
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, seq_ref, seq_qry);

                if (outPutFraged) {
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;
                    g_num_mutex.unlock();
                }
            }
            mafStrand = "+";
        }
        else if (lastStrand == NEGATIVE && alignmentMatch.getStrand() == NEGATIVE && lastAlignmentMatch.getRefEndPos() < alignmentMatch.getRefStartPos() && lastAlignmentMatch.getQueryStartPos() > alignmentMatch.getQueryEndPos() ) {
            if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() != startQuery) {
                endQuery = alignmentMatch.getQueryEndPos() + 1;
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                std::string _alignment_q = seq_qry;
                std::string _alignment_d = std::string(seq_qry.size(), '-');
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, std::string(), seq_qry);
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_qry.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_qry.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }
                alignmentScore += thiScore;

                if (outPutFraged) {
                    int64_t negative_startQuery = size_target_sq - seq_qry.size() - (startQuery-1);
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << seq_qry.size() << "\t-\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1  << "\t" << "FILLING" << "\t" << "0" << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryEndPos() == startQuery) {
                endRef = alignmentMatch.getRefStartPos() - 1;
                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);

                std::string _alignment_q = std::string(seq_ref.size(), '-');
                std::string _alignment_d = seq_ref;
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, seq_ref, std::string());
                int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_ref.size();
                int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_ref.size();
                if (thiScore < thiScore2) {
                    thiScore = thiScore2;
                }

                alignmentScore += thiScore;

                if (outPutFraged) {
                    g_num_mutex.lock();
                    int64_t negative_startQuery = size_target_sq - 0 - (startQuery-1);
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << 0 << "\t-\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << "FILLING" << "\t" << "0" << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1  << std::endl;
                    g_num_mutex.unlock();
                }
            }
            else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() == startQuery) {

            }
            else {
                endRef = alignmentMatch.getRefStartPos() - 1;
                endQuery = alignmentMatch.getQueryEndPos() + 1;

                std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
                std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                std::string _alignment_q;
                std::string _alignment_d;
                std::string alignMethod;
                std::shared_ptr<const ParallelGapResult> parallelGap =
                        gapScheduler.resultBeforeAnchor(anchorIndex);
                int64_t thiScore;
                if (parallelGap) {
                    _alignment_q = parallelGap->queryAlignment;
                    _alignment_d = parallelGap->referenceAlignment;
                    alignMethod = parallelGap->method;
                    thiScore = parallelGap->score;
                } else {
                    thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                }
                if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                    if (!anchorwave::ungappedSequenceEquals(
                                _alignment_d, seq_ref) ||
                        !anchorwave::ungappedSequenceEquals(
                                _alignment_q, seq_qry)) {
//                        std::cout << "seq_ref:" << seq_ref << std::endl;
//                        std::cout << "seq_qry:" << seq_qry << std::endl;
//                        std::cout << "_alignment_d1:" << _alignment_d << std::endl;
//                        std::cout << "_alignment_q1:" << _alignment_q << std::endl;

                        alignMethod = "SLIDING_WINDOW";
                        thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//                        std::cout << "_alignment_d2:" << _alignment_d << std::endl;
//                        std::cout << "_alignment_q2:" << _alignment_q << std::endl;
//                        std::cout << "line 818\t"  << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << alignMethod << "\t" << thiScore << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;

                    }

                    assert(anchorwave::ungappedSequenceEquals(
                            _alignment_d, seq_ref));
                    assert(anchorwave::ungappedSequenceEquals(
                            _alignment_q, seq_qry));
                }

                alignmentScore += thiScore;
                alignmentBuffer.append(
                        _alignment_d, _alignment_q, seq_ref, seq_qry);

                if (outPutFraged) {
                    int64_t negative_startQuery = size_target_sq - seq_qry.size() - (endQuery-1);
                    g_num_mutex.lock();
                    ofragfile << "a\tscore=" << thiScore << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                              << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << negative_startQuery << "\t" << std::setw(9) << seq_qry.size() << "\t-\t" << size_target_sq << "\t" << _alignment_q << std::endl
                              << std::endl;
                    g_num_mutex.unlock();
                }
                if(oMethodBed){
                    g_num_mutex.lock();
                    oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << "-" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;
                    g_num_mutex.unlock();
                }
            }
        }
        else {
            if (outPutMaf && !alignmentBuffer.empty()) {
                const std::size_t refLength =
                        alignmentBuffer.referenceLength();
                const std::size_t queryLength =
                        alignmentBuffer.queryLength();
                int32_t tm = mafQueryStart - queryLength + 1;
                if (lastStrand == POSITIVE) {
                    tm = mafQueryStart;
                }

                int32_t this_tm = tm;
                if (mafStrand == "-") {
                    this_tm = size_target_sq - queryLength - tm;
                }

                g_num_mutex.lock();
                omaffile << "a\tscore=" << alignmentScore << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << mafRefStart << "\t" << std::setw(9) << refLength << "\t+\t" << size_ref_sq << "\t";
                alignmentBuffer.writeReference(omaffile);
                omaffile << std::endl
                         << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << this_tm << "\t" << std::setw(9) << queryLength << "\t" << mafStrand << "\t" << size_target_sq << "\t";
                alignmentBuffer.writeQuery(omaffile);
                omaffile << std::endl << std::endl;
                g_num_mutex.unlock();
            }

            alignmentScore = 0;
            alignmentBuffer.reset();

            mafRefStart = alignmentMatch.getRefStartPos() - 1;
            mafQueryStart = alignmentMatch.getQueryStartPos() - 1;
            if (NEGATIVE == alignmentMatch.getStrand()) {
                mafQueryStart = alignmentMatch.getQueryEndPos() - 1;
            }
        }

        {
            if (POSITIVE == alignmentMatch.getStrand()) {
                mafStrand = "+";
            } else {
                mafStrand = "-";
            }

            startRef = alignmentMatch.getRefStartPos();
            startQuery = alignmentMatch.getQueryStartPos();
            endRef = alignmentMatch.getRefEndPos();
            endQuery = alignmentMatch.getQueryEndPos();

            std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
            std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

            std::string _alignment_q;
            std::string _alignment_d;
            std::string alignMethod;
            int64_t thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
            if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                if (!anchorwave::ungappedSequenceEquals(
                            _alignment_d, seq_ref) ||
                    !anchorwave::ungappedSequenceEquals(
                            _alignment_q, seq_qry)) {
//                    std::cout << "seq_ref:" << seq_ref << std::endl;
//                    std::cout << "seq_qry:" << seq_qry << std::endl;
//                    std::cout << "_alignment_d1:" << _alignment_d << std::endl;
//                    std::cout << "_alignment_q1:" << _alignment_q << std::endl;
                    alignMethod="SLIDING_WINDOW";
                    thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
//                    std::cout << "_alignment_d2:" << _alignment_d << std::endl;
//                    std::cout << "_alignment_q2:" << _alignment_q << std::endl;
//                    std::cout   << "line 936\t" << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << alignMethod << "\t" << thiScore << "\t" << mafStrand << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;

                }

                assert(anchorwave::ungappedSequenceEquals(
                        _alignment_d, seq_ref));
                assert(anchorwave::ungappedSequenceEquals(
                        _alignment_q, seq_qry));
            }

            alignmentScore += thiScore;
            alignmentBuffer.append(
                    _alignment_d, _alignment_q, seq_ref, seq_qry);

            if (outPutFraged) {
                int64_t this_startQuery = (startQuery-1);
                if( mafStrand == "-"){
                    this_startQuery = size_target_sq - seq_qry.size() - (startQuery-1);
                }
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << this_startQuery << "\t" << std::setw(9) << seq_qry.size() << "\t" << mafStrand << "\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
            if(oMethodBed){
                g_num_mutex.lock();
                oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << mafStrand << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 + seq_qry.size() << std::endl;
                g_num_mutex.unlock();
            }
        }

        startRef = alignmentMatch.getRefEndPos() + 1;
        startQuery = alignmentMatch.getQueryEndPos() + 1;
        if (alignmentMatch.getStrand() == NEGATIVE) {
            startQuery = alignmentMatch.getQueryStartPos() - 1;
        }

        lastStrand = alignmentMatch.getStrand();
        lastAlignmentMatch = alignmentMatch;
    }


    // last align
    if (!hasInversion) {
        endRef = getSequenceSizeFromPath2(map_ref[refChr]);
        endQuery = getSequenceSizeFromPath2(map_qry[queryChr]);

        if (startRef > endRef && startQuery <= endQuery) {
            std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

            std::string _alignment_q = seq_qry;
            std::string _alignment_d = std::string(seq_qry.size(), '-');
            alignmentBuffer.append(
                    _alignment_d, _alignment_q, std::string(), seq_qry);
            int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * seq_qry.size();
            int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * seq_qry.size();
            if (thiScore < thiScore2) {
                thiScore = thiScore2;
            }

            alignmentScore += thiScore;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
            if(oMethodBed){
                g_num_mutex.lock();
                oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1  << "\t" << "FILLING" << "\t" << "0" << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1+seq_qry.size() << std::endl;
                g_num_mutex.unlock();
            }
        }
        else if (startRef <= endRef && startQuery > endQuery) {
            std::string refSeq = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);

            std::string _alignment_q = std::string(refSeq.size(), '-');
            std::string _alignment_d = refSeq;
            alignmentBuffer.append(
                    _alignment_d, _alignment_q, refSeq, std::string());

            int64_t thiScore = openGapPenalty1 + extendGapPenalty1 * refSeq.size();
            int64_t thiScore2 = openGapPenalty2 + extendGapPenalty2 * refSeq.size();
            if (thiScore < thiScore2) {
                thiScore = thiScore2;
            }
            alignmentScore += thiScore;

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << refSeq.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << 0 << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
            if(oMethodBed){
                g_num_mutex.lock();
                oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + refSeq.size() << "\t" << "FILLING" << "\t" << "0" << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1  << std::endl;
                g_num_mutex.unlock();
            }
        }
        else if (startRef > endRef && startQuery > endQuery) {

        }
        else {
            std::string seq_ref = getSubsequence3(map_ref, fd_ref, refChr, startRef, endRef);
            std::string seq_qry = getSubsequence3(map_qry, fd_qry, queryChr, startQuery, endQuery);

            std::string _alignment_q;
            std::string _alignment_d;
            std::string alignMethod;
            int64_t thiScore = alignSlidingWindow(seq_qry, seq_ref, _alignment_q, _alignment_d, alignMethod, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
            if (alignMethod.compare("BANDED_MINIMAP2")==0 && checkResult) {
                if (!anchorwave::ungappedSequenceEquals(
                            _alignment_d, seq_ref) ||
                    !anchorwave::ungappedSequenceEquals(
                            _alignment_q, seq_qry)) {
//                    std::cout << "seq_ref:" << seq_ref << std::endl;
//                    std::cout << "seq_qry:" << seq_qry << std::endl;
//                    std::cout << "_alignment_d1:" << _alignment_d << std::endl;
//                    std::cout << "_alignment_q1:" << _alignment_q << std::endl;
                    thiScore = alignSlidingWindowNW(seq_qry, seq_ref, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);
                    alignMethod="SLIDING_WINDOW";
//                    std::cout << "_alignment_d2:" << _alignment_d << std::endl;
//                    std::cout << "_alignment_q2:" << _alignment_q << std::endl;
//                    std::cout << "line 1077\t" << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << alignMethod << "\t" << thiScore << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 +seq_qry.size() << std::endl;

                }
                assert(anchorwave::ungappedSequenceEquals(
                        _alignment_d, seq_ref));
                assert(anchorwave::ungappedSequenceEquals(
                        _alignment_q, seq_qry));
            }

            alignmentScore += thiScore;
            alignmentBuffer.append(
                    _alignment_d, _alignment_q, seq_ref, seq_qry);

            if (outPutFraged) {
                g_num_mutex.lock();
                ofragfile << "a\tscore=" << thiScore << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << startRef - 1 << "\t" << std::setw(9) << seq_ref.size() << "\t+\t" << size_ref_sq << "\t" << _alignment_d << std::endl
                          << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << startQuery - 1 << "\t" << std::setw(9) << seq_qry.size() << "\t+\t" << size_target_sq << "\t" << _alignment_q << std::endl
                          << std::endl;
                g_num_mutex.unlock();
            }
            if(oMethodBed){
                g_num_mutex.lock();
                oMethodBedfile << refChr << "\t" << startRef - 1 << "\t" << startRef - 1 + seq_ref.size() << "\t" << anchorwave::alignmentMethodBedLabel(alignMethod) << "\t" << thiScore << "\t" << "+" << "\t" << queryChr << "\t" << startQuery-1 << "\t" << startQuery-1 +seq_qry.size() << std::endl;
                g_num_mutex.unlock();
            }
        }
        assert(alignmentBuffer.referenceLength() == size_ref_sq);
        assert(alignmentBuffer.queryLength() == size_target_sq);
    }

    if (outPutMaf) {
        const std::size_t refLength = alignmentBuffer.referenceLength();
        const std::size_t queryLength = alignmentBuffer.queryLength();
        int32_t tm = mafQueryStart - queryLength + 1;
        if (lastStrand == POSITIVE) {
            tm = mafQueryStart;
        }
        int32_t this_tm = tm;
        if (mafStrand == "-") {
            this_tm = size_target_sq - queryLength - tm;
        }
        g_num_mutex.lock();
        omaffile << "a\tscore=" << alignmentScore << std::endl
                 << "s\t" << std::left << std::setw(chrWidth) << refChr << "\t" << std::right << std::setw(9) << mafRefStart << "\t" << std::setw(9) << refLength << "\t+\t" << size_ref_sq << "\t";
        alignmentBuffer.writeReference(omaffile);
        omaffile << std::endl
                 << "s\t" << std::left << std::setw(chrWidth) << queryChr << "\t" << std::right << std::setw(9) << this_tm << "\t" << std::setw(9) << queryLength << "\t" + mafStrand + "\t" << size_target_sq << "\t";
        alignmentBuffer.writeQuery(omaffile);
        omaffile << std::endl << std::endl;
        g_num_mutex.unlock();
    }

    close(fd_ref);
    close(fd_qry);

}


void genomeAlignmentAndVariantCalling(std::map<std::string, std::vector<AlignmentMatch>> &map_v_am,
                                      const std::string &path_ref_GenomeSequence, const std::string &path_target_GenomeSequence,
                                      const int32_t &windowWidth, const std::string &outPutMafFile,
                                      const std::string &outPutFragedFile, const std::string &outPutBedFile, const int32_t &matchingScore, const int32_t &mismatchingPenalty,
                                      const int32_t &openGapPenalty1, const int32_t &extendGapPenalty1, const int32_t &openGapPenalty2, const int32_t &extendGapPenalty2,
                                      const int &maxThread,
                                      const uint64_t &maxProcessMemoryBytes) {

    anchorwave::resetAlignmentSelectionTelemetry();

    bool outPutMaf = false;
    bool outPutFraged = false;
    bool oMethodBed = false;
    if (outPutMafFile.size() > 0) {
        outPutMaf = true;
    }

    if (outPutFragedFile.size() > 0) {
        outPutFraged = true;
    }
    if (outPutBedFile.size() > 0) {
        oMethodBed = true;
    }
    std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
    readFastaFile(path_ref_GenomeSequence, map_ref);

    std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
    readFastaFile(path_target_GenomeSequence, map_qry);

    anchorwave::AlignmentResourcePlan resourcePlan;
    try {
        resourcePlan = anchorwave::makeAlignmentResourcePlan(
                maxThread, maxProcessMemoryBytes, windowWidth,
                anchorwave::currentProcessResidentBytes());
    } catch (const std::exception &error) {
        std::cerr << "cannot schedule alignments: " << error.what() << std::endl;
        return;
    }
    anchorwave::AlignmentMemoryScheduler memoryScheduler(resourcePlan);
    anchorwave::ScopedAlignmentMemoryScheduler memorySchedulerScope(
            memoryScheduler);
    anchorwave::AnchorTaskExecutor alignmentExecutor(
            resourcePlan.effectiveMaxThreads);
    std::atomic<uint64_t> parallelGapTaskCount(0);
    reportAlignmentResourcePlan(resourcePlan);

    long unsigned int chrWidth = 4;
    std::string refFileName;
    std::string queryFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(path_ref_GenomeSequence, delim, elems);
    refFileName = elems.back();
    split(path_target_GenomeSequence, delim, elems);
    queryFileName = elems.back();

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_ref.begin(); itchr != map_ref.end(); ++itchr) {
        if ((refFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (refFileName + "." + itchr->first).size();
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator itchr = map_qry.begin(); itchr != map_qry.end(); ++itchr) {
        if ((queryFileName + "." + itchr->first).size() > chrWidth) {
            chrWidth = (queryFileName + "." + itchr->first).size();
        }
    }

    std::ofstream omaffile;
    std::ofstream ofragfile;
    if (outPutMaf) {
        omaffile.open(outPutMafFile);
//        omaffile << "##maf version=1" << std::endl;
    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon));
    if (ltm->tm_mday < 10) {
        filedate = filedate + "0" + std::to_string((ltm->tm_mday));
    } else {
        filedate = filedate + std::to_string((ltm->tm_mday));
    }

    if (outPutFraged) {
        ofragfile.open(outPutFragedFile);
//        ofragfile << "##maf version=1" << std::endl;
    }
    std::ofstream oMethodBedfile;

    if (oMethodBed){
        oMethodBedfile.open(outPutBedFile);
        oMethodBedfile << "# FILLING: no alignment was performed" << std::endl;
        oMethodBedfile << "# WAVEFRONT: exact 2-piece affine alignment was performed using a WFA mode (Singletrack/high/medium/low)" << std::endl;
        oMethodBedfile << "# MINIMAP2: exact 2-piece affine alignment was performed using the ksw_extd2 approach implemented in minimap2, without setting band" << std::endl;
        oMethodBedfile << "# BANDED_MINIMAP2: alignment was performed using the ksw_extd2 approach implemented in minimap2, with band setting" << std::endl;
        oMethodBedfile << "# SLIDING_WINDOW: alignment was performed using a sliding window approach" << std::endl;
        oMethodBedfile << "# WAVEFRONT and MINIMAP2 are exact dynamic-programming alignments. BANDED_MINIMAP2 and SLIDING_WINDOW share the fallback-quality tier; the selector chooses the result predicted to have the higher alignment score" << std::endl;
    }

    for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
        if (it->second.size() > 0) {
            std::vector<AlignmentMatch> *matches = &it->second;
            alignmentExecutor.submit(
                    anchorwave::anchorTaskEstimatedCost(matches->size()),
                    [&map_ref, &map_qry, matches, &chrWidth, &refFileName,
                     &queryFileName, &outPutMaf, &outPutFraged, &oMethodBed,
                     &omaffile, &ofragfile, &oMethodBedfile, &windowWidth,
                     &matchingScore, &mismatchingPenalty, &openGapPenalty1,
                     &extendGapPenalty1, &openGapPenalty2, &extendGapPenalty2,
                     &alignmentExecutor, &parallelGapTaskCount]() {
                        genomeAlignmentAndVariantCallingSingleThread(
                                map_ref, map_qry, *matches, chrWidth,
                                refFileName, queryFileName, outPutMaf,
                                outPutFraged, oMethodBed, omaffile, ofragfile,
                                oMethodBedfile, windowWidth, matchingScore,
                                mismatchingPenalty, openGapPenalty1,
                                extendGapPenalty1, openGapPenalty2,
                                extendGapPenalty2, alignmentExecutor,
                                parallelGapTaskCount);
                    });
        }
    }

    alignmentExecutor.waitForIdle();
    std::cerr << "AnchorWave resource scheduler: peak_active_threads="
              << alignmentExecutor.peakActiveTasks()
              << ", parallel_inter_anchor_tasks="
              << parallelGapTaskCount.load(std::memory_order_relaxed)
              << std::endl;
    reportAlignmentMemoryScheduler(memoryScheduler);
    reportAlignmentSelectionTelemetry();

    if (outPutMaf) {
        omaffile.close();
    }

    if (outPutFraged) {
        ofragfile.close();
    }
}
