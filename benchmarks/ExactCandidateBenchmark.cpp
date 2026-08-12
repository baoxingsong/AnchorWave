#include "src/myImportandFunction/AlignmentAlgorithmSelector.h"
#include "src/myImportandFunction/WfaAlignment.h"
#include "src/myImportandFunction/alignSlidingWindow.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int32_t alignment_minimap2(
        const std::string &query, const std::string &reference,
        std::string &queryAlignment, std::string &referenceAlignment,
        const int32_t &matchingScore, int32_t mismatchingPenalty,
        int32_t openGapPenalty1, int32_t extendGapPenalty1,
        int32_t openGapPenalty2, int32_t extendGapPenalty2,
        int32_t &endPositionQuery, int32_t &endPositionReference);

int32_t minimap2_alignment(
        const std::string &query, const std::string &reference,
        std::string &queryAlignment, std::string &referenceAlignment,
        const int32_t &matchingScore, int32_t mismatchingPenalty,
        int32_t openGapPenalty1, int32_t extendGapPenalty1,
        int32_t openGapPenalty2, int32_t extendGapPenalty2);

namespace {

struct SequencePair {
    std::string query;
    std::string reference;
};

struct Options {
    std::string referenceFasta;
    std::string queryFasta;
    std::string mafPath;
    std::size_t mafBlock = 0;
    std::string algorithm;
    uint64_t memoryBudgetBytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
    int wfaThreads = 1;
    int minOffsetsPerThread = 500;
    int biWfaLeafScore = 0;
    int memoryProbeScoreInterval = 3000;
    int windowWidth = 100000;
};

std::string removeGaps(std::string sequence) {
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '-'),
                   sequence.end());
    return sequence;
}

std::string alignmentHash(const std::string &query,
                          const std::string &reference) {
    uint64_t hash = 1469598103934665603ULL;
    const auto update = [&hash](unsigned char value) {
        hash ^= static_cast<uint64_t>(value);
        hash *= 1099511628211ULL;
    };
    for (const unsigned char value : query) {
        update(value);
    }
    update(0xff);
    for (const unsigned char value : reference) {
        update(value);
    }
    std::ostringstream stream;
    stream << std::hex << std::setw(16) << std::setfill('0') << hash;
    return stream.str();
}

bool readSingleFasta(const std::string &path, std::string &sequence) {
    std::ifstream input(path);
    if (!input.good()) {
        return false;
    }
    sequence.clear();
    std::string line;
    bool sawHeader = false;
    bool sawSequence = false;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            if (sawHeader) {
                throw std::runtime_error(
                        "benchmark FASTA must contain exactly one record: " +
                        path);
            }
            sawHeader = true;
            continue;
        }
        if (!sawHeader) {
            throw std::runtime_error("FASTA sequence precedes header: " + path);
        }
        for (const char base : line) {
            if (base != ' ' && base != '\t') {
                sequence.push_back(base);
            }
        }
        sawSequence = true;
    }
    return sawHeader && sawSequence;
}

bool readMafBlock(const std::string &path,
                  std::size_t requestedBlock,
                  SequencePair &pair) {
    std::ifstream input(path);
    if (!input.good() || requestedBlock == 0) {
        return false;
    }

    std::size_t block = 0;
    std::vector<std::string> alignedRows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.size() >= 2 && line[0] == 'a' &&
            (line[1] == ' ' || line[1] == '\t')) {
            ++block;
            alignedRows.clear();
            continue;
        }
        if (block != requestedBlock || line.size() < 2 || line[0] != 's' ||
            (line[1] != ' ' && line[1] != '\t')) {
            continue;
        }
        std::istringstream fields(line);
        std::string record;
        std::string source;
        std::string start;
        std::string size;
        std::string strand;
        std::string sourceSize;
        std::string alignedSequence;
        if (!(fields >> record >> source >> start >> size >> strand >>
              sourceSize >> alignedSequence)) {
            return false;
        }
        alignedRows.push_back(std::move(alignedSequence));
        if (alignedRows.size() == 2) {
            pair.reference = removeGaps(std::move(alignedRows[0]));
            pair.query = removeGaps(std::move(alignedRows[1]));
            return true;
        }
    }
    return false;
}

const char *statusName(anchorwave::WfaAlignmentStatus status) {
    switch (status) {
        case anchorwave::WfaAlignmentStatus::Completed: return "completed";
        case anchorwave::WfaAlignmentStatus::MemoryLimit:
            return "memory_limit";
        case anchorwave::WfaAlignmentStatus::TimeLimit:
            return "time_limit";
        case anchorwave::WfaAlignmentStatus::Failed: return "failed";
    }
    return "unknown";
}

using WfaAligner = anchorwave::WfaAlignmentStatus (*)(
        const std::string &,
        const std::string &,
        uint64_t,
        int32_t,
        int32_t,
        int32_t,
        int32_t,
        int32_t,
        anchorwave::WfaAlignmentResult &,
        const anchorwave::WfaExecutionOptions &);

void printHeader() {
    std::cout
            << "candidate\tstatus\tscore\tmemory_bytes\tmemory_peak_bytes\t"
            << "seconds\tsequence_reconstruction\tquery_length\t"
            << "reference_length\talignment_hash_fnv1a64\twfa_status\t"
            << "configured_threads\t"
            << "actual_threads\twfa_parallel_support\t"
            << "configured_biwfa_leaf_score\tbiwfa_leaf_alignments\t"
            << "biwfa_max_leaf_score\testimated_score\tconservative_score\t"
            << "profile_uncertainty\tsketch_jaccard\tlength_ratio\t"
            << "predicted_memory_bytes\tpredicted_minutes_p50\t"
            << "predicted_minutes_p90\tpredicted_success_probability\t"
            << "predicted_score_loss\tpredicted_score_loss_p90\t"
            << "predicted_memory_feasible\tpredicted_time_feasible\t"
            << "certified_optimal_score\tcertified_initial_band\t"
            << "certified_maximum_band\tcertified_final_band\t"
            << "certified_traceback_attempts\n";
}

const anchorwave::AlgorithmCostEstimate *findEstimate(
        const anchorwave::AlignmentSelectionPlan &plan,
        anchorwave::AlignmentCandidate candidate) {
    for (const anchorwave::AlgorithmCostEstimate &estimate : plan.estimates) {
        if (estimate.candidate == candidate) {
            return &estimate;
        }
    }
    return nullptr;
}

void printPrediction(
        const anchorwave::AlignmentSelectionPlan &plan,
        const anchorwave::AlgorithmCostEstimate *estimate) {
    std::cout << '\t' << plan.profile.estimatedScore
              << '\t' << plan.profile.conservativeScore
              << '\t' << std::setprecision(9) << plan.profile.uncertainty
              << '\t' << plan.profile.sketchJaccard
              << '\t' << plan.profile.lengthRatio;
    if (estimate == nullptr) {
        std::cout << "\t\t\t\t\t\t\t\t";
        return;
    }
    std::cout << '\t' << estimate->memoryBytes
              << '\t' << estimate->estimatedMinutes
              << '\t' << estimate->estimatedMinutesP90
              << '\t' << estimate->successProbability
              << '\t' << estimate->predictedScoreLoss
              << '\t' << estimate->predictedScoreLossP90
              << '\t' << (estimate->memoryFeasible ? "yes" : "no")
              << '\t' << (estimate->timeFeasible ? "yes" : "no");
}

void printPrediction(const anchorwave::AlignmentSelectionPlan &plan,
                     anchorwave::AlignmentCandidate candidate) {
    printPrediction(plan, findEstimate(plan, candidate));
}

void runWfa(const char *name,
            WfaAligner aligner,
            const SequencePair &pair,
            const Options &options,
            const anchorwave::AlignmentSelectionPlan &plan,
            const anchorwave::AlgorithmCostEstimate *estimate) {
    anchorwave::WfaExecutionOptions execution;
    execution.maxNumThreads = options.wfaThreads;
    execution.minOffsetsPerThread = options.minOffsetsPerThread;
    execution.biWfaLeafScore = options.biWfaLeafScore;
    execution.memoryProbeScoreInterval = options.memoryProbeScoreInterval;

    anchorwave::WfaAlignmentResult result;
    anchorwave::WfaAlignmentStatus status = anchorwave::WfaAlignmentStatus::Failed;
    const auto begin = std::chrono::steady_clock::now();
    try {
        status = aligner(pair.query, pair.reference, options.memoryBudgetBytes,
                         -6, -8, -2, -75, -1, result, execution);
    } catch (const std::bad_alloc &) {
        status = anchorwave::WfaAlignmentStatus::MemoryLimit;
    } catch (const std::exception &error) {
        std::cerr << "candidate " << name << " failed: " << error.what()
                  << '\n';
        status = anchorwave::WfaAlignmentStatus::Failed;
    }
    const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
    const bool reconstructed =
            status == anchorwave::WfaAlignmentStatus::Completed &&
            removeGaps(result.queryAlignment) == pair.query &&
            removeGaps(result.referenceAlignment) == pair.reference;
    std::cout << name << '\t' << statusName(status) << '\t'
              << result.score << '\t' << result.memoryUsedBytes << '\t'
              << result.memoryPeakBytes << '\t' << std::fixed
              << std::setprecision(6) << seconds << '\t'
              << (reconstructed ? "yes" : "no") << '\t'
              << pair.query.size() << '\t' << pair.reference.size() << '\t'
              << alignmentHash(result.queryAlignment,
                               result.referenceAlignment) << '\t'
              << result.wfaStatus << '\t' << result.configuredMaxThreads
              << '\t' << result.maxThreadsUsed << '\t'
              << (anchorwave::wfaParallelSupportEnabled() ? "yes" : "no")
              << '\t' << result.configuredBiWfaLeafScore << '\t'
              << result.biWfaLeafAlignments << '\t'
              << result.biWfaMaxLeafScore;
    printPrediction(plan, estimate);
    std::cout << "\t\t\t\t\t\n";
}

void runKsw2Full(const SequencePair &pair,
                 const anchorwave::AlignmentSelectionPlan &plan) {
    std::string queryAlignment;
    std::string referenceAlignment;
    int64_t score = 0;
    std::string status = "completed";
    const auto begin = std::chrono::steady_clock::now();
    try {
        score = alignSlidingWindow_minimap2(
                pair.query, pair.reference,
                static_cast<int64_t>(pair.query.size()),
                static_cast<int64_t>(pair.reference.size()), queryAlignment,
                referenceAlignment, -1, -6, -8, -2, -75, -1);
    } catch (const std::bad_alloc &) {
        status = "memory_limit";
    } catch (const std::exception &error) {
        std::cerr << "candidate KSW2_FULL failed: " << error.what() << '\n';
        status = "failed";
    }
    const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
    const bool reconstructed = status == "completed" &&
            removeGaps(queryAlignment) == pair.query &&
            removeGaps(referenceAlignment) == pair.reference;
    std::cout << "KSW2_FULL\t" << status << '\t' << score
              << "\t0\t0\t" << std::fixed << std::setprecision(6)
              << seconds << '\t' << (reconstructed ? "yes" : "no")
              << '\t' << pair.query.size() << '\t' << pair.reference.size()
              << '\t' << alignmentHash(queryAlignment, referenceAlignment)
              << "\t0\t1\t1\tno\t0\t0\t0";
    printPrediction(plan, anchorwave::AlignmentCandidate::Ksw2Full);
    std::cout << "\t\t\t\t\t\n";
}

void runApproximate(const SequencePair &pair,
                    const anchorwave::AlignmentSelectionPlan &plan,
                    bool banded) {
    std::string queryAlignment;
    std::string referenceAlignment;
    int64_t score = 0;
    std::string status = "completed";
    const auto candidate = banded
            ? anchorwave::AlignmentCandidate::Ksw2Banded
            : anchorwave::AlignmentCandidate::SlidingWindow;
    const auto begin = std::chrono::steady_clock::now();
    try {
        if (banded) {
            if (plan.ksw2BandWidth <= 0) {
                throw std::runtime_error("selector found no reachable KSW2 band");
            }
            score = alignSlidingWindow_minimap2(
                    pair.query, pair.reference,
                    static_cast<int64_t>(pair.query.size()),
                    static_cast<int64_t>(pair.reference.size()),
                    queryAlignment, referenceAlignment,
                    plan.ksw2BandWidth, -6, -8, -2, -75, -1);
        } else {
            std::string mutableQuery = pair.query;
            std::string mutableReference = pair.reference;
            score = alignSlidingWindowNW(
                    mutableQuery, mutableReference,
                    queryAlignment, referenceAlignment,
                    plan.slidingWindowWidth, 1,
                    -6, -8, -2, -75, -1);
        }
    } catch (const std::bad_alloc &) {
        status = "memory_limit";
    } catch (const std::exception &error) {
        std::cerr << "candidate "
                  << (banded ? "KSW2_BANDED" : "SLIDING_WINDOW")
                  << " failed: " << error.what() << '\n';
        status = "failed";
    }
    const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
    const bool reconstructed = status == "completed" &&
            removeGaps(queryAlignment) == pair.query &&
            removeGaps(referenceAlignment) == pair.reference;
    std::cout << (banded ? "KSW2_BANDED" : "SLIDING_WINDOW")
              << '\t' << status << '\t' << score
              << "\t0\t0\t" << std::fixed << std::setprecision(6)
              << seconds << '\t' << (reconstructed ? "yes" : "no")
              << '\t' << pair.query.size() << '\t' << pair.reference.size()
              << '\t' << alignmentHash(queryAlignment, referenceAlignment)
              << "\t0\t1\t1\tno\t0\t0\t0";
    printPrediction(plan, candidate);
    std::cout << "\t\t\t\t\t\n";
}

void runKsw2ScoreCertified(
        const SequencePair &pair,
        const anchorwave::AlignmentSelectionPlan &plan) {
    std::string queryAlignment;
    std::string referenceAlignment;
    anchorwave::Ksw2CertifiedResult result;
    const auto begin = std::chrono::steady_clock::now();
    try {
        result = alignScoreCertifiedKsw2(
                pair.query, pair.reference,
                queryAlignment, referenceAlignment,
                plan.ksw2CertifiedInitialBandWidth,
                plan.ksw2CertifiedMaximumBandWidth,
                0.0, -6, -8, -2, -75, -1);
    } catch (const std::bad_alloc &) {
        result.status = anchorwave::Ksw2CertifiedStatus::Failed;
    }
    const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
    const char *status = result.status ==
            anchorwave::Ksw2CertifiedStatus::Completed
            ? "completed"
            : result.status == anchorwave::Ksw2CertifiedStatus::TimeLimit
              ? "time_limit"
              : result.status ==
                        anchorwave::Ksw2CertifiedStatus::NotCertified
                ? "not_certified" : "failed";
    const bool reconstructed = result.status ==
            anchorwave::Ksw2CertifiedStatus::Completed &&
            removeGaps(queryAlignment) == pair.query &&
            removeGaps(referenceAlignment) == pair.reference;
    std::cout << "KSW2_SCORE_CERTIFIED\t" << status << '\t'
              << result.score << "\t0\t0\t" << std::fixed
              << std::setprecision(6) << seconds << '\t'
              << (reconstructed ? "yes" : "no") << '\t'
              << pair.query.size() << '\t' << pair.reference.size() << '\t'
              << alignmentHash(queryAlignment, referenceAlignment)
              << "\t0\t1\t1\tno\t0\t0\t0";
    printPrediction(plan,
                    anchorwave::AlignmentCandidate::Ksw2ScoreCertified);
    std::cout << '\t' << result.optimalScore
              << '\t' << plan.ksw2CertifiedInitialBandWidth
              << '\t' << plan.ksw2CertifiedMaximumBandWidth
              << '\t' << result.finalBandWidth
              << '\t' << result.tracebackAttempts << '\n';
}

int32_t legacyTwoPassSemiglobal(
        const std::string &query, const std::string &reference,
        std::string &queryAlignment, std::string &referenceAlignment) {
    int8_t matrix[25];
    for (int row = 0; row < 5; ++row) {
        for (int column = 0; column < 5; ++column) {
            matrix[row * 5 + column] = row == column ? 0 : -6;
        }
    }
    uint8_t encoding[256];
    std::memset(encoding, 4, sizeof(encoding));
    encoding['A'] = encoding['a'] = 0;
    encoding['C'] = encoding['c'] = 1;
    encoding['G'] = encoding['g'] = 2;
    encoding['T'] = encoding['t'] = 3;
    std::vector<uint8_t> encodedQuery(query.size());
    std::vector<uint8_t> encodedReference(reference.size());
    int flags = KSW_EZ_SCORE_ONLY;
    for (std::size_t index = 0; index < query.size(); ++index) {
        encodedQuery[index] = encoding[static_cast<uint8_t>(query[index])];
        if (encodedQuery[index] == 4) flags |= KSW_EZ_GENERIC_SC;
    }
    for (std::size_t index = 0; index < reference.size(); ++index) {
        encodedReference[index] =
                encoding[static_cast<uint8_t>(reference[index])];
        if (encodedReference[index] == 4) flags |= KSW_EZ_GENERIC_SC;
    }
    ksw_extz_t result;
    std::memset(&result, 0, sizeof(result));
#if defined(__AVX512BW__)
    ksw_extd2_avx512(
#elif defined(__AVX2__)
    ksw_extd2_avx2(
#else
    ksw_extd2_sse(
#endif
            nullptr, static_cast<int>(encodedQuery.size()),
            encodedQuery.data(),
            static_cast<int>(encodedReference.size()),
            encodedReference.data(), 5, matrix,
            8, 2, 75, 1, -1, -1, 0, flags, &result);
    const bool queryEndWins = result.mqe > result.mte;
    const int32_t queryEnd = queryEndWins
            ? static_cast<int32_t>(query.size()) : result.mte_q + 1;
    const int32_t referenceEnd = queryEndWins
            ? result.mqe_t + 1 : static_cast<int32_t>(reference.size());
    const int32_t score = queryEndWins ? result.mqe : result.mte;
    std::free(result.cigar);
    const int32_t tracebackScore = minimap2_alignment(
            query.substr(0, static_cast<std::size_t>(queryEnd)),
            reference.substr(0, static_cast<std::size_t>(referenceEnd)),
            queryAlignment, referenceAlignment,
            1, -6, -8, -2, -75, -1);
    if (tracebackScore != score) {
        throw std::runtime_error(
                "legacy semiglobal score/traceback mismatch");
    }
    return score;
}

void runSlidingBlock(const SequencePair &pair, bool legacy) {
    std::string queryAlignment;
    std::string referenceAlignment;
    int32_t queryEnd = 0;
    int32_t referenceEnd = 0;
    int32_t score = 0;
    const auto begin = std::chrono::steady_clock::now();
    if (legacy) {
        score = legacyTwoPassSemiglobal(
                pair.query, pair.reference,
                queryAlignment, referenceAlignment);
    } else {
        score = alignment_minimap2(
                pair.query, pair.reference,
                queryAlignment, referenceAlignment,
                1, -6, -8, -2, -75, -1,
                queryEnd, referenceEnd);
    }
    const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
    const bool reconstructed =
            removeGaps(queryAlignment) ==
                    pair.query.substr(0, removeGaps(queryAlignment).size()) &&
            removeGaps(referenceAlignment) ==
                    pair.reference.substr(
                            0, removeGaps(referenceAlignment).size());
    std::cout << (legacy ? "SLIDING_BLOCK_TWO_PASS" :
                            "SLIDING_BLOCK_ONE_PASS")
              << "\tcompleted\t" << score << "\t0\t0\t"
              << std::fixed << std::setprecision(6) << seconds << '\t'
              << (reconstructed ? "yes" : "no") << '\t'
              << pair.query.size() << '\t' << pair.reference.size() << '\t'
              << alignmentHash(queryAlignment, referenceAlignment)
              << "\t0\t1\t1\tno\t0\t0\t0"
              << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n";
}

void usage(std::ostream &stream) {
    stream
            << "usage:\n"
            << "  anchorwave_exact_candidate_benchmark --reference REF.fa "
            << "--query QUERY.fa --algorithm NAME [options]\n"
            << "  anchorwave_exact_candidate_benchmark --maf FILE "
            << "--maf-block N --algorithm NAME [options]\n"
            << "  anchorwave_exact_candidate_benchmark FRAGMENT_MAF BLOCK "
            << "[MEMORY_BUDGET_BYTES]  # legacy: runs every candidate\n\n"
            << "NAME: ksw2-score-certified, ksw2-full, singletrack, wfa-high, wfa-medium, "
            << "wfa-low, biwfa, ksw2-banded, sliding-window, sliding-block-one-pass, "
            << "sliding-block-two-pass\n"
            << "options:\n"
            << "  --memory-bytes N              per-attempt WFA limit\n"
            << "  --wfa-threads N               inner WFA threads (default 1)\n"
            << "  --min-offsets-per-thread N    WFA parallel grain (default 500)\n"
            << "  --biwfa-leaf-score N          0=derive from memory, 250=legacy\n"
            << "  --memory-probe-score-interval N\n"
            << "  --window N                    selector/sliding window (default 100000)\n";
}

uint64_t parseUnsigned(const std::string &value, const char *name) {
    if (value.empty() || value.front() == '-') {
        throw std::invalid_argument(std::string("invalid ") + name + ": " +
                                    value);
    }
    std::size_t consumed = 0;
    const uint64_t result = std::stoull(value, &consumed);
    if (consumed != value.size()) {
        throw std::invalid_argument(std::string("invalid ") + name + ": " +
                                    value);
    }
    return result;
}

int parsePositiveInt(const std::string &value, const char *name,
                     bool allowZero) {
    const uint64_t parsed = parseUnsigned(value, name);
    if ((!allowZero && parsed == 0) ||
        parsed > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
        throw std::out_of_range(std::string("invalid ") + name + ": " +
                                value);
    }
    return static_cast<int>(parsed);
}

Options parseOptions(int argc, char **argv) {
    Options options;
    if (argc >= 3 && argv[1][0] != '-') {
        if (argc > 4) {
            throw std::invalid_argument("too many legacy arguments");
        }
        options.mafPath = argv[1];
        options.mafBlock = static_cast<std::size_t>(
                parseUnsigned(argv[2], "MAF block"));
        options.algorithm = "all";
        if (argc == 4) {
            options.memoryBudgetBytes = parseUnsigned(argv[3], "memory bytes");
        }
        return options;
    }

    for (int index = 1; index < argc; ++index) {
        const std::string flag = argv[index];
        if (flag == "--help" || flag == "-h") {
            usage(std::cout);
            std::exit(EXIT_SUCCESS);
        }
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + flag);
        }
        const std::string value = argv[++index];
        if (flag == "--reference") {
            options.referenceFasta = value;
        } else if (flag == "--query") {
            options.queryFasta = value;
        } else if (flag == "--maf") {
            options.mafPath = value;
        } else if (flag == "--maf-block") {
            options.mafBlock = static_cast<std::size_t>(
                    parseUnsigned(value, "MAF block"));
        } else if (flag == "--algorithm") {
            options.algorithm = value;
        } else if (flag == "--memory-bytes") {
            options.memoryBudgetBytes = parseUnsigned(value, "memory bytes");
        } else if (flag == "--wfa-threads") {
            options.wfaThreads = parsePositiveInt(value, "WFA threads", false);
        } else if (flag == "--min-offsets-per-thread") {
            options.minOffsetsPerThread =
                    parsePositiveInt(value, "minimum offsets per thread", false);
        } else if (flag == "--biwfa-leaf-score") {
            options.biWfaLeafScore =
                    parsePositiveInt(value, "BiWFA leaf score", true);
        } else if (flag == "--memory-probe-score-interval") {
            options.memoryProbeScoreInterval =
                    parsePositiveInt(value, "memory probe interval", false);
        } else if (flag == "--window") {
            options.windowWidth = parsePositiveInt(value, "window", false);
        } else {
            throw std::invalid_argument("unknown option: " + flag);
        }
    }
    const bool usesMaf = !options.mafPath.empty() && options.mafBlock > 0;
    const bool usesFasta = !options.referenceFasta.empty() &&
                           !options.queryFasta.empty();
    if ((!usesMaf && !usesFasta) || options.algorithm.empty()) {
        throw std::invalid_argument(
                "a FASTA pair or --maf/--maf-block, plus --algorithm, "
                "is required");
    }
    return options;
}

void runSelected(const Options &options, const SequencePair &pair) {
    const anchorwave::AlignmentSelectionPlan plan =
            anchorwave::makeAlignmentSelectionPlan(
                    pair.query, pair.reference, options.windowWidth,
                    -6, -8, -2, -75, -1, 0.0,
                    options.memoryBudgetBytes, options.memoryBudgetBytes);
    if (options.algorithm == "ksw2-score-certified" ||
        options.algorithm == "all") {
        runKsw2ScoreCertified(pair, plan);
    }
    if (options.algorithm == "ksw2-full" || options.algorithm == "all") {
        runKsw2Full(pair, plan);
    }
    if (options.algorithm == "ksw2-banded" || options.algorithm == "all") {
        runApproximate(pair, plan, true);
    }
    if (options.algorithm == "sliding-window" || options.algorithm == "all") {
        runApproximate(pair, plan, false);
    }
    if (options.algorithm == "singletrack" || options.algorithm == "all") {
        runWfa("WFA_SINGLETRACK_HIGH", anchorwave::alignWithSingletrackWfa,
               pair, options, plan,
               findEstimate(plan,
                       anchorwave::AlignmentCandidate::SingletrackWfa));
    }
    if (options.algorithm == "wfa-high" || options.algorithm == "all") {
        runWfa("WFA_HIGH", anchorwave::alignWithStandardWfa, pair, options,
               plan, findEstimate(
                       plan, anchorwave::AlignmentCandidate::StandardWfa));
    }
    if (options.algorithm == "wfa-medium" || options.algorithm == "all") {
        runWfa("WFA_MEDIUM", anchorwave::alignWithMediumWfa, pair, options,
               plan, findEstimate(
                       plan, anchorwave::AlignmentCandidate::MediumWfa));
    }
    if (options.algorithm == "wfa-low" || options.algorithm == "all") {
        runWfa("WFA_LOW", anchorwave::alignWithLowWfa, pair, options, plan,
               findEstimate(plan,
                       anchorwave::AlignmentCandidate::LowWfa));
    }
    if (options.algorithm == "biwfa" || options.algorithm == "all") {
        runWfa("BIWFA", anchorwave::alignWithBiWfa, pair, options, plan,
               nullptr);
    }
    if (options.algorithm == "sliding-block-one-pass") {
        runSlidingBlock(pair, false);
    }
    if (options.algorithm == "sliding-block-two-pass") {
        runSlidingBlock(pair, true);
    }
    if (options.algorithm != "ksw2-score-certified" &&
        options.algorithm != "ksw2-full" &&
        options.algorithm != "ksw2-banded" &&
        options.algorithm != "sliding-window" &&
        options.algorithm != "singletrack" &&
        options.algorithm != "wfa-high" &&
        options.algorithm != "wfa-medium" &&
        options.algorithm != "wfa-low" && options.algorithm != "biwfa" &&
        options.algorithm != "sliding-block-one-pass" &&
        options.algorithm != "sliding-block-two-pass" &&
        options.algorithm != "all") {
        throw std::invalid_argument("unknown algorithm: " + options.algorithm);
    }
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parseOptions(argc, argv);
        SequencePair pair;
        if (!options.mafPath.empty()) {
            if (!readMafBlock(options.mafPath, options.mafBlock, pair)) {
                throw std::runtime_error(
                        "could not read two sequence rows from requested MAF block");
            }
        } else if (!readSingleFasta(options.referenceFasta, pair.reference) ||
                   !readSingleFasta(options.queryFasta, pair.query)) {
            throw std::runtime_error("could not read single-record FASTA pair");
        }
        printHeader();
        runSelected(options, pair);
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        usage(std::cerr);
        return EXIT_FAILURE;
    }
}
