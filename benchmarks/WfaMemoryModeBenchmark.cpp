#include <cstdio>

extern "C" {
#include "WFA2-lib/wavefront/wfa.h"
}

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

struct BenchmarkResult {
    double seconds = 0.0;
    uint64_t peakMemoryBytes = 0;
    int completed = 0;
    int64_t scoreChecksum = 0;
};

wavefront_aligner_t *makeAligner(wavefront_memory_t memoryMode,
                                 uint64_t memoryBudgetBytes) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = 6;
    attributes.affine2p_penalties.gap_opening1 = 8;
    attributes.affine2p_penalties.gap_extension1 = 2;
    attributes.affine2p_penalties.gap_opening2 = 75;
    attributes.affine2p_penalties.gap_extension2 = 1;
    attributes.alignment_form.span = alignment_end2end;
    attributes.alignment_scope = compute_alignment;
    attributes.memory_mode = memoryMode;
    attributes.heuristic.strategy = wf_heuristic_none;
    attributes.system.max_memory_abort = memoryBudgetBytes;
    attributes.system.verbose = 0;
    attributes.system.max_num_threads = 1;
    return wavefront_aligner_new(&attributes);
}

std::string randomDna(std::size_t length, std::mt19937 &generator) {
    static const char bases[] = "ACGT";
    std::uniform_int_distribution<int> base(0, 3);
    std::string sequence;
    sequence.reserve(length);
    for (std::size_t i = 0; i < length; ++i) {
        sequence.push_back(bases[base(generator)]);
    }
    return sequence;
}

std::string mutate(const std::string &reference,
                   double errorRate,
                   std::mt19937 &generator) {
    static const char bases[] = "ACGT";
    std::uniform_real_distribution<double> event(0.0, 1.0);
    std::uniform_int_distribution<int> base(0, 3);
    std::string query;
    query.reserve(reference.size() +
                  static_cast<std::size_t>(reference.size() * errorRate));
    for (char nucleotide : reference) {
        const double draw = event(generator);
        if (draw < errorRate / 3.0) {
            continue;
        }
        if (draw < 2.0 * errorRate / 3.0) {
            char replacement = nucleotide;
            while (replacement == nucleotide) {
                replacement = bases[base(generator)];
            }
            query.push_back(replacement);
            continue;
        }
        if (draw < errorRate) {
            query.push_back(bases[base(generator)]);
        }
        query.push_back(nucleotide);
    }
    return query;
}

BenchmarkResult run(
        const std::vector<std::pair<std::string, std::string>> &pairs,
        int repetitions,
        wavefront_memory_t memoryMode,
        uint64_t memoryBudgetBytes) {
    BenchmarkResult result;
    wavefront_aligner_t *const aligner =
            makeAligner(memoryMode, memoryBudgetBytes);
    if (aligner == nullptr) {
        return result;
    }

    const auto begin = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition) {
        for (const auto &pair : pairs) {
            const std::string &query = pair.first;
            const std::string &reference = pair.second;
            const int status = wavefront_align(
                    aligner,
                    reference.data(), static_cast<int>(reference.size()),
                    query.data(), static_cast<int>(query.size()));
            result.peakMemoryBytes = std::max<uint64_t>(
                    result.peakMemoryBytes,
                    aligner->align_status.memory_used);
            if (status == WF_STATUS_ALG_COMPLETED) {
                ++result.completed;
                result.scoreChecksum += aligner->cigar->score;
            }
        }
    }
    const auto end = std::chrono::steady_clock::now();
    result.seconds =
            std::chrono::duration<double>(end - begin).count();
    wavefront_aligner_delete(aligner);
    return result;
}

void printResult(const char *mode,
                 std::size_t length,
                 double errorRate,
                 std::size_t pairs,
                 int repetitions,
                 const BenchmarkResult &result) {
    const int total = static_cast<int>(pairs) * repetitions;
    const double rate = result.seconds == 0.0
                        ? 0.0
                        : total / result.seconds;
    std::cout << mode << ',' << length << ',' << errorRate << ','
              << pairs << ',' << repetitions << ',' << result.completed << ','
              << result.scoreChecksum << ',' << result.peakMemoryBytes << ','
              << std::fixed << std::setprecision(6) << result.seconds << ','
              << std::setprecision(2) << rate << '\n';
}

}  // namespace

int main(int argc, char **argv) {
    const std::size_t length = argc > 1
                               ? static_cast<std::size_t>(std::stoull(argv[1]))
                               : 1000;
    const double errorRate = argc > 2 ? std::stod(argv[2]) : 0.05;
    const int pairCount = argc > 3 ? std::stoi(argv[3]) : 20;
    const int repetitions = argc > 4 ? std::stoi(argv[4]) : 5;
    const uint64_t memoryBudgetBytes = argc > 5
                                       ? std::stoull(argv[5])
                                       : 10000000000ULL;
    if (length == 0 || errorRate < 0.0 || errorRate >= 1.0 ||
        pairCount <= 0 || repetitions <= 0 || memoryBudgetBytes == 0 ||
        length > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        std::cerr << "usage: anchorwave_wfa_mode_benchmark "
                  << "[length] [error-rate] [pairs] [repetitions] "
                  << "[memory-budget-bytes]\n";
        return EXIT_FAILURE;
    }

    std::mt19937 generator(20260808);
    std::vector<std::pair<std::string, std::string>> pairs;
    pairs.reserve(static_cast<std::size_t>(pairCount));
    for (int i = 0; i < pairCount; ++i) {
        std::string reference = randomDna(length, generator);
        std::string query = mutate(reference, errorRate, generator);
        pairs.emplace_back(std::move(query), std::move(reference));
    }

    // Warm up every exact WFA memory mode before measuring it.
    run(pairs, 1, wavefront_memory_high, memoryBudgetBytes);
    run(pairs, 1, wavefront_memory_med, memoryBudgetBytes);
    run(pairs, 1, wavefront_memory_low, memoryBudgetBytes);
    run(pairs, 1, wavefront_memory_ultralow, memoryBudgetBytes / 3);

    const BenchmarkResult standard = run(
            pairs, repetitions, wavefront_memory_high, memoryBudgetBytes);
    const BenchmarkResult medium = run(
            pairs, repetitions, wavefront_memory_med, memoryBudgetBytes);
    const BenchmarkResult low = run(
            pairs, repetitions, wavefront_memory_low, memoryBudgetBytes);
    const BenchmarkResult bidirectional = run(
            pairs, repetitions, wavefront_memory_ultralow,
            memoryBudgetBytes / 3);

    std::cout << "mode,length,error_rate,pairs,repetitions,completed,"
              << "score_checksum,peak_memory_bytes,seconds,alignments_per_second\n";
    printResult("WFA", length, errorRate, pairs.size(), repetitions,
                standard);
    printResult("WFA_MEDIUM", length, errorRate, pairs.size(), repetitions,
                medium);
    printResult("WFA_LOW", length, errorRate, pairs.size(), repetitions,
                low);
    printResult("BIWFA", length, errorRate, pairs.size(), repetitions,
                bidirectional);

    if (standard.completed != medium.completed ||
        standard.completed != low.completed ||
        standard.completed != bidirectional.completed ||
        standard.scoreChecksum != medium.scoreChecksum ||
        standard.scoreChecksum != low.scoreChecksum ||
        standard.scoreChecksum != bidirectional.scoreChecksum) {
        std::cerr << "exact WFA memory modes differ\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
