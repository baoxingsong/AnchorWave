//
// Created by song on 8/4/18.
//

#include "TransferGffWithNucmerResult.h"
#include "ReadSamUtils.h"
#include "NovelAnchorThreshold.h"

#include <atomic>
#include <chrono>
#include <memory>
#include <unordered_map>
#include <utility>

namespace {

using TranscriptIndex =
        std::unordered_map<std::string, const Transcript *>;

struct NovelAnchorRequest {
    std::string refChr;
    std::string queryChr;
    std::string refSequence;
    std::string querySequence;
    uint32_t startRef = 0;
    uint32_t startQuery = 0;
    STRAND strand = POSITIVE;
};

struct NovelAnchorResult {
    uint32_t blacklistStartRef = 0;
    bool shouldBlacklist = false;
    bool found = false;
    AlignmentMatch anchor;
};

// The quota pipeline has slightly different blacklist and negative-strand
// scoring semantics from the non-quota novel-anchor search.  Keep those
// details in a separate result type so gap searches can run concurrently,
// while their effects are still committed in the original scan order.
struct QuotaNovelAnchorResult {
    uint32_t requestStartRef = 0;
    uint32_t discoveredStartRef = 0;
    bool validForwardMapping = false;
    bool found = false;
    AlignmentMatch anchor;
};

uint64_t saturatingAddCost(uint64_t first, uint64_t second) {
    return second > std::numeric_limits<uint64_t>::max() - first
           ? std::numeric_limits<uint64_t>::max() : first + second;
}

uint64_t quotaNovelAnchorBlockEstimatedCost(
        const std::vector<AlignmentMatch> &matches) {
    if (matches.size() <= 1) {
        return 0;
    }
    uint64_t cost = 0;
    for (std::size_t index = 1; index < matches.size(); ++index) {
        const AlignmentMatch &previous = matches[index - 1];
        const AlignmentMatch &current = matches[index];
        if (previous.getStrand() != current.getStrand() ||
            current.getRefStartPos() <= previous.getRefEndPos()) {
            continue;
        }
        const uint64_t referenceLength =
                current.getRefStartPos() - previous.getRefEndPos() - 1;
        uint64_t queryLength = 0;
        if (current.getStrand() == POSITIVE &&
            current.getQueryStartPos() > previous.getQueryEndPos()) {
            queryLength = current.getQueryStartPos() -
                          previous.getQueryEndPos() - 1;
        } else if (current.getStrand() == NEGATIVE &&
                   previous.getQueryStartPos() > current.getQueryEndPos()) {
            queryLength = previous.getQueryStartPos() -
                          current.getQueryEndPos() - 1;
        }
        if (queryLength != 0) {
            cost = saturatingAddCost(
                    cost, anchorwave::anchorGapTaskEstimatedCost(
                                  referenceLength, queryLength));
        }
    }
    return std::max(cost,
                    anchorwave::anchorTaskEstimatedCost(matches.size()));
}

NovelAnchorResult alignNovelAnchorGap(
        const NovelAnchorRequest &request,
        int32_t matchingScore, int32_t mismatchingPenalty,
        int32_t openGapPenalty, int32_t extendGapPenalty,
        int k, int w, int isHpc, int64_t windowWidth,
        double minimumSimilarity) {
    NovelAnchorResult result;
    result.blacklistStartRef = request.startRef;

    mm_idxopt_t indexOptions;
    mm_mapopt_t mappingOptions;
    mm_set_opt(0, &indexOptions, &mappingOptions);
    mappingOptions.flag |= MM_F_CIGAR;
    mappingOptions.flag |= MM_F_NO_PRINT_2ND;
    mappingOptions.bw = windowWidth / 5;
    mappingOptions.flag |= MM_F_NO_LJOIN;
    mappingOptions.a = matchingScore;
    mappingOptions.b = mismatchingPenalty;
    mappingOptions.q = openGapPenalty;
    mappingOptions.e = extendGapPenalty;
    mappingOptions.q2 = openGapPenalty;
    mappingOptions.e2 = extendGapPenalty;
    // Preserve the legacy effective setting: the second assignment replaced
    // windowWidth/5 before mapping began.
    mappingOptions.bw = k;
    mappingOptions.mid_occ = 2;
    mappingOptions.min_cnt = 2;

    // Both minimap2 entry points accept const input and the request owns these
    // strings for the complete call. Avoid two full gap-sized staging copies.
    const char *referenceSequences[1] = {request.refSequence.c_str()};
    const char *sequenceNames[1] = {"temp"};
    char queryName[] = "temp";
    mm_idx_t *index = mm_idx_str(
            w, k, isHpc, 2, 1, referenceSequences, sequenceNames);
    mm_mapopt_update(&mappingOptions, index);
    mm_tbuf_t *threadBuffer = mm_tbuf_init();
    int regionCount = 0;
    mm_reg1_t *regions = mm_map(
            index, static_cast<int>(request.querySequence.size()),
            request.querySequence.c_str(), &regionCount, threadBuffer,
            &mappingOptions, queryName);

    if (regionCount > 0 && regions[0].rev == 0) {
        mm_reg1_t *region = &regions[0];
        const int32_t newAnchorRefEnd = region->re - 1;
        const int32_t newAnchorQueryEnd = region->qe - 1;
        const double length = region->re - region->rs + 1.0;
        double matchingBases = 0.0;
        for (uint32_t cigarIndex = 0;
             cigarIndex < region->p->n_cigar; ++cigarIndex) {
            if ("MIDNSH"[region->p->cigar[cigarIndex] & 0xf] == 'M') {
                matchingBases += region->p->cigar[cigarIndex] >> 4;
            }
        }

        if (matchingBases / length > minimumSimilarity) {
            const uint32_t refStart = request.startRef + region->rs;
            const uint32_t refEnd = request.startRef + newAnchorRefEnd;
            uint32_t queryStart = 0;
            uint32_t queryEnd = 0;
            if (request.strand == POSITIVE) {
                queryStart = request.startQuery + region->qs;
                queryEnd = request.startQuery + newAnchorQueryEnd;
            } else {
                queryStart = request.startQuery - newAnchorQueryEnd;
                queryEnd = request.startQuery - region->qs;
            }
            const std::string alignmentName =
                    "localAlignment_" + request.refChr + "_" +
                    std::to_string(refStart) + "_" + std::to_string(refEnd);
            double score = region->score / 2;
            score /= ((region->re - region->rs) * matchingScore);
            result.anchor = AlignmentMatch(
                    request.refChr, request.queryChr,
                    refStart, refEnd, queryStart, queryEnd,
                    score, request.strand, alignmentName, alignmentName);
            result.found = true;
        }
    } else {
        result.shouldBlacklist = true;
    }

    for (int regionIndex = 0; regionIndex < regionCount; ++regionIndex) {
        free(regions[regionIndex].p);
    }
    free(regions);
    mm_tbuf_destroy(threadBuffer);
    mm_idx_destroy(index);
    return result;
}

QuotaNovelAnchorResult alignQuotaNovelAnchorGap(
        const NovelAnchorRequest &request,
        int32_t matchingScore, int32_t mismatchingPenalty,
        int32_t openGapPenalty, int32_t extendGapPenalty,
        int k, int w, int isHpc, int64_t windowWidth,
        double minimumSimilarity) {
    QuotaNovelAnchorResult result;
    result.requestStartRef = request.startRef;

    mm_idxopt_t indexOptions;
    mm_mapopt_t mappingOptions;
    mm_set_opt(0, &indexOptions, &mappingOptions);
    mappingOptions.flag |= MM_F_CIGAR;
    mappingOptions.flag |= MM_F_NO_PRINT_2ND;
    mappingOptions.bw = windowWidth / 5;
    mappingOptions.flag |= MM_F_NO_LJOIN;
    mappingOptions.a = matchingScore;
    mappingOptions.b = mismatchingPenalty;
    mappingOptions.q = openGapPenalty;
    mappingOptions.e = extendGapPenalty;
    mappingOptions.q2 = openGapPenalty;
    mappingOptions.e2 = extendGapPenalty;
    // Preserve the legacy effective setting: this assignment replaced the
    // earlier windowWidth/5 value before minimap2 was invoked.
    mappingOptions.bw = k;
    mappingOptions.mid_occ = 2;
    mappingOptions.min_cnt = 2;

    const char *referenceSequences[1] = {request.refSequence.c_str()};
    const char *sequenceNames[1] = {"temp"};
    char queryName[] = "temp";
    mm_idx_t *index = mm_idx_str(
            w, k, isHpc, 2, 1, referenceSequences, sequenceNames);
    mm_mapopt_update(&mappingOptions, index);
    mm_tbuf_t *threadBuffer = mm_tbuf_init();
    int regionCount = 0;
    mm_reg1_t *regions = mm_map(
            index, static_cast<int>(request.querySequence.size()),
            request.querySequence.c_str(), &regionCount, threadBuffer,
            &mappingOptions, queryName);

    if (regionCount > 0 && regions[0].rev == 0) {
        result.validForwardMapping = true;
        mm_reg1_t *region = &regions[0];
        result.discoveredStartRef = request.startRef + region->rs;
        const int32_t newAnchorRefEnd = region->re - 1;
        const int32_t newAnchorQueryEnd = region->qe - 1;
        const double length = region->re - region->rs + 1.0;
        double matchingBases = 0.0;
        for (uint32_t cigarIndex = 0;
             cigarIndex < region->p->n_cigar; ++cigarIndex) {
            if ("MIDNSH"[region->p->cigar[cigarIndex] & 0xf] == 'M') {
                matchingBases += region->p->cigar[cigarIndex] >> 4;
            }
        }

        if (matchingBases / length > minimumSimilarity) {
            const uint32_t refStart = result.discoveredStartRef;
            const uint32_t refEnd = request.startRef + newAnchorRefEnd;
            uint32_t queryStart = 0;
            uint32_t queryEnd = 0;
            double scoreDenominator = 0.0;
            if (request.strand == POSITIVE) {
                queryStart = request.startQuery + region->qs;
                queryEnd = request.startQuery + newAnchorQueryEnd;
                scoreDenominator =
                        std::abs(region->re - region->rs) * matchingScore;
            } else {
                queryStart = request.startQuery - newAnchorQueryEnd;
                queryEnd = request.startQuery - region->qs;
                // Preserve the quota implementation's historical +1 on the
                // negative strand; changing it alters anchor scores/output.
                scoreDenominator =
                        std::abs(region->re - region->rs + 1) * matchingScore;
            }
            const std::string alignmentName =
                    "localAlignment_" + request.refChr + "_" +
                    std::to_string(refStart) + "_" +
                    std::to_string(refEnd);
            double score = region->score / 2;
            score /= scoreDenominator;
            result.anchor = AlignmentMatch(
                    request.refChr, request.queryChr,
                    refStart, refEnd, queryStart, queryEnd,
                    score, request.strand, alignmentName, alignmentName);
            result.found = true;
        }
    }

    for (int regionIndex = 0; regionIndex < regionCount; ++regionIndex) {
        free(regions[regionIndex].p);
    }
    free(regions);
    mm_tbuf_destroy(threadBuffer);
    mm_idx_destroy(index);
    return result;
}

}  // namespace

void readSam(std::vector<AlignmentMatch> &alignmentMatchsMapT, std::ifstream &infile,
             const TranscriptIndex &transcriptHashMap, int &expectCopy, const double &minimumSimilarity,
             double &secondarySimilarity, std::set<std::string> &blackGeneList, const std::string &anchorSequenceFile,
             int32_t &matchingScore, int32_t &mismatchingPenalty, int32_t &openGapPenalty1,
             int32_t &extendGapPenalty1, int &k, bool &H, int &w,
             std::map<std::string, std::tuple<std::string, long, long, int> > &queryGenome) {

    // The previous implementation loaded the anchor FASTA and fetched query
    // sequence for every M operation to compute a local score that was never
    // consumed. Keep the parameter for API compatibility, but avoid that I/O.
    (void)anchorSequenceFile;

    using anchorwave::read_sam_detail::CdsCoordinateIndex;
    using anchorwave::read_sam_detail::CigarSummary;
    using anchorwave::read_sam_detail::GenomicInterval;
    using anchorwave::read_sam_detail::SamFields;

    std::string line;
    SamFields fields;
    std::unordered_map<const Transcript *, CdsCoordinateIndex> coordinateIndexes;
    std::unordered_map<std::string, std::pair<std::string, int32_t>> lastAlignment;
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> geneScores;

    while (std::getline(infile, line)) { // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if (line.compare(0, 3, "@PG") == 0) {
            std::vector<std::string> elements;
            std::vector<std::string> elements2;
            char seperator = ' ';
            char seperator2 = ',';
            split(line, seperator, elements);
            if (elements[0].find("ID:minimap2") != std::string::npos) {
                std::cout << "using parameters detected from the input SAM file for novel anchors identification" << std::endl;
                for (size_t i = 0; i < elements.size(); ++i) {
                    std::string element = elements[i];
                    if (element == "-k" and (i + 1) < elements.size()) {
                        k = std::stoi(elements[i + 1]);
                        w = 0.666 * k;
                    }
                    if (element == "-A" and (i + 1) < elements.size()) {
                        matchingScore = std::stoi(elements[i + 1]);
                    }
                    if (element == "-B" and (i + 1) < elements.size()) {
                        mismatchingPenalty = std::stoi(elements[i + 1]);
                    }
                    if (element == "-O" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        openGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-E" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        extendGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-H") {
                        H = true;
                    }
                }
            }
        }

        if (line.size() > 1 && line[0] != '@') { //ignore the header
            if (!anchorwave::read_sam_detail::parseSamFields(line, fields)) {
                continue;
            }

            const auto transcriptIt = transcriptHashMap.find(fields.queryName);
            if (fields.referenceName != "*" && transcriptIt != transcriptHashMap.end() &&
                queryGenome.find(fields.referenceName) != queryGenome.end()) { // ignore non-mapping records
                const Transcript &transcript = *transcriptIt->second;
                const std::string &geneName = fields.queryName;
                const std::string &queryChr = fields.referenceName;
                const std::string &databaseChr = transcript.getChromeSomeName();
                const int32_t queryStart = fields.position;
                const int samFlag = fields.flag;
                // Preserve the legacy orientation test exactly. For the
                // unpaired minimap2 records used here this is normally 0/16.
                const bool reverseAlignment = samFlag % 32 != 0;

                CigarSummary cigar;
                if (!anchorwave::read_sam_detail::summarizeCigar(
                        line, fields.cigarBegin, fields.cigarEnd, cigar)) {
                    continue;
                }
                const int32_t queryEnd = queryStart + cigar.referenceSpan - 1;
                std::unordered_map<const Transcript *, CdsCoordinateIndex>::iterator coordinateIt =
                    coordinateIndexes.find(&transcript);
                if (coordinateIt == coordinateIndexes.end()) {
                    coordinateIt = coordinateIndexes.emplace(
                        &transcript, CdsCoordinateIndex::build(transcript)).first;
                }
                const CdsCoordinateIndex &coordinateIndex = coordinateIt->second;
                const GenomicInterval databaseInterval =
                    anchorwave::read_sam_detail::alignedTranscriptInterval(
                        coordinateIndex, cigar.headClipping, cigar.tailClipping, reverseAlignment);
                const int32_t databaseStart = databaseInterval.start;
                const int32_t databaseEnd = databaseInterval.end;
                assert(databaseStart < databaseEnd);

                std::unordered_map<std::string, std::pair<std::string, int32_t>>::iterator lastIt =
                    lastAlignment.find(geneName);
                if (lastIt != lastAlignment.end() && lastIt->second.first == queryChr &&
                    std::min(std::abs(lastIt->second.second - queryEnd),
                             std::abs(lastIt->second.second - queryStart)) <
                        std::abs(transcript.getPStart() - transcript.getPEnd())) {
                    blackGeneList.insert(geneName);
                } // remove those genes generated weired alignment

                lastAlignment[geneName] = std::make_pair(queryChr, queryEnd);

                const double thisScore = static_cast<double>(cigar.numberOfAlignedBases) /
                                         static_cast<double>(coordinateIndex.cdsLength());
                if (thisScore > minimumSimilarity) {
                    const bool sameOrientation =
                        reverseAlignment == (transcript.getStrand() == NEGATIVE);
                    alignmentMatchsMapT.emplace_back(
                        databaseChr, queryChr, databaseStart, databaseEnd, queryStart, queryEnd,
                        thisScore, sameOrientation ? POSITIVE : NEGATIVE, geneName, geneName);
                    anchorwave::read_sam_detail::retainTopCopyScores(
                            geneScores[geneName][queryChr], thisScore,
                            expectCopy);
                }
            }
        }
    }

    if (expectCopy > 0) {
        for (auto &geneScore : geneScores) {
            for (auto &chromosomeScore : geneScore.second) {
                std::vector<double> &scores = chromosomeScore.second;
                if (anchorwave::read_sam_detail::secondaryCopyExceeds(
                            scores, expectCopy, secondarySimilarity)) {
                    blackGeneList.insert(geneScore.first);
                }
            }
        }
    }

    alignmentMatchsMapT.erase(
        std::remove_if(alignmentMatchsMapT.begin(), alignmentMatchsMapT.end(),
                       [&blackGeneList](const AlignmentMatch &match) {
                           return blackGeneList.find(match.getReferenceGeneName()) != blackGeneList.end();
                       }),
        alignmentMatchsMapT.end());

    if (alignmentMatchsMapT.empty()) {
        std::cout << "there is no match anchor found in the input sam file" << std::endl;
        std::exit(1);
    }
}

void readSam(std::vector<AlignmentMatch> &alignmentMatchsMapT, std::ifstream &infile, const TranscriptIndex &transcriptHashMap, int &expectCopy, const double &minimumSimilarity,
             double &secondarySimilarity, std::set<std::string> &blackGeneList,
             int32_t &matchingScore, int32_t &mismatchingPenalty, int32_t &openGapPenalty1,
             int32_t &extendGapPenalty1, int &k, bool &H, int &w) {

    using anchorwave::read_sam_detail::CdsCoordinateIndex;
    using anchorwave::read_sam_detail::CigarSummary;
    using anchorwave::read_sam_detail::GenomicInterval;
    using anchorwave::read_sam_detail::SamFields;

    std::string line;
    SamFields fields;
    std::unordered_map<const Transcript *, CdsCoordinateIndex> coordinateIndexes;
    std::unordered_map<std::string, std::pair<std::string, int32_t>> lastAlignment;
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> geneScores;

    while (std::getline(infile, line)) { // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if (line.compare(0, 3, "@PG") == 0) {
            std::vector<std::string> elements;
            std::vector<std::string> elements2;
            char seperator = ' ';
            char seperator2 = ',';
            split(line, seperator, elements);

            if (elements[0].find("ID:minimap2") != std::string::npos) {
//                std::cout << "using parameters detected from the input SAM file for novel anchors identification" << std::endl;
                for (size_t i = 0; i < elements.size(); ++i) {
                    std::string element = elements[i];
                    if (element == "-k" and (i + 1) < elements.size()) {
                        k = std::stoi(elements[i + 1]);
                        w = 0.666 * k;
                    }
                    if (element == "-A" and (i + 1) < elements.size()) {
                        matchingScore = std::stoi(elements[i + 1]);
                    }
                    if (element == "-B" and (i + 1) < elements.size()) {
                        mismatchingPenalty = std::stoi(elements[i + 1]);
                    }
                    if (element == "-O" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        openGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-E" and (i + 1) < elements.size()) {
                        split(elements[i + 1], seperator2, elements2);
                        extendGapPenalty1 = std::stoi(elements2[0]);
                    }
                    if (element == "-H") {
                        H = true;
                    }
                }
            }
        }

        if (line.size() > 1 && line[0] != '@') { //ignore the header
            if (!anchorwave::read_sam_detail::parseSamFields(line, fields)) {
                continue;
            }

            const auto transcriptIt = transcriptHashMap.find(fields.queryName);
            if (fields.referenceName != "*" && transcriptIt != transcriptHashMap.end()) { // ignore non-mapping records
                const Transcript &transcript = *transcriptIt->second;
                const std::string &geneName = fields.queryName;
                const std::string &queryChr = fields.referenceName;
                const std::string &databaseChr = transcript.getChromeSomeName();
                const int32_t queryStart = fields.position;
                const int samFlag = fields.flag;
                // Preserve the legacy orientation test exactly. For the
                // unpaired minimap2 records used here this is normally 0/16.
                const bool reverseAlignment = samFlag % 32 != 0;

                CigarSummary cigar;
                if (!anchorwave::read_sam_detail::summarizeCigar(
                        line, fields.cigarBegin, fields.cigarEnd, cigar)) {
                    continue;
                }
                const int32_t queryEnd = queryStart + cigar.referenceSpan - 1;

                std::unordered_map<const Transcript *, CdsCoordinateIndex>::iterator coordinateIt =
                    coordinateIndexes.find(&transcript);
                if (coordinateIt == coordinateIndexes.end()) {
                    coordinateIt = coordinateIndexes.emplace(
                        &transcript, CdsCoordinateIndex::build(transcript)).first;
                }
                const CdsCoordinateIndex &coordinateIndex = coordinateIt->second;
                const GenomicInterval databaseInterval =
                    anchorwave::read_sam_detail::alignedTranscriptInterval(
                        coordinateIndex, cigar.headClipping, cigar.tailClipping, reverseAlignment);
                const int32_t databaseStart = databaseInterval.start;
                const int32_t databaseEnd = databaseInterval.end;
                assert(databaseStart < databaseEnd);

                std::unordered_map<std::string, std::pair<std::string, int32_t>>::iterator lastIt =
                    lastAlignment.find(geneName);
                if (lastIt != lastAlignment.end() && lastIt->second.first == queryChr &&
                    std::min(std::abs(lastIt->second.second - queryEnd),
                             std::abs(lastIt->second.second - queryStart)) <
                        std::abs(transcript.getPStart() - transcript.getPEnd())) {
                    blackGeneList.insert(geneName);
                } // remove those genes generated weired alignment

                lastAlignment[geneName] = std::make_pair(queryChr, queryEnd);

                const double thisScore = static_cast<double>(cigar.numberOfAlignedBases) /
                                         static_cast<double>(coordinateIndex.cdsLength());
                if (thisScore > minimumSimilarity) {
                    const bool sameOrientation =
                        reverseAlignment == (transcript.getStrand() == NEGATIVE);
                    alignmentMatchsMapT.emplace_back(
                        databaseChr, queryChr, databaseStart, databaseEnd, queryStart, queryEnd,
                        thisScore, sameOrientation ? POSITIVE : NEGATIVE, geneName, geneName);
                    anchorwave::read_sam_detail::retainTopCopyScores(
                            geneScores[geneName][queryChr], thisScore,
                            expectCopy);
                }
            }
        }
    }

    if (expectCopy > 0) {
        for (auto &geneScore : geneScores) {
            for (auto &chromosomeScore : geneScore.second) {
                std::vector<double> &scores = chromosomeScore.second;
                if (anchorwave::read_sam_detail::secondaryCopyExceeds(
                            scores, expectCopy, secondarySimilarity)) {
                    blackGeneList.insert(geneScore.first);
                }
            }
        }
    }

    alignmentMatchsMapT.erase(
        std::remove_if(alignmentMatchsMapT.begin(), alignmentMatchsMapT.end(),
                       [&blackGeneList](const AlignmentMatch &match) {
                           return blackGeneList.find(match.getReferenceGeneName()) != blackGeneList.end();
                       }),
        alignmentMatchsMapT.end());

    if (alignmentMatchsMapT.empty()) {
        std::cout << "there is no match anchor found in the input sam file" << std::endl;
        std::exit(1);
    }
}

void setupAnchorsWithSpliceAlignmentResult(const std::string &gffFilePath, const std::string &cdsSequenceFile, const std::string &samFile, std::map<std::string, std::vector<AlignmentMatch>> &map_v_am,
                                           double &inversion_PENALTY, double &MIN_ALIGNMENT_SCORE, bool &considerInversion, const int &minExon, const int64_t &windowWidth, const double &minimumSimilarity, const double &minimumSimilarity2,
                                           std::map<std::string, std::tuple<std::string, long, long, int> > &map_ref,
                                           std::map<std::string, std::tuple<std::string, long, long, int> > &map_qry,
                                           int &expectedCopies, double &maximumSimilarity,
                                           const std::string &referenceSamFilePath, const int64_t &wfaSize3, const bool &searchForNewAnchors, const bool &exonModel,
                                           const int &maxThreads) {
    const uint64_t minimumNovelAnchorArea =
            anchorwave::novelAnchorMinimumArea(wfaSize3);
    std::ifstream infile(samFile);
    if (!infile.good()) {
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit(1);
    }

    // those are default parameter from minimap2
    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = 4;
    int32_t openGapPenalty1 = 4;
    int32_t extendGapPenalty1 = 2;
    int k = 15;
    int w = 0.666 * k;
    bool H = false;

    anchorwave::AnchorTaskExecutor anchorExecutor(maxThreads);
    std::cerr << "AnchorWave anchor scheduler: requested_threads="
              << maxThreads << ", worker_count="
              << anchorExecutor.workerCount() << std::endl;

    //read genome and gff file begin
    std::map<std::string, std::vector<Transcript> > map_v_ts;
//    if (exonModel) {
//        readGffFile(gffFilePath, map_v_ts, "exon", minExon);
//    }
//    else {
//        readGffFile(gffFilePath, map_v_ts, "CDS", minExon);
//    }

    std::string regex = "([\\s\\S]*)Parent=([a-zA-Z0-9.:_-]+)";
    if (exonModel) {
        readGffFile_exon(gffFilePath, map_v_ts, regex, minExon);
    }
    else {
        readGffFile1(gffFilePath, map_v_ts, regex, minExon);
    }


    std::set<std::string> set_rm_chr;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = map_v_ts.begin(); it != map_v_ts.end(); ++it) {
        if (map_ref.find(it->first) == map_ref.end()) {
            set_rm_chr.insert(it->first);
        }
        if (map_qry.find(it->first) == map_qry.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_ref.begin(); it != map_ref.end(); ++it) {
        if (map_v_ts.find(it->first) == map_v_ts.end()) {
            set_rm_chr.insert(it->first);
        }
        if (map_qry.find(it->first) == map_qry.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_qry.begin(); it != map_qry.end(); ++it) {
        if (map_v_ts.find(it->first) == map_v_ts.end()) {
            set_rm_chr.insert(it->first);
        }
        if (map_ref.find(it->first) == map_ref.end()) {
            set_rm_chr.insert(it->first);
        }
    }

    for (std::string chr: set_rm_chr) {
        std::cerr << "There is not enough anchors found on " << chr << std::endl;
        if (map_v_ts.find(chr) != map_v_ts.end()) {
            map_v_ts.erase(chr);
        }
        if (map_ref.find(chr) != map_ref.end()) {
            map_ref.erase(chr);
        }
        if (map_qry.find(chr) != map_qry.end()) {
            map_qry.erase(chr);
        }
    }

    TranscriptIndex transcriptHashMap;
    for (auto &chromosomeTranscripts : map_v_ts) {
        for (Transcript &transcript : chromosomeTranscripts.second) {
            transcriptHashMap.emplace(transcript.getName(), &transcript);
        }
    }

    //read genome and gff file end

    // set gene black list by reading reference gff file begin
    std::set<std::string> blackGeneList;
    if (referenceSamFilePath.size() > 0) {
        std::ifstream infileReferencSam(referenceSamFilePath);
        if (!infileReferencSam.good()) {
            std::cerr << "error in opening sam file " << referenceSamFilePath << std::endl;
            exit(1);
        }

        std::vector<AlignmentMatch> alignmentMatchsMapT0;
        std::cout << "reading reference sam begin" << std::endl;
        readSam(alignmentMatchsMapT0, infileReferencSam, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w);
        std::cout << "reading reference sam done" << std::endl;
        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;

        for (AlignmentMatch &orthologPair2: alignmentMatchsMapT0) {
            if (alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end()) {
                alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
            }

            if (orthologPair2.getRefChr() == orthologPair2.getQueryChr() && orthologPair2.getStrand() == POSITIVE) {
                alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(
                        std::move(orthologPair2));
            }
        }
        std::vector<AlignmentMatch>().swap(alignmentMatchsMapT0);

        bool keepTandemDuplication = false;
        std::vector<std::pair<std::string, std::vector<AlignmentMatch>>>
                referenceResults;
        referenceResults.reserve(alignmentMatchsMapT.size());
        for (const auto &chromosomeMatches : alignmentMatchsMapT) {
            referenceResults.emplace_back(
                    chromosomeMatches.first, std::vector<AlignmentMatch>());
        }
        for (std::size_t resultIndex = 0;
             resultIndex < referenceResults.size(); ++resultIndex) {
            const auto chromosomeIt = alignmentMatchsMapT.find(
                    referenceResults[resultIndex].first);
            const uint64_t estimatedCost = anchorwave::anchorTaskEstimatedCost(
                    chromosomeIt->second.size());
            std::vector<AlignmentMatch> taskMatches =
                    std::move(chromosomeIt->second);
            anchorExecutor.submit(
                    estimatedCost,
                    [&, resultIndex,
                     matches = std::move(taskMatches)]() mutable {
                        myAlignmentMatchSort(
                                matches, inversion_PENALTY,
                                MIN_ALIGNMENT_SCORE, keepTandemDuplication,
                                false);
                        std::vector<AlignmentMatch> sortedAlignmentMatchs;
                        if (matches.size() > 1) {
                            double taskScoreThreshold = MIN_ALIGNMENT_SCORE;
                            longestPath(
                                    matches, sortedAlignmentMatchs,
                                    keepTandemDuplication,
                                    taskScoreThreshold);
                        } else {
                            sortedAlignmentMatchs = std::move(matches);
                        }
                        referenceResults[resultIndex].second =
                                std::move(sortedAlignmentMatchs);
                    });
        }
        anchorExecutor.waitForIdle();
        std::cerr << "AnchorWave anchor scheduler: reference_tasks="
                  << referenceResults.size() << std::endl;
        for (auto &result : referenceResults) {
            map_v_am[result.first] = std::move(result.second);
        }

        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
            for (size_t rangeIndex = 0; rangeIndex < it->second.size(); ++rangeIndex) {
                if (it->second[rangeIndex].getRefStartPos() != it->second[rangeIndex].getQueryStartPos() || it->second[rangeIndex].getRefEndPos() != it->second[rangeIndex].getQueryEndPos()) {
                    blackGeneList.insert(it->second[rangeIndex].getReferenceGeneName());
                }
            }
        }

        infileReferencSam.close();
        map_v_am.clear();
    }

    // set gene black list by reading reference gff file end

    std::vector<AlignmentMatch> alignmentMatchsMapT0;
    std::cout << "reading qry sam begin" << std::endl;
    readSam(alignmentMatchsMapT0, infile, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, cdsSequenceFile, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w, map_qry);
    std::cout << "reading qry sam end" << std::endl;

    std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;
    for (AlignmentMatch &orthologPair2: alignmentMatchsMapT0) {
        if (alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end()) {
            alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
        }

        if (orthologPair2.getRefChr() == orthologPair2.getQueryChr() && map_v_ts.find(orthologPair2.getRefChr()) != map_v_ts.end() &&
                map_ref.find(orthologPair2.getRefChr()) != map_ref.end() && map_qry.find(orthologPair2.getRefChr()) != map_qry.end()) {
            if (!considerInversion && orthologPair2.getStrand() == NEGATIVE) {

            } else {
                alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(
                        std::move(orthologPair2));
            }
        }
    }
    std::vector<AlignmentMatch>().swap(alignmentMatchsMapT0);

    bool keepTandemDuplication = false;
    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_ref.begin(); it != map_ref.end(); ++it) {
        if ( alignmentMatchsMapT.find( it->first ) == alignmentMatchsMapT.end() ){
            std::cerr << "There is not enough anchors found on " << it->first << std::endl;
        }
    }
    std::vector<std::string> queryTaskChromosomes;
    for (const auto &chromosomeMatches : alignmentMatchsMapT) {
        if (chromosomeMatches.second.size() < 3) {
            std::cerr << "There is not enough anchors found on "
                      << chromosomeMatches.first << std::endl;
        } else {
            queryTaskChromosomes.push_back(chromosomeMatches.first);
        }
    }

    std::vector<std::pair<std::string, std::vector<AlignmentMatch>>>
            queryResults;
    queryResults.reserve(queryTaskChromosomes.size());
    for (const std::string &chromosome : queryTaskChromosomes) {
        queryResults.emplace_back(chromosome, std::vector<AlignmentMatch>());
    }

    // minimap2 exposes this verbosity setting as a process-global variable.
    // Set it once before worker threads begin instead of writing it inside
    // concurrent chromosome tasks.
    mm_verbose = 2;
    std::atomic<std::size_t> novelGapTaskCount(0);
    const std::size_t maximumPendingGapsPerChromosome =
            static_cast<std::size_t>(
                    std::max(1, std::min(8, anchorExecutor.workerCount())));
    for (std::size_t taskIndex = 0;
         taskIndex < queryTaskChromosomes.size(); ++taskIndex) {
        const auto chromosomeIt = alignmentMatchsMapT.find(
                queryTaskChromosomes[taskIndex]);
        const uint64_t estimatedCost = anchorwave::anchorTaskEstimatedCost(
                chromosomeIt->second.size());
        std::vector<AlignmentMatch> taskMatches =
                std::move(chromosomeIt->second);
        anchorExecutor.submit(
                estimatedCost,
                [&, taskIndex, matches = std::move(taskMatches)]() mutable {

        std::vector<AlignmentMatch> temp;

        myAlignmentMatchSort(matches, inversion_PENALTY,
                             MIN_ALIGNMENT_SCORE, keepTandemDuplication,
                             considerInversion);
        if (matches.size() > 1) {
            double taskScoreThreshold = MIN_ALIGNMENT_SCORE;
            longestPath(matches, temp, keepTandemDuplication,
                        taskScoreThreshold);
        } else {
            temp = matches;
        }

        int is_hpc = 0; // no, do not use  homopolymer-compressed (HPC) minimizers.
        if (H) {
            is_hpc = 1;
        }

        bool changed = false;
        if (searchForNewAnchors) {
            changed = true;
        }

        std::set<int32_t> blackList;
        while (changed) {
            std::vector<AlignmentMatch> temp2;
            changed = false;

            anchorwave::AnchorTaskGroup gapTaskGroup;
            std::vector<std::shared_ptr<NovelAnchorResult>> gapResults;
            gapResults.reserve(temp.size() + 1);
            auto scheduleNovelAnchor =
                    [&](std::string refSequence, std::string querySequence,
                        uint32_t gapStartRef, uint32_t gapStartQuery,
                        STRAND strand) {
                NovelAnchorRequest request;
                request.refChr = matches[0].getRefChr();
                request.queryChr = request.refChr;
                request.refSequence = std::move(refSequence);
                request.querySequence = std::move(querySequence);
                request.startRef = gapStartRef;
                request.startQuery = gapStartQuery;
                request.strand = strand;

                const uint64_t estimatedCost =
                        anchorwave::anchorGapTaskEstimatedCost(
                                request.refSequence.size(),
                                request.querySequence.size());
                auto output = std::make_shared<NovelAnchorResult>();
                gapResults.push_back(output);
                novelGapTaskCount.fetch_add(1, std::memory_order_relaxed);
                anchorExecutor.submit(
                        gapTaskGroup, estimatedCost,
                        [request = std::move(request), output,
                         matchingScore, mismatchingPenalty,
                         openGapPenalty1, extendGapPenalty1,
                         k, w, is_hpc, windowWidth,
                         minimumSimilarity2]() {
                            *output = alignNovelAnchorGap(
                                    request, matchingScore,
                                    mismatchingPenalty, openGapPenalty1,
                                    extendGapPenalty1, k, w, is_hpc,
                                    windowWidth, minimumSimilarity2);
                        });
                // Bound queued sequence data per chromosome. The parent helps
                // execute global work when the small pending window is full,
                // so nested parallelism does not trade speed for unbounded RSS.
                anchorExecutor.waitUntilGroupSizeAtMost(
                        gapTaskGroup, maximumPendingGapsPerChromosome);
            };

            size_t startRef = 1;
            size_t startQuery = 1;
            size_t endRef;
            size_t endQuery;

            std::string refChr = matches[0].getRefChr();
            std::string queryChr = refChr;

            STRAND lastStrand = POSITIVE;

            bool hasInversion = false;
            int32_t temp_size = temp.size();
            myAlignmentMatchSort(temp, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, false);
            for (int32_t m = 0; m < temp_size; ++m) {
                AlignmentMatch alignmentMatch = temp[m];
                if (alignmentMatch.getStrand() == NEGATIVE) {
                    hasInversion = true;
                }

                if (lastStrand == POSITIVE && alignmentMatch.getStrand() == POSITIVE) {
                    if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() != startQuery) {
                        endQuery = alignmentMatch.getQueryStartPos() - 1;
                    } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryStartPos() == startQuery) {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                    } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryStartPos() == startQuery) {

                    } else {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        endQuery = alignmentMatch.getQueryStartPos() - 1;

                        std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                        std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                        if (anchorwave::novelAnchorAreaExceeds(
                                    static_cast<uint64_t>(refSeq.size()),
                                    static_cast<uint64_t>(querySeq.size()),
                                    minimumNovelAnchorArea) &&
                            refSeq.size() > static_cast<std::size_t>(k) &&
                            querySeq.size() > static_cast<std::size_t>(k) &&
                            blackList.find(startRef) == blackList.end()) {
                            scheduleNovelAnchor(
                                    std::move(refSeq), std::move(querySeq),
                                    static_cast<uint32_t>(startRef),
                                    static_cast<uint32_t>(startQuery),
                                    POSITIVE);
                        }
                    }
                } else if (lastStrand == NEGATIVE && alignmentMatch.getStrand() == NEGATIVE
                           && alignmentMatch.getRefStartPos() > temp[m-1].getRefEndPos()
                           && alignmentMatch.getQueryEndPos() < temp[m-1].getQueryStartPos()   ) {
                    if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() != startQuery) {

                    } else if (alignmentMatch.getRefStartPos() != startRef && alignmentMatch.getQueryEndPos() == startQuery) {

                    } else if (alignmentMatch.getRefStartPos() == startRef && alignmentMatch.getQueryEndPos() == startQuery) {

                    } else {
                        endRef = alignmentMatch.getRefStartPos() - 1;
                        endQuery = alignmentMatch.getQueryEndPos() + 1;

                        std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                        std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery, alignmentMatch.getStrand());

                        if (anchorwave::novelAnchorAreaExceeds(
                                    static_cast<uint64_t>(refSeq.size()),
                                    static_cast<uint64_t>(querySeq.size()),
                                    minimumNovelAnchorArea) &&
                            refSeq.size() > static_cast<std::size_t>(k) &&
                            querySeq.size() > static_cast<std::size_t>(k) &&
                            blackList.find(startRef) == blackList.end()) {
                            scheduleNovelAnchor(
                                    std::move(refSeq), std::move(querySeq),
                                    static_cast<uint32_t>(startRef),
                                    static_cast<uint32_t>(startQuery),
                                    NEGATIVE);
                        }
                    }
                }

                startRef = alignmentMatch.getRefEndPos() + 1;
                startQuery = alignmentMatch.getQueryEndPos() + 1;
                if (alignmentMatch.getStrand() == NEGATIVE) {
                    startQuery = alignmentMatch.getQueryStartPos() - 1;
                }
                lastStrand = alignmentMatch.getStrand();
            }

            if (!hasInversion) {
                endRef = getSequenceSizeFromPath2(map_ref.at(refChr));
                endQuery = getSequenceSizeFromPath2(map_qry.at(queryChr));

                if (startRef > endRef && startQuery <= endQuery) {

                } else if (startRef <= endRef && startQuery > endQuery) {

                } else if (startRef > endRef && startQuery > endQuery) {

                } else { // have bug by aligning b73 with mo17
                    std::string refSeq = getSubsequence2(map_ref, refChr, startRef, endRef);
                    std::string querySeq = getSubsequence2(map_qry, queryChr, startQuery, endQuery);

                    if (anchorwave::novelAnchorAreaExceeds(
                                static_cast<uint64_t>(refSeq.size()),
                                static_cast<uint64_t>(querySeq.size()),
                                minimumNovelAnchorArea) &&
                        refSeq.size() > static_cast<std::size_t>(k) &&
                        querySeq.size() > static_cast<std::size_t>(k) &&
                        blackList.find(startRef) == blackList.end()) {
                        scheduleNovelAnchor(
                                std::move(refSeq), std::move(querySeq),
                                static_cast<uint32_t>(startRef),
                                static_cast<uint32_t>(startQuery), POSITIVE);
                    }
                }
            }

            anchorExecutor.waitForGroup(gapTaskGroup);
            for (const auto &gapResult : gapResults) {
                if (gapResult->shouldBlacklist) {
                    blackList.insert(gapResult->blacklistStartRef);
                }
                if (gapResult->found) {
                    changed = true;
                    temp2.push_back(gapResult->anchor);
                }
            }

            for (AlignmentMatch &alignmentMatch : temp2) {
                temp.push_back(std::move(alignmentMatch));
            }
        }

        if (temp.size() > 0) {
            myAlignmentMatchSort(temp, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);

            std::vector<AlignmentMatch> sortedAlignmentMatchs;

            if (temp.size() > 1) {
                double taskScoreThreshold = MIN_ALIGNMENT_SCORE;
                longestPath(temp, sortedAlignmentMatchs,
                            keepTandemDuplication, taskScoreThreshold);
            } else {
                sortedAlignmentMatchs = std::move(temp);
            }
            queryResults[taskIndex].second = std::move(sortedAlignmentMatchs);
        }
                });
    }

    anchorExecutor.waitForIdle();
    std::cerr << "AnchorWave anchor scheduler: query_tasks="
              << queryResults.size() << ", peak_active_tasks="
              << anchorExecutor.peakActiveTasks()
              << ", novel_gap_tasks=" << novelGapTaskCount.load()
              << ", max_pending_gaps_per_chromosome="
              << maximumPendingGapsPerChromosome
              << std::endl;
    for (auto &result : queryResults) {
        if (!result.second.empty()) {
            map_v_am[result.first] = std::move(result.second);
        }
    }

    infile.close();
}

void setupAnchorsWithSpliceAlignmentResultQuota(const std::string &gffFilePath, const std::string &samFile, const std::string &cdsSequenceFile, std::vector<std::vector<AlignmentMatch>> &alignmentMatchsMap,
                                                double &INDEL_SCORE, double &GAP_OPEN_PENALTY, double &MIN_ALIGNMENT_SCORE, int &MAX_DIST_BETWEEN_MATCHES, int &refMaximumTimes, int &queryMaximumTimes,
                                                double &calculateIndelDistance, const int &minExon, const int64_t &windowWidth, const double &minimumSimilarity,
                                                const double &minimumSimilarity2,
                                                std::map<std::string, std::tuple<std::string, long, long, int> > &map_ref,
                                                std::map<std::string, std::tuple<std::string, long, long, int> > &map_qry,
                                                int &expectedCopies, const int64_t &wfaSize3,
                                                double &maximumSimilarity, const std::string &referenceSamFilePath, bool &searchForNewAnchors, const bool &exonModel,
                                                const int &maxThreads) {
    const uint64_t minimumNovelAnchorArea =
            anchorwave::novelAnchorMinimumArea(wfaSize3);

    // they are default parameter from minimap2
    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = 4;
    int32_t openGapPenalty1 = 4;
    int32_t extendGapPenalty1 = 2;
    int k = 15;
    int w = 0.666 * k;
    bool H = false;

    std::ifstream infile(samFile);
    if (!infile.good()) {
        std::cerr << "error in opening sam file " << samFile << std::endl;
        exit(1);
    }

    // read reference genome and gff begin
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
//    if (exonModel) {
//        readGffFile(gffFilePath, transcriptHashSet, "exon", minExon);
//    } else {
//        readGffFile(gffFilePath, transcriptHashSet, "CDS", minExon);
//    }

    std::string regex = "([\\s\\S]*)Parent=([a-zA-Z0-9.:_-]+)";
    if (exonModel) {
        readGffFile_exon(gffFilePath, transcriptHashSet, regex, minExon);
    }
    else {
        readGffFile1(gffFilePath, transcriptHashSet, regex, minExon);
    }

    // clean reference genome and annotation chromosomes
    std::set<std::string> toRemoveChrs;
    for (std::map<std::string, std::vector<Transcript> >::iterator it = transcriptHashSet.begin();
         it != transcriptHashSet.end(); ++it) {
        if (map_ref.find(it->first) == map_ref.end()) {
            toRemoveChrs.insert(it->first);
        }
    }

    for (std::map<std::string, std::tuple<std::string, long, long, int> >::iterator it = map_ref.begin(); it != map_ref.end(); ++it) {
        if (transcriptHashSet.find(it->first) == transcriptHashSet.end()) {
            toRemoveChrs.insert(it->first);
        }
    }

    for (std::string chr: toRemoveChrs) {
        if (transcriptHashSet.find(chr) != transcriptHashSet.end()) {
            transcriptHashSet.erase(chr);
        }

        if (map_ref.find(chr) != map_ref.end()) {
            map_ref.erase(chr);
        }
    }

    TranscriptIndex transcriptHashMap;
    for (auto &chromosomeTranscripts : transcriptHashSet) {
        for (Transcript &transcript : chromosomeTranscripts.second) {
            transcriptHashMap.emplace(transcript.getName(), &transcript);
        }
    }

    std::set<std::string> blackGeneList;
    if (referenceSamFilePath.size() > 0) {
        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap000;
        std::ifstream infileReferencSam(referenceSamFilePath);
        if (!infileReferencSam.good()) {
            std::cerr << "error in opening sam file " << referenceSamFilePath << std::endl;
            exit(1);
        }

        std::vector<AlignmentMatch> alignmentMatchsMapT0;
        readSam(alignmentMatchsMapT0, infileReferencSam, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w);
        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;

        for (AlignmentMatch &orthologPair2: alignmentMatchsMapT0) {
            if (alignmentMatchsMapT.find(orthologPair2.getRefChr()) == alignmentMatchsMapT.end()) {
                alignmentMatchsMapT[orthologPair2.getRefChr()] = std::vector<AlignmentMatch>();
            }
            if (orthologPair2.getRefChr() == orthologPair2.getQueryChr() && orthologPair2.getStrand() == POSITIVE) {
                alignmentMatchsMapT[orthologPair2.getRefChr()].push_back(
                        std::move(orthologPair2));
            }
        }
        std::vector<AlignmentMatch>().swap(alignmentMatchsMapT0);

        bool keepTandemDuplication = false;
        double inversion_PENALTY = -1;
        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMapT.begin(); it != alignmentMatchsMapT.end(); ++it) {
            myAlignmentMatchSort(it->second, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, false); // since here the considerInversion is false, so the parameters of inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, would not be used
            std::vector<AlignmentMatch> sortedAlignmentMatchs;
            if (it->second.size() > 1) {
                longestPath(it->second, sortedAlignmentMatchs, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
            } else {
                sortedAlignmentMatchs = std::move(it->second);
            }

            alignmentMatchsMap000[it->first] =
                    std::move(sortedAlignmentMatchs);
        }

        for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMap000.begin(); it != alignmentMatchsMap000.end(); ++it) {
            for (size_t rangeIndex = 0; rangeIndex < it->second.size(); ++rangeIndex) {
                if (it->second[rangeIndex].getRefStartPos() != it->second[rangeIndex].getQueryStartPos() ||
                    it->second[rangeIndex].getRefEndPos() != it->second[rangeIndex].getQueryEndPos()) {
                    blackGeneList.insert(it->second[rangeIndex].getReferenceGeneName());
                }
            }
        }

        infileReferencSam.close();
        alignmentMatchsMap.clear();
        alignmentMatchsMap000.clear();
    }

    {
        std::vector<AlignmentMatch> alignmentMatchsMapT;
        readSam(alignmentMatchsMapT, infile, transcriptHashMap, expectedCopies, minimumSimilarity, maximumSimilarity, blackGeneList, cdsSequenceFile, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, k, H, w, map_qry);

        // begin setting index, they are necessary in the longest path approach
        std::map<std::string, int64_t> queryIndex;
        std::map<std::string, std::map<int64_t, AlignmentMatch>> queryIndexMap; // chr, index, AlignmentMatch
        myOrthologPairsSortQueryQuotaV2(alignmentMatchsMapT);
        for (size_t ii = 0; ii < alignmentMatchsMapT.size(); ++ii) {
            if (queryIndex.find(alignmentMatchsMapT[ii].getQueryChr()) == queryIndex.end()) {
                queryIndex[alignmentMatchsMapT[ii].getQueryChr()] = 0;
                queryIndexMap[alignmentMatchsMapT[ii].getQueryChr()] = std::map<int64_t, AlignmentMatch>();
            } else {
                queryIndex[alignmentMatchsMapT[ii].getQueryChr()] = queryIndex[alignmentMatchsMapT[ii].getQueryChr()] + 1;
            }

            alignmentMatchsMapT[ii].setQueryId(queryIndex[alignmentMatchsMapT[ii].getQueryChr()]);
            queryIndexMap[alignmentMatchsMapT[ii].getQueryChr()][queryIndex[alignmentMatchsMapT[ii].getQueryChr()]] = alignmentMatchsMapT[ii];
        }

//        std::cout << "query id setting done" << std::endl;

        myOrthologPairsSortQuotaV2(alignmentMatchsMapT);
        std::map<std::string, int64_t> refIndex; // key is chr
        std::map<std::string, std::map<int64_t, AlignmentMatch>> refIndexMap; // chr, index, AlignmentMatch
        for (size_t ii = 0; ii < alignmentMatchsMapT.size(); ++ii) {
            if (refIndex.find(alignmentMatchsMapT[ii].getRefChr()) == refIndex.end()) {
                refIndex[alignmentMatchsMapT[ii].getRefChr()] = 0;
                refIndexMap[alignmentMatchsMapT[ii].getRefChr()] = std::map<int64_t, AlignmentMatch>();
            }
            else {
                refIndex[alignmentMatchsMapT[ii].getRefChr()] = refIndex[alignmentMatchsMapT[ii].getRefChr()] + 1;
            }
            alignmentMatchsMapT[ii].setRefId(refIndex[alignmentMatchsMapT[ii].getRefChr()]);
            refIndexMap[alignmentMatchsMapT[ii].getRefChr()][refIndex[alignmentMatchsMapT[ii].getRefChr()]] = alignmentMatchsMapT[ii];
        }

        // index setting end
//        std::cout << "reference id setting done" << std::endl;
        orthologPairSortRefForAccelerate(alignmentMatchsMapT);
        longestPathQuotav2(alignmentMatchsMapT, alignmentMatchsMap, refIndexMap, queryIndexMap, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES, refMaximumTimes, queryMaximumTimes, calculateIndelDistance, false, maxThreads);
        for (size_t i = 0; i < alignmentMatchsMap.size(); ++i) {
            std::vector<size_t> keepIndexs;
            std::set<size_t> keepIndexset;
            for (size_t j = 0; j < alignmentMatchsMap[i].size(); ++j) {
                assert(alignmentMatchsMap[i][j].getStrand() == alignmentMatchsMap[i][0].getStrand());
                if (alignmentMatchsMap[i][j].getStrand() == POSITIVE) {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryEndPos() <
                               alignmentMatchsMap[i][j].getQueryStartPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() <
                               alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        std::cout << "line 1329" << std::endl;
                    }
                } else {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryStartPos() >
                               alignmentMatchsMap[i][j].getQueryEndPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() <
                               alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        std::cout << "line 1342" << std::endl;
                    }
                }
            }

            int this_size = alignmentMatchsMap[i].size();
            for (int jj = this_size - 1; jj >= 0; jj--) {
                if (keepIndexset.find(jj) == keepIndexset.end()) {
                    alignmentMatchsMap[i].erase(alignmentMatchsMap[i].begin() + jj);
                }
            }
            bool keepTandemDuplication = false;
            double inversion_PENALTY = -1;
            bool considerInversion = false;
            myAlignmentMatchSort(alignmentMatchsMap[i], inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion); // since here the considerInversion is false, so the parameters of inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, would not be used
        }
    }

    int is_hpc = 0; // no, do not use  homopolymer-compressed (HPC) minimizers.
    if (H) {
        is_hpc = 1;
    }

    mm_verbose = 2; // minimap2 verbosity is process-global; set before workers.
    if (searchForNewAnchors) {
        const auto novelStageStart = std::chrono::steady_clock::now();
        anchorwave::AnchorTaskExecutor novelAnchorExecutor(maxThreads);
        std::size_t novelAnchorBlockTasks = 0;
        std::atomic<std::size_t> novelAnchorGapTasks{0};
        std::vector<double> novelBlockSeconds(alignmentMatchsMap.size(), 0.0);
        for (size_t i = 0; i < alignmentMatchsMap.size(); ++i) {
            if (alignmentMatchsMap[i].size() <= 1) {
                continue;
            }
            ++novelAnchorBlockTasks;
            const uint64_t estimatedCost =
                    quotaNovelAnchorBlockEstimatedCost(alignmentMatchsMap[i]);
            novelAnchorExecutor.submit(estimatedCost, [&, i]() {
            const auto blockStart = std::chrono::steady_clock::now();
            std::set<int32_t> blackList;
            bool changed = true;
            while (alignmentMatchsMap[i].size() > 1 && changed) {
                changed = false;
                const STRAND strand = alignmentMatchsMap[i][0].getStrand();
                const std::string refChr = alignmentMatchsMap[i][0].getRefChr();
                const std::string queryChr = alignmentMatchsMap[i][0].getQueryChr();
                const int32_t snapshotSize =
                        static_cast<int32_t>(alignmentMatchsMap[i].size());
                size_t startRef = alignmentMatchsMap[i][0].getRefStartPos();
                size_t queryCursor = strand == POSITIVE
                                     ? alignmentMatchsMap[i][0].getQueryStartPos()
                                     : alignmentMatchsMap[i][0].getQueryEndPos();

                anchorwave::AnchorTaskGroup gapGroup;
                std::vector<std::shared_ptr<QuotaNovelAnchorResult>> results;
                results.reserve(static_cast<std::size_t>(snapshotSize));
                const std::size_t maximumPendingGaps = std::max<std::size_t>(
                        1, std::min<std::size_t>(8,
                                novelAnchorExecutor.workerCount()));

                for (int32_t matchIndex = 0;
                     matchIndex < snapshotSize; ++matchIndex) {
                    const AlignmentMatch orthologPair =
                            alignmentMatchsMap[i][matchIndex];
                    size_t endRef = 0;
                    size_t queryStart = 0;
                    size_t queryEnd = 0;
                    bool hasGap = false;

                    if (strand == POSITIVE) {
                        if (orthologPair.getRefStartPos() == startRef &&
                            orthologPair.getQueryStartPos() != queryCursor) {
                            queryCursor = orthologPair.getQueryStartPos() - 1;
                        } else if (orthologPair.getRefStartPos() != startRef &&
                                   orthologPair.getQueryStartPos() == queryCursor) {
                            endRef = orthologPair.getRefStartPos() - 1;
                        } else if (orthologPair.getRefStartPos() != startRef ||
                                   orthologPair.getQueryStartPos() != queryCursor) {
                            endRef = orthologPair.getRefStartPos() - 1;
                            queryStart = queryCursor;
                            queryEnd = orthologPair.getQueryStartPos() - 1;
                            hasGap = true;
                        }
                    } else {
                        if (orthologPair.getRefStartPos() == startRef &&
                            orthologPair.getQueryEndPos() != queryCursor) {
                            queryStart = orthologPair.getQueryEndPos() + 1;
                        } else if (orthologPair.getRefStartPos() != startRef &&
                                   orthologPair.getQueryEndPos() == queryCursor) {
                            endRef = orthologPair.getRefStartPos() - 1;
                        } else if (orthologPair.getRefStartPos() != startRef ||
                                   orthologPair.getQueryEndPos() != queryCursor) {
                            endRef = orthologPair.getRefStartPos() - 1;
                            queryStart = orthologPair.getQueryEndPos() + 1;
                            queryEnd = queryCursor;
                            hasGap = true;
                        }
                    }

                    if (hasGap) {
                        std::string refSequence = getSubsequence2(
                                map_ref, refChr, startRef, endRef);
                        std::string querySequence = strand == POSITIVE
                                ? getSubsequence2(map_qry, queryChr,
                                                  queryStart, queryEnd)
                                : getSubsequence2(map_qry, queryChr,
                                                  queryStart, queryEnd, strand);
                        auto result =
                                std::make_shared<QuotaNovelAnchorResult>();
                        result->requestStartRef =
                                static_cast<uint32_t>(startRef);
                        if (anchorwave::novelAnchorAreaExceeds(
                                    static_cast<uint64_t>(refSequence.size()),
                                    static_cast<uint64_t>(querySequence.size()),
                                    minimumNovelAnchorArea) &&
                            refSequence.size() > static_cast<std::size_t>(k) &&
                            querySequence.size() > static_cast<std::size_t>(k) &&
                            blackList.find(startRef) == blackList.end()) {
                            NovelAnchorRequest request;
                            request.refChr = refChr;
                            request.queryChr = queryChr;
                            request.refSequence = std::move(refSequence);
                            request.querySequence = std::move(querySequence);
                            request.startRef = static_cast<uint32_t>(startRef);
                            request.startQuery = static_cast<uint32_t>(
                                    strand == POSITIVE ? queryStart : queryEnd);
                            request.strand = strand;
                            const uint64_t gapCost =
                                    anchorwave::anchorGapTaskEstimatedCost(
                                            request.refSequence.size(),
                                            request.querySequence.size());
                            novelAnchorExecutor.submit(
                                    gapGroup, gapCost,
                                    [&, request = std::move(request), result]() {
                                *result = alignQuotaNovelAnchorGap(
                                        request, matchingScore,
                                        mismatchingPenalty, openGapPenalty1,
                                        extendGapPenalty1, k, w, is_hpc,
                                        windowWidth, minimumSimilarity2);
                                novelAnchorGapTasks.fetch_add(
                                        1, std::memory_order_relaxed);
                            });
                            novelAnchorExecutor.waitUntilGroupSizeAtMost(
                                    gapGroup, maximumPendingGaps - 1);
                        }
                        // Even gaps that do not launch minimap2 participate in
                        // ordered commit: the serial code blacklisted them at
                        // this exact point in the scan.  Applying that effect
                        // early can suppress an earlier concurrent result.
                        results.push_back(std::move(result));
                    }

                    startRef = orthologPair.getRefEndPos() + 1;
                    queryCursor = strand == POSITIVE
                                  ? orthologPair.getQueryEndPos() + 1
                                  : orthologPair.getQueryStartPos() - 1;
                }

                novelAnchorExecutor.waitForGroup(gapGroup);
                for (const auto &result : results) {
                    if (blackList.find(result->requestStartRef) !=
                        blackList.end()) {
                        continue;
                    }
                    if (!result->validForwardMapping) {
                        blackList.insert(result->requestStartRef);
                    } else if (blackList.find(result->discoveredStartRef) ==
                               blackList.end()) {
                        if (result->found) {
                            alignmentMatchsMap[i].push_back(result->anchor);
                            changed = true;
                        }
                        blackList.insert(result->discoveredStartRef);
                    }
                }
                bool keepTandemDuplication = false;
                double inversion_PENALTY = -1;
                bool considerInversion = false;
                myAlignmentMatchSort(alignmentMatchsMap[i], inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);
            }
            novelBlockSeconds[i] = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - blockStart).count();
            });
        }
        novelAnchorExecutor.waitForIdle();
        const double novelWallSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - novelStageStart).count();
        double novelTaskSeconds = 0.0;
        double criticalBlockSeconds = 0.0;
        for (double seconds : novelBlockSeconds) {
            novelTaskSeconds += seconds;
            criticalBlockSeconds = std::max(criticalBlockSeconds, seconds);
        }
        std::cerr << "AnchorWave quota novel-anchor scheduler: worker_count="
                  << novelAnchorExecutor.workerCount()
                  << ", peak_active_tasks="
                  << novelAnchorExecutor.peakActiveTasks()
                  << ", block_tasks=" << novelAnchorBlockTasks
                  << ", gap_tasks="
                  << novelAnchorGapTasks.load(std::memory_order_relaxed)
                  << ", wall_seconds=" << novelWallSeconds
                  << ", task_seconds=" << novelTaskSeconds
                  << ", critical_block_seconds=" << criticalBlockSeconds
                  << std::endl;
    }

    {
        for (size_t i = 0; i < alignmentMatchsMap.size(); ++i) {
            std::vector<size_t> keepIndexs;
            std::set<size_t> keepIndexset;
            for (size_t j = 0; j < alignmentMatchsMap[i].size(); ++j) {
                if (alignmentMatchsMap[i][j].getStrand() == POSITIVE) {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryEndPos() < alignmentMatchsMap[i][j].getQueryStartPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() < alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        // here should never run.
                        std::cout << "removed" << alignmentMatchsMap[i][j].getRefChr() << "\t"
                                  << alignmentMatchsMap[i][j].getRefStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getRefEndPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryChr() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryEndPos() << "\t"
                                  << "+" << "\t"
                                  << alignmentMatchsMap[i][j].getReferenceGeneName() << "\t"
                                  << alignmentMatchsMap[i][j].getScore() << std::endl;
                    }
                } else {
                    if (keepIndexs.size() == 0) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else if (alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getQueryStartPos() > alignmentMatchsMap[i][j].getQueryEndPos() &&
                               alignmentMatchsMap[i][keepIndexs[keepIndexs.size() - 1]].getRefEndPos() < alignmentMatchsMap[i][j].getRefStartPos()) {
                        keepIndexs.push_back(j);
                        keepIndexset.insert(j);
                    } else {
                        std::cout << "removed" << alignmentMatchsMap[i][j].getRefChr() << "\t"
                                  << alignmentMatchsMap[i][j].getRefStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getRefEndPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryChr() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryStartPos() << "\t"
                                  << alignmentMatchsMap[i][j].getQueryEndPos() << "\t"
                                  << "-" << "\t"
                                  << alignmentMatchsMap[i][j].getReferenceGeneName() << "\t"
                                  << alignmentMatchsMap[i][j].getScore() << std::endl;
                    }
                }
            }

            int this_size = alignmentMatchsMap[i].size();
            for (int jj = this_size - 1; jj >= 0; jj--) {
                if (keepIndexset.find(jj) == keepIndexset.end()) {
                    alignmentMatchsMap[i].erase(alignmentMatchsMap[i].begin() + jj);
                }
            }

            bool keepTandemDuplication = false;
            double inversion_PENALTY = -1;
            bool considerInversion = false;

            // since here the considerInversion is false, so the parameters of inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, would not be used
            myAlignmentMatchSort(alignmentMatchsMap[i], inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);
        }
    }
}
