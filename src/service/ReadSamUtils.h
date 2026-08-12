#pragma once

#include "../model/Transcript.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace anchorwave {
namespace read_sam_detail {

inline void retainTopCopyScores(std::vector<double> &scores, double score,
                                int expectedCopies) {
    if (expectedCopies <= 0) {
        return;
    }
    const std::size_t limit =
            static_cast<std::size_t>(expectedCopies) + 1;
    if (scores.size() < limit) {
        scores.push_back(score);
        return;
    }
    const auto smallest = std::min_element(scores.begin(), scores.end());
    if (score > *smallest) {
        *smallest = score;
    }
}

inline bool secondaryCopyExceeds(const std::vector<double> &scores,
                                 int expectedCopies,
                                 double similarityFraction) {
    if (expectedCopies <= 0 ||
        scores.size() <= static_cast<std::size_t>(expectedCopies)) {
        return false;
    }
    const double best = *std::max_element(scores.begin(), scores.end());
    const double extra = *std::min_element(scores.begin(), scores.end());
    return best > 0.0 && extra / best > similarityFraction;
}

struct SamFields {
    std::string queryName;
    int flag = 0;
    std::string referenceName;
    int32_t position = 0;
    size_t cigarBegin = 0;
    size_t cigarEnd = 0;
};

inline bool parseNonNegativeInt(const std::string &text, const size_t begin, const size_t end, int32_t &value) {
    if (begin >= end) {
        return false;
    }

    int64_t parsed = 0;
    for (size_t i = begin; i < end; ++i) {
        const char c = text[i];
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + (c - '0');
        if (parsed > std::numeric_limits<int32_t>::max()) {
            return false;
        }
    }
    value = static_cast<int32_t>(parsed);
    return true;
}

// Parse only the first six mandatory SAM columns. SAM lines often contain many
// optional tags, so splitting the entire line creates avoidable strings.
inline bool parseSamFields(const std::string &line, SamFields &fields) {
    size_t fieldBegin = 0;
    size_t fieldEnd = line.find('\t', fieldBegin);
    if (fieldEnd == std::string::npos) {
        return false;
    }
    fields.queryName.assign(line, fieldBegin, fieldEnd - fieldBegin);

    fieldBegin = fieldEnd + 1;
    fieldEnd = line.find('\t', fieldBegin);
    int32_t flag = 0;
    if (fieldEnd == std::string::npos || !parseNonNegativeInt(line, fieldBegin, fieldEnd, flag)) {
        return false;
    }
    fields.flag = flag;

    fieldBegin = fieldEnd + 1;
    fieldEnd = line.find('\t', fieldBegin);
    if (fieldEnd == std::string::npos) {
        return false;
    }
    fields.referenceName.assign(line, fieldBegin, fieldEnd - fieldBegin);

    fieldBegin = fieldEnd + 1;
    fieldEnd = line.find('\t', fieldBegin);
    if (fieldEnd == std::string::npos || !parseNonNegativeInt(line, fieldBegin, fieldEnd, fields.position)) {
        return false;
    }

    // Skip MAPQ.
    fieldBegin = fieldEnd + 1;
    fieldEnd = line.find('\t', fieldBegin);
    if (fieldEnd == std::string::npos) {
        return false;
    }

    fields.cigarBegin = fieldEnd + 1;
    fields.cigarEnd = line.find('\t', fields.cigarBegin);
    if (fields.cigarEnd == std::string::npos) {
        fields.cigarEnd = line.size();
    }
    return fields.cigarBegin < fields.cigarEnd;
}

struct CigarSummary {
    int32_t headClipping = 0;
    int32_t tailClipping = 0;
    int32_t numberOfAlignedBases = 0;
    int32_t referenceSpan = 0;
};

// Scan the CIGAR directly in the SAM line. This avoids a vector allocation and
// one string allocation plus stoi call for every CIGAR operation.
inline bool summarizeCigar(const std::string &line, const size_t begin, const size_t end, CigarSummary &summary) {
    if (begin >= end || line[begin] == '*') {
        return false;
    }

    summary = CigarSummary();
    size_t cursor = begin;
    while (cursor < end) {
        const size_t lengthBegin = cursor;
        while (cursor < end && line[cursor] >= '0' && line[cursor] <= '9') {
            ++cursor;
        }

        int32_t length = 0;
        if (lengthBegin == cursor || cursor >= end ||
            !parseNonNegativeInt(line, lengthBegin, cursor, length)) {
            return false;
        }

        const char operation = line[cursor++];
        const bool isLastOperation = cursor == end;
        if (isLastOperation && (operation == 'H' || operation == 'S')) {
            summary.tailClipping += length;
            continue;
        }

        switch (operation) {
            case 'H':
            case 'S':
                summary.headClipping += length;
                break;
            case 'M':
            case '=':
            case 'X':
                summary.referenceSpan += length;
                summary.numberOfAlignedBases += length;
                break;
            case 'D':
            case 'N':
                summary.referenceSpan += length;
                break;
            case 'I':
            case 'P':
                break;
            default:
                return false;
        }
    }
    return true;
}

struct CdsSegment {
    int32_t cumulativeEnd = 0;
    int32_t firstGenomicPosition = 0;
    int32_t direction = 1;
};

class CdsCoordinateIndex {
public:
    static CdsCoordinateIndex build(const Transcript &transcript) {
        CdsCoordinateIndex index;
        const std::vector<GenomeBasicFeature> &cds = transcript.getCdsVector();
        int32_t cumulativeEnd = 0;

        if (transcript.getStrand() == POSITIVE) {
            for (size_t i = 0; i < cds.size(); ++i) {
                cumulativeEnd += cds[i].getEnd() - cds[i].getStart() + 1;
                index.segments_.push_back({cumulativeEnd, cds[i].getStart(), 1});
            }
        } else {
            for (size_t i = cds.size(); i > 0; --i) {
                const GenomeBasicFeature &feature = cds[i - 1];
                cumulativeEnd += feature.getEnd() - feature.getStart() + 1;
                index.segments_.push_back({cumulativeEnd, feature.getEnd(), -1});
            }
        }

        index.cdsLength_ = cumulativeEnd;
        return index;
    }

    int32_t cdsLength() const {
        return cdsLength_;
    }

    int32_t genomicPosition(const int32_t cdsPosition) const {
        if (cdsPosition < 1 || cdsPosition > cdsLength_) {
            // This matches the old positionsMap[invalid_position] behavior in
            // release builds while keeping the invalid input visible in debug.
            assert(false && "CDS position is outside the transcript");
            return 0;
        }

        const std::vector<CdsSegment>::const_iterator segment = std::lower_bound(
            segments_.begin(), segments_.end(), cdsPosition,
            [](const CdsSegment &candidate, const int32_t position) {
                return candidate.cumulativeEnd < position;
            });
        assert(segment != segments_.end());
        const int32_t previousEnd = segment == segments_.begin() ? 0 : (segment - 1)->cumulativeEnd;
        return segment->firstGenomicPosition +
               segment->direction * (cdsPosition - previousEnd - 1);
    }

private:
    std::vector<CdsSegment> segments_;
    int32_t cdsLength_ = 0;
};

struct GenomicInterval {
    int32_t start = 0;
    int32_t end = 0;
};

inline GenomicInterval alignedTranscriptInterval(const CdsCoordinateIndex &index,
                                                  const int32_t headClipping,
                                                  const int32_t tailClipping,
                                                  const bool reverseAlignment) {
    int32_t firstPosition = 0;
    int32_t lastPosition = 0;
    if (reverseAlignment) {
        firstPosition = index.genomicPosition(index.cdsLength() - headClipping);
        lastPosition = index.genomicPosition(1 + tailClipping);
    } else {
        firstPosition = index.genomicPosition(1 + headClipping);
        lastPosition = index.genomicPosition(index.cdsLength() - tailClipping);
    }

    GenomicInterval interval;
    interval.start = std::min(firstPosition, lastPosition);
    interval.end = std::max(firstPosition, lastPosition);
    return interval;
}

}  // namespace read_sam_detail
}  // namespace anchorwave
