#pragma once

#include <cstddef>
#include <cstdio>
#include <iosfwd>
#include <string>

namespace anchorwave {

// Compare an aligned row with its source without materialising an ungapped
// copy. This keeps validation memory proportional to the current interval.
bool ungappedSequenceEquals(const std::string &alignment,
                            const std::string &source);

// Accumulates only alignment lengths when whole-MAF output is disabled. When
// enabled, aligned rows are spooled to anonymous temporary files and copied to
// the final stream at block completion. The in-memory footprint is bounded by
// the small stdio transfer buffer instead of the chromosome/block length.
class AlignmentBlockBuffer {
public:
    explicit AlignmentBlockBuffer(bool spoolRows);
    ~AlignmentBlockBuffer();

    AlignmentBlockBuffer(const AlignmentBlockBuffer &) = delete;
    AlignmentBlockBuffer &operator=(const AlignmentBlockBuffer &) = delete;

    void append(const std::string &referenceAlignment,
                const std::string &queryAlignment,
                const std::string &referenceSource,
                const std::string &querySource);

    void reset();
    bool empty() const;
    std::size_t referenceLength() const;
    std::size_t queryLength() const;
    std::size_t alignmentColumns() const;

    void writeReference(std::ostream &output);
    void writeQuery(std::ostream &output);

private:
    static std::FILE *openSpool();
    static void appendTo(std::FILE *file, const std::string &text);
    static void copyTo(std::FILE *file, std::ostream &output);
    void closeSpools();

    bool spoolRows_;
    std::FILE *referenceSpool_;
    std::FILE *querySpool_;
    std::size_t referenceLength_;
    std::size_t queryLength_;
    std::size_t alignmentColumns_;
};

}  // namespace anchorwave
