#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace anchorwave {
namespace trio {

class MafFormatError : public std::runtime_error {
public:
    explicit MafFormatError(const std::string &message) : std::runtime_error(message) {}
};

// Coordinates follow the MAF specification: start0 is measured on the strand
// shown in text. forwardStart0() converts the covered interval to the forward
// strand and columnForwardPosition0() maps a displayed base to its source base.
struct MafRow {
    std::string taxon;
    std::string source;
    int64_t start0 = 0;
    int64_t size = 0;
    char strand = '+';
    int64_t sourceSize = 0;
    std::string text;
    std::size_t lineNumber = 0;

    int64_t forwardStart0() const;
    int64_t forwardEnd0() const;
    int64_t columnForwardPosition0(std::size_t column) const;
};

struct MafBlock {
    std::map<std::string, std::string> attributes;
    std::vector<MafRow> rows;
    std::size_t lineNumber = 0;
    std::size_t blockIndex = 0;
};

struct MafReadOptions {
    bool requireExactlyTwoRows = true;
    bool requireDnaAlphabet = true;
    bool allowEmptyBlocks = false;
};

// Reads AnchorWave pairwise MAF output without assigning biological reference
// semantics to either row. Row 0 receives leftTaxon and row 1 rightTaxon; the
// orientation is retained only as input provenance.
class MafEvidenceReader {
public:
    static std::vector<MafBlock> read(
        std::istream &input,
        const std::string &sourceName,
        const std::string &leftTaxon,
        const std::string &rightTaxon,
        const MafReadOptions &options = MafReadOptions());

    static std::vector<MafBlock> readFile(
        const std::string &path,
        const std::string &leftTaxon,
        const std::string &rightTaxon,
        const MafReadOptions &options = MafReadOptions());

};

}  // namespace trio
}  // namespace anchorwave
