#include "src/trio/io/MafEvidenceReader.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

std::string location(const std::string &source, std::size_t line) {
    std::ostringstream out;
    out << source << ':' << line;
    return out.str();
}

std::vector<std::string> fields(const std::string &line) {
    std::istringstream input(line);
    std::vector<std::string> result;
    std::string value;
    while (input >> value) {
        result.push_back(value);
    }
    return result;
}

int64_t parseInt64(const std::string &value, const std::string &where,
                   const std::string &fieldName) {
    if (value.empty()) {
        throw MafFormatError(where + ": empty " + fieldName);
    }
    errno = 0;
    char *end = nullptr;
    const long long parsed = std::strtoll(value.c_str(), &end, 10);
    if (errno == ERANGE || end == value.c_str() || *end != '\0') {
        throw MafFormatError(where + ": invalid " + fieldName + " '" + value + "'");
    }
    if (parsed < std::numeric_limits<int64_t>::min() ||
        parsed > std::numeric_limits<int64_t>::max()) {
        throw MafFormatError(where + ": out-of-range " + fieldName + " '" + value + "'");
    }
    return static_cast<int64_t>(parsed);
}

bool isDnaCharacter(char value) {
    const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
    switch (upper) {
        case 'A': case 'C': case 'G': case 'T': case 'U':
        case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
        case 'B': case 'D': case 'H': case 'V': case 'N': case 'X':
        case '.': case '*': case '-':
            return true;
        default:
            return false;
    }
}

void validateRow(const MafRow &row, const std::string &source,
                 const MafReadOptions &options) {
    const std::string where = location(source, row.lineNumber);
    if (row.source.empty()) {
        throw MafFormatError(where + ": empty source name");
    }
    if (row.start0 < 0 || row.size < 0 || row.sourceSize < 0) {
        throw MafFormatError(where + ": MAF coordinates and sizes must be non-negative");
    }
    if (row.start0 > row.sourceSize || row.size > row.sourceSize - row.start0) {
        throw MafFormatError(where + ": interval exceeds source size");
    }
    if (row.strand != '+' && row.strand != '-') {
        throw MafFormatError(where + ": strand must be '+' or '-'");
    }
    int64_t bases = 0;
    for (char value : row.text) {
        if (value != '-') {
            ++bases;
        }
        if (options.requireDnaAlphabet && !isDnaCharacter(value)) {
            std::ostringstream message;
            message << where << ": invalid DNA/MAF character '" << value << "'";
            throw MafFormatError(message.str());
        }
    }
    if (bases != row.size) {
        std::ostringstream message;
        message << where << ": row declares size " << row.size
                << " but alignment text contains " << bases << " non-gap bases";
        throw MafFormatError(message.str());
    }
}

void finalizeBlock(MafBlock &block, std::vector<MafBlock> &blocks,
                   const std::string &source, const std::string &leftTaxon,
                   const std::string &rightTaxon, const MafReadOptions &options) {
    if (block.lineNumber == 0) {
        return;
    }
    if (block.rows.empty() && !options.allowEmptyBlocks) {
        throw MafFormatError(location(source, block.lineNumber) +
                             ": alignment block contains no sequence rows");
    }
    if (options.requireExactlyTwoRows && block.rows.size() != 2) {
        std::ostringstream message;
        message << location(source, block.lineNumber)
                << ": expected exactly two sequence rows, found " << block.rows.size();
        throw MafFormatError(message.str());
    }
    if (!block.rows.empty()) {
        const std::size_t width = block.rows.front().text.size();
        for (const MafRow &row : block.rows) {
            if (row.text.size() != width) {
                throw MafFormatError(location(source, row.lineNumber) +
                                     ": rows in one MAF block have different alignment widths");
            }
        }
    }
    if (block.rows.size() >= 1) {
        block.rows[0].taxon = leftTaxon;
    }
    if (block.rows.size() >= 2) {
        block.rows[1].taxon = rightTaxon;
    }
    block.blockIndex = blocks.size();
    blocks.push_back(block);
    block = MafBlock();
}

}  // namespace

int64_t MafRow::forwardStart0() const {
    return strand == '+' ? start0 : sourceSize - start0 - size;
}

int64_t MafRow::forwardEnd0() const {
    return forwardStart0() + size;
}

int64_t MafRow::columnForwardPosition0(std::size_t column) const {
    if (column >= text.size()) {
        throw std::out_of_range("MAF column is outside the row text");
    }
    if (text[column] == '-') {
        return -1;
    }
    int64_t precedingBases = 0;
    for (std::size_t i = 0; i < column; ++i) {
        if (text[i] != '-') {
            ++precedingBases;
        }
    }
    if (strand == '+') {
        return start0 + precedingBases;
    }
    return sourceSize - start0 - 1 - precedingBases;
}

std::vector<MafBlock> MafEvidenceReader::read(
    std::istream &input, const std::string &sourceName,
    const std::string &leftTaxon, const std::string &rightTaxon,
    const MafReadOptions &options) {
    if (leftTaxon.empty() || rightTaxon.empty()) {
        throw MafFormatError(sourceName + ": pairwise taxa must not be empty");
    }
    if (leftTaxon == rightTaxon) {
        throw MafFormatError(sourceName + ": pairwise taxa must be different");
    }

    std::vector<MafBlock> blocks;
    MafBlock current;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        const std::vector<std::string> tokens = fields(line);
        if (tokens.empty()) {
            finalizeBlock(current, blocks, sourceName, leftTaxon, rightTaxon, options);
            continue;
        }
        if (tokens[0][0] == '#') {
            continue;
        }
        if (tokens[0] == "a") {
            finalizeBlock(current, blocks, sourceName, leftTaxon, rightTaxon, options);
            current.lineNumber = lineNumber;
            for (std::size_t i = 1; i < tokens.size(); ++i) {
                const std::size_t equals = tokens[i].find('=');
                if (equals == std::string::npos || equals == 0) {
                    throw MafFormatError(location(sourceName, lineNumber) +
                                         ": malformed alignment attribute '" + tokens[i] + "'");
                }
                current.attributes[tokens[i].substr(0, equals)] = tokens[i].substr(equals + 1);
            }
            continue;
        }
        if (tokens[0] == "s") {
            if (current.lineNumber == 0) {
                throw MafFormatError(location(sourceName, lineNumber) +
                                     ": sequence row appears before an 'a' block line");
            }
            if (tokens.size() != 7) {
                std::ostringstream message;
                message << location(sourceName, lineNumber)
                        << ": expected 7 fields in sequence row, found " << tokens.size();
                throw MafFormatError(message.str());
            }
            MafRow row;
            row.source = tokens[1];
            row.start0 = parseInt64(tokens[2], location(sourceName, lineNumber), "start");
            row.size = parseInt64(tokens[3], location(sourceName, lineNumber), "size");
            if (tokens[4].size() != 1) {
                throw MafFormatError(location(sourceName, lineNumber) + ": malformed strand");
            }
            row.strand = tokens[4][0];
            row.sourceSize = parseInt64(tokens[5], location(sourceName, lineNumber), "source size");
            row.text = tokens[6];
            row.lineNumber = lineNumber;
            validateRow(row, sourceName, options);
            current.rows.push_back(row);
            continue;
        }
        throw MafFormatError(location(sourceName, lineNumber) +
                             ": unsupported MAF record type '" + tokens[0] + "'");
    }
    finalizeBlock(current, blocks, sourceName, leftTaxon, rightTaxon, options);
    if (!input.eof() && input.fail()) {
        throw MafFormatError(sourceName + ": I/O error while reading MAF");
    }
    return blocks;
}

std::vector<MafBlock> MafEvidenceReader::readFile(
    const std::string &path, const std::string &leftTaxon,
    const std::string &rightTaxon, const MafReadOptions &options) {
    std::ifstream input(path.c_str());
    if (!input) {
        throw MafFormatError(path + ": unable to open MAF file");
    }
    return read(input, path, leftTaxon, rightTaxon, options);
}

}  // namespace trio
}  // namespace anchorwave
