#include "AlignmentBlockBuffer.h"

#include <array>
#include <cerrno>
#include <cstring>
#include <ostream>
#include <stdexcept>

namespace anchorwave {

bool ungappedSequenceEquals(const std::string &alignment,
                            const std::string &source) {
    std::size_t sourceOffset = 0;
    for (char base : alignment) {
        if (base == '-') {
            continue;
        }
        if (sourceOffset == source.size() || source[sourceOffset] != base) {
            return false;
        }
        ++sourceOffset;
    }
    return sourceOffset == source.size();
}

AlignmentBlockBuffer::AlignmentBlockBuffer(bool spoolRows)
        : spoolRows_(spoolRows),
          referenceSpool_(nullptr),
          querySpool_(nullptr),
          referenceLength_(0),
          queryLength_(0),
          alignmentColumns_(0) {
    if (spoolRows_) {
        referenceSpool_ = openSpool();
        try {
            querySpool_ = openSpool();
        } catch (...) {
            std::fclose(referenceSpool_);
            referenceSpool_ = nullptr;
            throw;
        }
    }
}

AlignmentBlockBuffer::~AlignmentBlockBuffer() {
    closeSpools();
}

void AlignmentBlockBuffer::append(
        const std::string &referenceAlignment,
        const std::string &queryAlignment,
        const std::string &referenceSource,
        const std::string &querySource) {
    if (referenceAlignment.size() != queryAlignment.size()) {
        throw std::runtime_error("alignment rows have different lengths");
    }
    if (!ungappedSequenceEquals(referenceAlignment, referenceSource) ||
        !ungappedSequenceEquals(queryAlignment, querySource)) {
        throw std::runtime_error(
                "aligned row does not reconstruct its source sequence");
    }

    referenceLength_ += referenceSource.size();
    queryLength_ += querySource.size();
    alignmentColumns_ += referenceAlignment.size();
    if (spoolRows_) {
        appendTo(referenceSpool_, referenceAlignment);
        appendTo(querySpool_, queryAlignment);
    }
}

void AlignmentBlockBuffer::reset() {
    closeSpools();
    referenceLength_ = 0;
    queryLength_ = 0;
    alignmentColumns_ = 0;
    if (spoolRows_) {
        referenceSpool_ = openSpool();
        try {
            querySpool_ = openSpool();
        } catch (...) {
            std::fclose(referenceSpool_);
            referenceSpool_ = nullptr;
            throw;
        }
    }
}

bool AlignmentBlockBuffer::empty() const {
    return alignmentColumns_ == 0;
}

std::size_t AlignmentBlockBuffer::referenceLength() const {
    return referenceLength_;
}

std::size_t AlignmentBlockBuffer::queryLength() const {
    return queryLength_;
}

std::size_t AlignmentBlockBuffer::alignmentColumns() const {
    return alignmentColumns_;
}

void AlignmentBlockBuffer::writeReference(std::ostream &output) {
    if (!spoolRows_) {
        throw std::logic_error("alignment rows were not spooled");
    }
    copyTo(referenceSpool_, output);
}

void AlignmentBlockBuffer::writeQuery(std::ostream &output) {
    if (!spoolRows_) {
        throw std::logic_error("alignment rows were not spooled");
    }
    copyTo(querySpool_, output);
}

std::FILE *AlignmentBlockBuffer::openSpool() {
    std::FILE *file = std::tmpfile();
    if (file == nullptr) {
        throw std::runtime_error(
                std::string("cannot create temporary alignment spool: ") +
                std::strerror(errno));
    }
    return file;
}

void AlignmentBlockBuffer::appendTo(std::FILE *file,
                                    const std::string &text) {
    if (!text.empty() &&
        std::fwrite(text.data(), 1, text.size(), file) != text.size()) {
        throw std::runtime_error(
                std::string("cannot write temporary alignment spool: ") +
                std::strerror(errno));
    }
}

void AlignmentBlockBuffer::copyTo(std::FILE *file, std::ostream &output) {
    if (std::fflush(file) != 0 || std::fseek(file, 0, SEEK_SET) != 0) {
        throw std::runtime_error(
                std::string("cannot rewind temporary alignment spool: ") +
                std::strerror(errno));
    }

    // Keep this below conservative worker-stack limits (notably macOS) while
    // still copying the spool in reasonably large sequential chunks.
    std::array<char, 64 * 1024> buffer;
    while (true) {
        const std::size_t count =
                std::fread(buffer.data(), 1, buffer.size(), file);
        if (count != 0) {
            output.write(buffer.data(), static_cast<std::streamsize>(count));
            if (!output.good()) {
                throw std::runtime_error("cannot write whole-MAF output");
            }
        }
        if (count != buffer.size()) {
            if (std::ferror(file)) {
                throw std::runtime_error(
                        std::string("cannot read temporary alignment spool: ") +
                        std::strerror(errno));
            }
            break;
        }
    }
}

void AlignmentBlockBuffer::closeSpools() {
    if (referenceSpool_ != nullptr) {
        std::fclose(referenceSpool_);
        referenceSpool_ = nullptr;
    }
    if (querySpool_ != nullptr) {
        std::fclose(querySpool_);
        querySpool_ = nullptr;
    }
}

}  // namespace anchorwave
