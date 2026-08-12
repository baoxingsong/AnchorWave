#include "CompressedInput.h"

#include <zlib.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h>

namespace anchorwave {
namespace io {
namespace {

std::runtime_error inputError(
        const std::string &path, const std::string &message) {
    return std::runtime_error(
            "cannot read input file '" + path + "': " + message);
}

std::string temporaryDirectory() {
    const char *configured = std::getenv("TMPDIR");
    if (configured != nullptr && configured[0] != '\0') {
        return configured;
    }
    return "/tmp";
}

void writeAll(int descriptor, const unsigned char *data, std::size_t size,
              const std::string &sourcePath) {
    std::size_t written = 0;
    while (written < size) {
        const ssize_t count = ::write(
                descriptor, data + written, size - written);
        if (count < 0) {
            if (errno == EINTR) {
                continue;
            }
            throw inputError(
                    sourcePath,
                    "cannot write decompressed temporary file: " +
                    std::string(std::strerror(errno)));
        }
        if (count == 0) {
            throw inputError(
                    sourcePath,
                    "cannot write decompressed temporary file");
        }
        written += static_cast<std::size_t>(count);
    }
}

std::string decompressToTemporaryFile(const std::string &path) {
    gzFile input = gzopen(path.c_str(), "rb");
    if (input == nullptr) {
        throw inputError(path, "zlib could not open the gzip stream");
    }

    std::string pattern = temporaryDirectory();
    if (pattern.empty() || pattern.back() != '/') {
        pattern.push_back('/');
    }
    pattern += "anchorwave-input-XXXXXX";
    std::vector<char> temporaryPath(pattern.begin(), pattern.end());
    temporaryPath.push_back('\0');
    const int output = mkstemp(temporaryPath.data());
    if (output < 0) {
        const std::string error = std::strerror(errno);
        gzclose(input);
        throw inputError(
                path, "cannot create a temporary file: " + error);
    }

    const std::string outputPath(temporaryPath.data());
    bool complete = false;
    try {
        std::vector<unsigned char> buffer(1024 * 1024);
        while (true) {
            const int count = gzread(
                    input, buffer.data(),
                    static_cast<unsigned int>(buffer.size()));
            if (count > 0) {
                writeAll(output, buffer.data(),
                         static_cast<std::size_t>(count), path);
                continue;
            }
            if (count < 0) {
                int errorNumber = Z_OK;
                const char *errorText = gzerror(input, &errorNumber);
                throw inputError(
                        path, errorText == nullptr
                                      ? "gzip decompression failed"
                                      : errorText);
            }
            break;
        }

        const int gzipCloseStatus = gzclose(input);
        input = nullptr;
        if (gzipCloseStatus != Z_OK) {
            throw inputError(path, "gzip integrity check failed");
        }
        if (::close(output) != 0) {
            throw inputError(
                    path, "cannot close decompressed temporary file: " +
                          std::string(std::strerror(errno)));
        }
        complete = true;
    } catch (...) {
        if (input != nullptr) {
            gzclose(input);
        }
        ::close(output);
        std::remove(outputPath.c_str());
        throw;
    }

    if (!complete) {
        std::remove(outputPath.c_str());
        throw inputError(path, "gzip decompression did not complete");
    }
    return outputPath;
}

class MaterializedInputRegistry {
public:
    ~MaterializedInputRegistry() {
        for (const std::string &path : temporaryFiles_) {
            std::remove(path.c_str());
        }
    }

    std::string resolve(const std::string &path) {
        if (!isGzipCompressed(path)) {
            return path;
        }

        // Input setup is normally serial. Keeping the lock during
        // decompression also prevents concurrent duplicate materializations.
        std::lock_guard<std::mutex> lock(mutex_);
        const auto cached = materialized_.find(path);
        if (cached != materialized_.end()) {
            return cached->second;
        }

        const std::string temporaryPath = decompressToTemporaryFile(path);
        materialized_[path] = temporaryPath;
        temporaryFiles_.push_back(temporaryPath);
        std::cerr << "AnchorWave input: detected gzip/BGZF stream '"
                  << path << "'; using a managed temporary file" << std::endl;
        return temporaryPath;
    }

private:
    std::mutex mutex_;
    std::unordered_map<std::string, std::string> materialized_;
    std::vector<std::string> temporaryFiles_;
};

MaterializedInputRegistry &registry() {
    static MaterializedInputRegistry instance;
    return instance;
}

}  // namespace

bool isGzipCompressed(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input.good()) {
        throw inputError(path, "cannot open file");
    }

    unsigned char signature[2] = {0, 0};
    input.read(reinterpret_cast<char *>(signature), sizeof(signature));
    return input.gcount() == static_cast<std::streamsize>(sizeof(signature)) &&
           signature[0] == 0x1f && signature[1] == 0x8b;
}

std::string materializeInputFile(const std::string &path) {
    return registry().resolve(path);
}

}  // namespace io
}  // namespace anchorwave
