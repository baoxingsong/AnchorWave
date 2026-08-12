#include "src/impl/readFastaFile.h"
#include "src/impl/readGffFile.h"
#include "src/impl/getSubsequence.h"
#include "src/io/CompressedInput.h"

#include "gtest/gtest.h"

#include <zlib.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <unistd.h>

namespace {

class TemporaryFiles {
public:
    ~TemporaryFiles() {
        for (const std::string &path : paths_) {
            std::remove(path.c_str());
        }
    }

    std::string create(const std::string &suffix = std::string()) {
        std::string pattern = "/tmp/anchorwave-compressed-test-XXXXXX";
        std::vector<char> writable(pattern.begin(), pattern.end());
        writable.push_back('\0');
        const int descriptor = mkstemp(writable.data());
        if (descriptor < 0) {
            throw std::runtime_error("mkstemp failed");
        }
        close(descriptor);

        std::string path(writable.data());
        if (!suffix.empty()) {
            const std::string renamed = path + suffix;
            if (std::rename(path.c_str(), renamed.c_str()) != 0) {
                std::remove(path.c_str());
                throw std::runtime_error("temporary rename failed");
            }
            path = renamed;
        }
        paths_.push_back(path);
        return path;
    }

private:
    std::vector<std::string> paths_;
};

void writePlain(const std::string &path, const std::string &content) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    ASSERT_TRUE(output.good());
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
    ASSERT_TRUE(output.good());
}

void writeGzip(const std::string &path, const std::string &content,
               const char *mode = "wb") {
    gzFile output = gzopen(path.c_str(), mode);
    ASSERT_NE(nullptr, output);
    ASSERT_EQ(static_cast<int>(content.size()),
              gzwrite(output, content.data(),
                      static_cast<unsigned int>(content.size())));
    ASSERT_EQ(Z_OK, gzclose(output));
}

std::string readAll(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    return std::string(
            (std::istreambuf_iterator<char>(input)),
            std::istreambuf_iterator<char>());
}

TEST(CompressedInputTest, DetectsCompressionByMagicNotFilename) {
    TemporaryFiles files;
    const std::string plainWithGzipSuffix = files.create(".gz");
    writePlain(plainWithGzipSuffix, "plain text\n");
    EXPECT_FALSE(anchorwave::io::isGzipCompressed(plainWithGzipSuffix));
    EXPECT_EQ(plainWithGzipSuffix,
              anchorwave::io::materializeInputFile(plainWithGzipSuffix));

    const std::string gzipWithoutSuffix = files.create(".fasta-data");
    writeGzip(gzipWithoutSuffix, ">chr1\nACGT\n");
    EXPECT_TRUE(anchorwave::io::isGzipCompressed(gzipWithoutSuffix));
    const std::string materialized =
            anchorwave::io::materializeInputFile(gzipWithoutSuffix);
    EXPECT_NE(gzipWithoutSuffix, materialized);
    EXPECT_EQ(">chr1\nACGT\n", readAll(materialized));
    EXPECT_EQ(materialized,
              anchorwave::io::materializeInputFile(gzipWithoutSuffix));
}

TEST(CompressedInputTest, ReadsConcatenatedGzipMembersLikeBgzfStreams) {
    TemporaryFiles files;
    const std::string path = files.create(".bgz");
    writeGzip(path, "first\n");
    writeGzip(path, "second\n", "ab");

    const std::string materialized =
            anchorwave::io::materializeInputFile(path);
    EXPECT_EQ("first\nsecond\n", readAll(materialized));
}

TEST(CompressedInputTest, RejectsAFileWithCorruptGzipPayload) {
    TemporaryFiles files;
    const std::string path = files.create(".broken");
    const char corrupt[] = {static_cast<char>(0x1f),
                            static_cast<char>(0x8b), 'x', 'x', 'x'};
    writePlain(path, std::string(corrupt, sizeof(corrupt)));
    EXPECT_TRUE(anchorwave::io::isGzipCompressed(path));
    EXPECT_THROW(anchorwave::io::materializeInputFile(path),
                 std::runtime_error);
}

TEST(CompressedInputTest, FastaAndGffReadersAcceptCompressedInput) {
    TemporaryFiles files;
    const std::string fasta = files.create(".payload");
    writeGzip(fasta, ">chr1 description\nACGT\nTGCA\n");

    std::map<std::string, std::tuple<std::string, long, long, int>> sequences;
    readFastaFile(fasta, sequences);
    ASSERT_EQ(1U, sequences.size());
    ASSERT_EQ(1U, sequences.count("chr1"));
    EXPECT_EQ(8, std::get<1>(sequences.at("chr1")));
    EXPECT_EQ(4, std::get<3>(sequences.at("chr1")));
    EXPECT_NE(fasta, std::get<0>(sequences.at("chr1")));
    EXPECT_EQ("CGTTGC", getSubsequence2(sequences, "chr1", 2, 7));

    const std::string gff = files.create(".annotation");
    writeGzip(
            gff,
            "##gff-version 3\n"
            "chr1\ttest\tCDS\t2\t5\t.\t+\t0\tID=cds1;Parent=tx1\n");
    std::map<std::string, std::vector<Transcript>> transcripts;
    readGffFile1(
            gff, transcripts,
            "([\\s\\S]*)Parent=([a-zA-Z0-9.:_-]+)", 0);
    ASSERT_EQ(1U, transcripts.count("chr1"));
    ASSERT_EQ(1U, transcripts.at("chr1").size());
    EXPECT_EQ("tx1", transcripts.at("chr1")[0].getName());
    EXPECT_EQ(2, transcripts.at("chr1")[0].getPStart());
    EXPECT_EQ(5, transcripts.at("chr1")[0].getPEnd());
}

}  // namespace
