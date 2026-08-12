#include "src/trio/io/MafEvidenceReader.h"

#include "gtest/gtest.h"

#include <sstream>

namespace anchorwave {
namespace trio {
namespace {

TEST(MafEvidenceReader, ReadsPositiveAndNegativeRows) {
    std::istringstream input(
        "##maf version=1\n"
        "a score=12 method=wfa2\n"
        "s chr1 2 5 + 20 AC-GTA\n"
        "s ctg7 3 6 - 30 ACGGTA\n\n");
    const std::vector<MafBlock> blocks =
        MafEvidenceReader::read(input, "fixture.maf", "I1", "O1");
    ASSERT_EQ(1u, blocks.size());
    ASSERT_EQ(2u, blocks[0].rows.size());
    EXPECT_EQ("I1", blocks[0].rows[0].taxon);
    EXPECT_EQ("O1", blocks[0].rows[1].taxon);
    EXPECT_EQ(2, blocks[0].rows[0].forwardStart0());
    EXPECT_EQ(21, blocks[0].rows[1].forwardStart0());
    EXPECT_EQ(26, blocks[0].rows[1].columnForwardPosition0(0));
    EXPECT_EQ(25, blocks[0].rows[1].columnForwardPosition0(1));
    EXPECT_EQ(-1, blocks[0].rows[0].columnForwardPosition0(2));
}

TEST(MafEvidenceReader, RejectsDeclaredSizeMismatch) {
    std::istringstream input(
        "a score=1\n"
        "s chr1 0 4 + 10 A-C\n"
        "s chr2 0 3 + 10 ATC\n");
    EXPECT_THROW(MafEvidenceReader::read(input, "bad.maf", "I1", "I2"),
                 MafFormatError);
}

TEST(MafEvidenceReader, RejectsCoordinateOverflow) {
    std::istringstream input(
        "a score=1\n"
        "s chr1 9 2 + 10 AC\n"
        "s chr2 0 2 + 10 AC\n");
    EXPECT_THROW(MafEvidenceReader::read(input, "bad.maf", "I1", "I2"),
                 MafFormatError);
}

TEST(MafEvidenceReader, RejectsDifferentAlignmentWidths) {
    std::istringstream input(
        "a score=1\n"
        "s chr1 0 3 + 10 ACT\n"
        "s chr2 0 3 + 10 A-CT\n");
    EXPECT_THROW(MafEvidenceReader::read(input, "bad.maf", "I1", "I2"),
                 MafFormatError);
}

TEST(MafEvidenceReader, RequiresBlockHeaderAndTwoRows) {
    std::istringstream noHeader("s chr1 0 1 + 1 A\n");
    EXPECT_THROW(MafEvidenceReader::read(noHeader, "bad.maf", "I1", "I2"),
                 MafFormatError);

    std::istringstream oneRow("a score=0\ns chr1 0 1 + 1 A\n");
    EXPECT_THROW(MafEvidenceReader::read(oneRow, "bad.maf", "I1", "I2"),
                 MafFormatError);
}

}  // namespace
}  // namespace trio
}  // namespace anchorwave
