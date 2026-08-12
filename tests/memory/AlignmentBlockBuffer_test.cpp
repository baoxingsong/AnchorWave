#include "src/impl/AlignmentBlockBuffer.h"

#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

TEST(AlignmentBlockBuffer, ValidatesWithoutSpoolingRows) {
    anchorwave::AlignmentBlockBuffer buffer(false);
    buffer.append("AC-GT", "A-TGT", "ACGT", "ATGT");

    EXPECT_EQ(4u, buffer.referenceLength());
    EXPECT_EQ(4u, buffer.queryLength());
    EXPECT_EQ(5u, buffer.alignmentColumns());
    EXPECT_FALSE(buffer.empty());
}

TEST(AlignmentBlockBuffer, SpoolsAndResetsRows) {
    anchorwave::AlignmentBlockBuffer buffer(true);
    buffer.append("AC-", "A-T", "AC", "AT");
    buffer.append("G", "G", "G", "G");

    std::ostringstream reference;
    std::ostringstream query;
    buffer.writeReference(reference);
    buffer.writeQuery(query);
    EXPECT_EQ("AC-G", reference.str());
    EXPECT_EQ("A-TG", query.str());

    buffer.reset();
    EXPECT_TRUE(buffer.empty());
    EXPECT_EQ(0u, buffer.referenceLength());
    buffer.append("T", "C", "T", "C");
    std::ostringstream resetReference;
    buffer.writeReference(resetReference);
    EXPECT_EQ("T", resetReference.str());
}

TEST(AlignmentBlockBuffer, RejectsInvalidRows) {
    anchorwave::AlignmentBlockBuffer buffer(false);
    EXPECT_THROW(buffer.append("AC", "A-", "AG", "A"),
                 std::runtime_error);
    EXPECT_THROW(buffer.append("AC", "A", "AC", "A"),
                 std::runtime_error);
}

TEST(AlignmentBlockBuffer, UngappedComparisonDoesNotNeedATemporaryCopy) {
    EXPECT_TRUE(anchorwave::ungappedSequenceEquals("A--CG-T", "ACGT"));
    EXPECT_FALSE(anchorwave::ungappedSequenceEquals("A--CG-A", "ACGT"));
    EXPECT_FALSE(anchorwave::ungappedSequenceEquals("A--CG", "ACGT"));
}
