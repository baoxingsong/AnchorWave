#include "src/model/AlignmentMatch.h"

#include <gtest/gtest.h>

TEST(AlignmentMatchMemory, InternsRepeatedChromosomeAndGeneStrings) {
    AlignmentMatch first("chr1", "chr1", 1, 10, 2, 11, 1.0,
                         POSITIVE, "gene1", "gene1");
    AlignmentMatch second("chr1", "chr1", 20, 30, 21, 31, 1.0,
                          POSITIVE, "gene1", "gene1");

    EXPECT_EQ(&first.getRefChr(), &first.getQueryChr());
    EXPECT_EQ(&first.getRefChr(), &second.getRefChr());
    EXPECT_EQ(&first.getReferenceGeneName(),
              &first.getQueryGeneName());
    EXPECT_EQ(&first.getReferenceGeneName(),
              &second.getReferenceGeneName());
}

TEST(AlignmentMatchMemory, CopiesDoNotDuplicateStrings) {
    AlignmentMatch original("chr2", "chr3", 1, 10, 2, 11, 1.0,
                            NEGATIVE, "gene2", "gene2");
    AlignmentMatch copy = original;

    EXPECT_EQ(&original.getRefChr(), &copy.getRefChr());
    EXPECT_EQ(&original.getReferenceGeneName(),
              &copy.getReferenceGeneName());
    // Four std::string objects alone were commonly 96-128 bytes. The compact
    // representation stores four pointers plus coordinates and scores.
    EXPECT_LT(sizeof(AlignmentMatch), 96u);
}
