#include "gtest/gtest.h"

#include <cstdint>

#include "minimap2/ksort.h"

#define anchorwave_radix_key(value) (value)
KRADIX_SORT_INIT(anchorwave_test_64, uint64_t, anchorwave_radix_key, 8)

TEST(MinimapRadixSortTest, AcceptsNullEmptyRange) {
    uint64_t *empty = nullptr;
    EXPECT_NO_THROW(radix_sort_anchorwave_test_64(empty, empty));
}

TEST(MinimapRadixSortTest, SortsSmallRanges) {
    uint64_t one[] = {7};
    radix_sort_anchorwave_test_64(one, one + 1);
    EXPECT_EQ(7U, one[0]);

    uint64_t values[] = {9, 1, 7, 3};
    radix_sort_anchorwave_test_64(values, values + 4);
    EXPECT_EQ(1U, values[0]);
    EXPECT_EQ(3U, values[1]);
    EXPECT_EQ(7U, values[2]);
    EXPECT_EQ(9U, values[3]);
}
