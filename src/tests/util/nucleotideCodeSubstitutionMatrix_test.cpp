//
// Created by baoxing on 2/25/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../util/util.h"
#include <string>
#include <iostream>
#include <sstream>
using namespace std;

TEST(nucleotideCodeSubstitutionMatrix_L_test, c1){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    for( int i=0; i<16; ++i ){
        for( int j=0; j<16; ++j ){
            int32_t a1 = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j];
            int32_t a2 = nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[j][i];
            ASSERT_EQ(a1, a2);
            if( i == j && i!=5 ){
                ASSERT_EQ(a1, 36);
            }
        }
    }
    ASSERT_EQ(0, 0);
}
