//
// Created by bs674 on 6/16/21.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include "../../service/service.h"
#include "../../util/myutil.h"
#include <string>
#include <iostream>
#include <map>
#include <regex>
#include <cstdlib>
#include "../../controlLayer.h"

const double PhylogeneticTreeNode::_smallest_edge_length = 1.0e-12;

TEST(PhylogeneticTree, c1){
    std::cout << "Starting..." << std::endl;
    PhylogeneticTree PhylogeneticTree;
    std::cout << "\nFinished!" << std::endl;

    ASSERT_EQ(0, 0);
}
