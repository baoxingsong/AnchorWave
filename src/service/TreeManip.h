//
// Created by bs674 on 6/16/21.
//

#ifndef ANCHORWAVE_TREEMANIP_H
#define ANCHORWAVE_TREEMANIP_H

#include <string>
#include <cstdint>
#include <vector>
#include <iostream>
#include <sstream>
#include "../model/model.h"

class TreeManip {
private:
    PhylogeneticTree::SharedPtr             _tree;
public:
    TreeManip();
    TreeManip(PhylogeneticTree::SharedPtr t);
    ~TreeManip();
    void                        setTree(PhylogeneticTree::SharedPtr t);
    PhylogeneticTree::SharedPtr             getTree();
    double                      calcTreeLength() const;
    unsigned                    countEdges() const;
    void                        scaleAllEdgeLengths(double scaler);
    void                        createTestTree();
    void                        clear();
    typedef std::shared_ptr< TreeManip > SharedPtr;
};


#endif //ANCHORWAVE_TREEMANIP_H
