//
// Created by bs674 on 6/16/21.
//
// https://stromtutorial.github.io/mac/steps/step-02/01-create-node-class.html

#ifndef ANCHORWAVE_PHYLOGENETICTREENODE_H
#define ANCHORWAVE_PHYLOGENETICTREENODE_H
#include <string>
#include <cstdint>
#include <vector>
#include <iostream>
#include <sstream>

class PhylogeneticTreeNode {
    friend class PhylogeneticTree;
private:
    void clear();
    PhylogeneticTreeNode * _left_child;
    PhylogeneticTreeNode * _right_sib;
    PhylogeneticTreeNode * _parent;
    int _number;
    std::string _name;
    double _edge_length;

public:
    PhylogeneticTreeNode();
    ~PhylogeneticTreeNode();

    PhylogeneticTreeNode * getParent()     {return _parent;}
    PhylogeneticTreeNode * getLeftChild()  {return _left_child;}
    PhylogeneticTreeNode * getRightSib()   {return _right_sib;}
    int getNumber() {return _number;}
    std::string getName() {return _name;}
    double              getEdgeLength() {return _edge_length;}
    void                setEdgeLength(double v);
    static const double _smallest_edge_length;
    typedef std::vector<PhylogeneticTreeNode> Vector;
    typedef std::vector<PhylogeneticTreeNode *> PtrVector;
};

inline PhylogeneticTreeNode::PhylogeneticTreeNode() {
    std::cout << "Creating Node object" << std::endl;
    clear();
}

inline PhylogeneticTreeNode::~PhylogeneticTreeNode() {
    std::cout << "Destroying Node object" << std::endl;
}

inline void PhylogeneticTreeNode::clear() {
    _left_child = 0;
    _right_sib = 0;
    _parent = 0;
    _number = -1;
    _name = "";
    _edge_length = _smallest_edge_length;
}

inline void PhylogeneticTreeNode::setEdgeLength(double v) {
    _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
}

#endif //ANCHORWAVE_PHYLOGENETICTREENODE_H
