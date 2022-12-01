#ifndef __DYNMATRIX
#define __DYNMATRIX

#include <string>
#include <vector>

struct DistanceNode {
    double distance;
    DistanceNode *nextInColumn;
    DistanceNode *nextInRow;
    // you can add a constructor
};

struct ClusterNode {
    std::string name;
    // you can add more data member
    int numClusters;
    ClusterNode *prev;
    ClusterNode *next;
    DistanceNode *row;
    DistanceNode *column;
    // you can add a constructor
};

void addCluster(ClusterNode *&head, ClusterNode *&tail, const std::string &name);
// adds a cluster (at the tail) and the corresponding row and column to data structure
// distance of all added DistanceNodes should be initialized to 0.0
// at the end, tail should point to the newly added ClusterNode

void removeCluster(ClusterNode *&head, ClusterNode *&tail, ClusterNode *toBeRemoved);
// removes a cluster pointed to by toBeRemoved and the corresponding row and column
// if toBeRemoved is the first or last cluster then head or tail needs to be updated

void findMinimum(ClusterNode *head, ClusterNode *&C, ClusterNode *&D);
// finds the minimum distance (between two different clusters) in the data structure 
// and returns the two clusters via C and D

void UPGMA(ClusterNode *&head, ClusterNode *&tail);
// Implements UPGMA Algorithm and calls combineCluster and useFormula

void combineCluster(ClusterNode *&head, ClusterNode *&tail, ClusterNode *&C, ClusterNode *&D, std::vector<double> values);

// Helper function to implement UPGMA Algorithm by removing and adding ClusterNode
std::vector<double> useFormula(ClusterNode *clusterOne, ClusterNode *clusterTwo, std::vector<double> firstValues, std::vector<double> secondValues);
// Helper function to implement UPGMA Algorithm by computing values using the formula

void printRowByRow(ClusterNode *head);

// Prints DynMatrix row by row to ensure correctly linked nodes
void printColumnByColumn(ClusterNode *head);
// Prints DynMatrix column by column to ensure correctly linked nodes


#endif
