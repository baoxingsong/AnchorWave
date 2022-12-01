#include "dynmatrix.h"
#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <iomanip>

double MIN_DISTANCE = 0.0;

void addCluster(ClusterNode *&head, ClusterNode *&tail, const std::string &name)
// adds a cluster (at the tail) and the corresponding row and column to data structure
// distance of all added DistanceNodes should be initialized to 0.0
// at the end, tail should point to the newly added ClusterNode
{
    // Initialize ClusterNode with name and the number of clusters inside
    ClusterNode *newNode = new ClusterNode;
    newNode->name = name;
    newNode->numClusters = 1; // Used to keep track of how many clusters have been combined into one ClusterNode

    if (head == NULL) { // If we are adding a cluster to an empty DynMatrix
        DistanceNode *firstDistance = new DistanceNode;
        firstDistance->distance = 0.0;

        head = tail = newNode;
        head->column = head->row = firstDistance;
    } else {
        // Reattach pointers of the added ClusterNode and set tail to newly created ClusterNode
        newNode->prev = tail;
        tail->next = newNode;
        tail = newNode;
        // Ensure tail->next is set to NULL
        tail->next = NULL;

        // Initialize first DistanceNodes
        DistanceNode *firstCol = new DistanceNode;
        DistanceNode *firstRow = new DistanceNode;
        firstCol->distance = 0.0;
        firstRow->distance = 0.0;

        // Ensure pointers are set to NULL to prevent pointing to wrong locations
        firstCol->nextInColumn = NULL;
        firstCol->nextInRow = NULL;

        firstRow->nextInColumn = NULL;
        firstRow->nextInRow = NULL;

        tail->column = firstCol;
        tail->row = firstRow;

        // Attach first DistanceNodes to first column and row of newly added ClusterNode
        DistanceNode *currCol = tail->column;
        DistanceNode *currRow = tail->row;

        DistanceNode *prevCol = tail->prev->column;
        DistanceNode *prevRow = tail->prev->row;

        prevCol->nextInRow = currCol;
        prevRow->nextInColumn = currRow;

        // Go through each DistanceNode and attach pointers like a singly-linked list
        while (prevCol != prevRow) {
            // Initialize DistanceNodes for use
            DistanceNode *newCol = new DistanceNode;
            DistanceNode *newRow = new DistanceNode;
            newCol->distance = 0.0;
            newRow->distance = 0.0;

            // Ensure pointers are set to NULL to prevent pointing to wrong locations
            newCol->nextInColumn = NULL;
            newCol->nextInRow = NULL;

            newRow->nextInColumn = NULL;
            newRow->nextInRow = NULL;

            currCol->nextInColumn = newCol;
            currRow->nextInRow = newRow;

            currCol = currCol->nextInColumn;
            currRow = currRow->nextInRow;

            prevCol = prevCol->nextInColumn;
            prevRow = prevRow->nextInRow;

            prevCol->nextInRow = currCol;
            prevRow->nextInColumn = currRow;

        }
        // Initialize and attach the last DistanceNode in newly added ClusterNode
        DistanceNode *edgeNode = new DistanceNode;
        edgeNode->distance = 0.0;

        currCol->nextInColumn = currRow->nextInRow = edgeNode;
    }
}

void removeCluster(ClusterNode *&head, ClusterNode *&tail, ClusterNode *toBeRemoved)
// removes a cluster pointed to by toBeRemoved and the corresponding row and column
{

    // If there is only one ClusterNode left in the list, delete the only node and set pointers to NULL for addCluster to function correctly
    if (toBeRemoved == head && toBeRemoved == tail) {
        ClusterNode *tempNode = head;

        tempNode->next = NULL;
        tempNode->prev = NULL;

        head = NULL;
        tail = NULL;
        delete tempNode;
        return;
    }
    if (toBeRemoved == head) {
        // Steps:
        // reassign head to head->next
        // go down ->nextInColumn until NULL, then go down ->nextInRow until NULL and remove pointers

        ClusterNode *tempHead = head;
        ClusterNode *tempHeadRow;
        DistanceNode *tempVerticalDistance;
        DistanceNode *tempHeadCol;
        head = head->next;
        head->prev->next = NULL;
        head->prev = NULL;

        for (tempVerticalDistance = tempHead->column; tempVerticalDistance != NULL; tempVerticalDistance = tempVerticalDistance->nextInColumn) {
            tempVerticalDistance->nextInRow = NULL;
        }

        tempHeadCol = head->column->nextInColumn;

        for (tempHeadRow = head; tempHeadRow != NULL; tempHeadRow = tempHeadRow->next) {
            tempHeadRow->column = tempHeadRow->column->nextInColumn;
            tempHeadRow->row = tempHeadCol;

            tempHeadCol = tempHeadCol->nextInColumn;
        }
        tempHead = NULL;
        return;
    } else if (toBeRemoved == tail) {
        // Steps:
        // go down ->column->nextInColumn until NULL and remove pointers; also set pointers pointing to tail list to NULL
        // go down ->row->nextInRow until NULL and remove pointers; also set pointers pointing to tail list NULL;

        tail = tail->prev;
        tail->next = NULL;
        DistanceNode *tempDistanceRow;
        DistanceNode *tempDistanceCol;

        for (tempDistanceRow = tail->row, tempDistanceCol = tail->column; tempDistanceRow != tempDistanceCol; tempDistanceRow = tempDistanceRow->nextInRow, tempDistanceCol = tempDistanceCol->nextInColumn) {
            tempDistanceCol->nextInRow = NULL;
            tempDistanceRow->nextInColumn = NULL;

        }
        tempDistanceCol->nextInColumn = tempDistanceCol->nextInRow = NULL;


    } else {
        // Steps:
        // reattach pointers to cluster node
        // Go down nextInRow and nextInColumn until pointers are equal and remove + reattach pointers
        // go to ->prev and ->next cluster node and traverse down list
        ClusterNode *prevNode = toBeRemoved->prev;
        ClusterNode *nextNode = toBeRemoved->next;

        toBeRemoved->prev->next = toBeRemoved->next;
        toBeRemoved->next->prev = toBeRemoved->prev;

        DistanceNode *tempColPrev = prevNode->column;
        DistanceNode *tempColNext = nextNode->column;

        DistanceNode *tempRowPrev = prevNode->row;
        DistanceNode *tempRowNext = nextNode->row;

        for (; tempColPrev != NULL; tempColPrev = tempColPrev->nextInColumn, tempColNext = tempColNext->nextInColumn) {
            tempColPrev->nextInRow = tempColNext;
        }

        for (; tempRowPrev != NULL; tempRowPrev = tempRowPrev->nextInRow, tempRowNext = tempRowNext->nextInRow) {
            tempRowPrev->nextInColumn = tempRowNext;
        }
    }

}


void findMinimum(ClusterNode *head, ClusterNode *&C, ClusterNode *&D)
// finds the minimum distance (between two different clusters) in the data structure 
// and returns the two clusters via C and D
{
    // Create temporary Cluster and Distance Nodes
    ClusterNode *curr = NULL;
    DistanceNode *tempStart = NULL;

    ClusterNode *secondCurr = NULL;

    ClusterNode *minClusterNode = NULL;
    DistanceNode *minDistanceNode = NULL;
    // Set minDistance to max Double value
    double minDistance = std::numeric_limits<double>::max();

    // Go through each DistanceNode and check for minimum value
    for (curr = head; curr != NULL; curr = curr->next) {
        for (tempStart = curr->column; tempStart != NULL; tempStart = tempStart->nextInColumn) {
            if (tempStart->distance < minDistance && tempStart->distance != 0) {
                minDistance = tempStart->distance;
                minClusterNode = curr;
                minDistanceNode = tempStart;
            }
        }
    }
    // set corresponding ClusterNode of minimum DistanceNode to given ClusterNode C
    C = minClusterNode;
    // Go through each DistanceNode and find matching minimum value
    for (secondCurr = minClusterNode->next; secondCurr != NULL; secondCurr = secondCurr->next) {
        for (tempStart = secondCurr->column; tempStart != NULL; tempStart = tempStart->nextInColumn) {
            if (tempStart->distance == minDistanceNode->distance) {
                // Set minimum distance value to global variable
                MIN_DISTANCE = minDistanceNode->distance;
                D = secondCurr;
                return;
            }
        }
    }
}

void UPGMA(ClusterNode *&head, ClusterNode *&tail)
// Implement the UPGMA Algorithm
{
    // Create temporary Cluster and Distance Nodes
    ClusterNode *tempHead = head;

    ClusterNode *clusterOne = NULL;
    ClusterNode *clusterTwo = NULL;

    DistanceNode *distanceOne;
    DistanceNode *distanceTwo;
    // Create vectors to hold distances of DistanceNodes
    std::vector<double> clusterOneDistances;
    std::vector<double> clusterTwoDistances;
    std::vector<double> resultantDistances;

    // Call findMinimum function to locate two ClusterNodes containing the minimum distance in DynMatrix
    // clusterOne and clusterTwo are the ClusterNodes containing the minimum distance
    findMinimum(tempHead, clusterOne, clusterTwo);
    std::cout << clusterOne->name << "\tline 247\t" << clusterTwo->name << std::endl;

    // Store the distances in the two ClusterNodes' clusterOne and clusterTwo in vectors
    for (distanceOne = clusterOne->column; distanceOne != NULL; distanceOne = distanceOne->nextInColumn) {
        clusterOneDistances.push_back(distanceOne->distance);
    }
    for (distanceTwo = clusterTwo->column; distanceTwo != NULL; distanceTwo = distanceTwo->nextInColumn) {
        clusterTwoDistances.push_back(distanceTwo->distance);
    }
    // Use useFormula function to compute average distances and store inside resultant vector
    resultantDistances = useFormula(clusterOne, clusterTwo, clusterOneDistances, clusterTwoDistances);

    // Call combineCluster function to remove the two ClusterNodes and add the new ClusterNode into the DynMatrix
    combineCluster(head, tail, clusterOne, clusterTwo, resultantDistances);

}

void combineCluster(ClusterNode *&head, ClusterNode *&tail, ClusterNode *&C, ClusterNode *&D, std::vector<double> values)
// Adds a cluster using a combination of C->name and D->name, then fills the DistanceNodes of the resulting ClusterNode using the values of the given vector parameter
{
    // Set variables to make new ClusterNode name
    std::string leftparen = "(";
    std::string comma = ",";
    std::string rightparen = ")";
    // Assign values of vector parameter to new vector
    std::vector<double> upgmaValues = values;
    // Construct new ClusterNode name
    std::string name = leftparen + C->name + comma + D->name + rightparen;

    DistanceNode *firstDNCol;
    DistanceNode *firstDNRow;
    ClusterNode *tempAddedClusterNode = NULL;
    int i = 0; // initialize vector index
    // Remove the two ClusterNodes containing the minimum distance
    removeCluster(head, tail, C);
    removeCluster(head, tail, D);

    // Add a new cluster to the matrix, with 0.0 distance values for rows and columns.
    addCluster(head, tail, name);
    // Update numClusters of newly added ClusterNode to keep track of how many ClusterNodes are in the newly added one
    tempAddedClusterNode = tail;
    tempAddedClusterNode->numClusters = C->numClusters + D->numClusters;

    // Fill distance values of DistanceNodes in newly added ClusterNode with values from given vector
    for (firstDNCol = tempAddedClusterNode->column, firstDNRow = tempAddedClusterNode->row; firstDNCol != firstDNRow; firstDNCol = firstDNCol->nextInColumn, firstDNRow = firstDNRow->nextInRow) {

        firstDNCol->distance = upgmaValues[i];
        firstDNRow->distance = upgmaValues[i];

        i++;
    }
    // Ensure pointers are NULL to prevent accessing wrongly linked locations
    firstDNCol->nextInColumn = firstDNRow->nextInRow = NULL;
    firstDNCol->nextInRow = firstDNRow->nextInColumn = NULL;

}

std::vector<double> useFormula(ClusterNode *clusterOne, ClusterNode *clusterTwo, std::vector<double> firstValues, std::vector<double> secondValues)
// Uses given formula to calculate values to be put into each DistanceNode of added ClusterNode
{
    // Initialize and set variables
    std::vector<double> result;
    DistanceNode *tempOne;
    DistanceNode *tempTwo;
    int numClusterOne = clusterOne->numClusters;
    int numClusterTwo = clusterTwo->numClusters;
    double numerator = 0.0;
    double denominator = 0.0;

    // Loop through each DistanceNode
    for (tempOne = clusterOne->column, tempTwo = clusterTwo->column; tempOne != NULL; tempOne = tempOne->nextInColumn, tempTwo = tempTwo->nextInColumn) {
        // Only compute values that are not the minimum distance or zero
        if ((tempOne->distance != MIN_DISTANCE && tempOne->distance != 0) && (tempTwo->distance != MIN_DISTANCE && tempTwo->distance != 0)) {
            // Compute numerator and denominator of formula
            numerator = (numClusterOne * tempOne->distance) + (numClusterTwo * tempTwo->distance);
            denominator = numClusterOne + numClusterTwo;
            // Push computed value into vector to be returned
            result.push_back(numerator / denominator);
        }
    }

    return result;
}


void printRowByRow(ClusterNode *head)
// Print DynMatrix row by row for debugging to ensure correctly linked Cluster and Distance Nodes
{
    while (head) {
        std::cout << std::setw(10) << head->name << ":\t";
        DistanceNode *curD = head->row;
        while (curD) {
            std::cout << curD->distance << "\t";
            curD = curD->nextInRow;
        }
        std::cout << std::endl;
        head = head->next;

    }
}

void printColumnByColumn(ClusterNode *head)
// Print DynMatrix column by column for debugging to ensure correctly linked Cluster and Distance Nodes
{
    while (head) {
        std::cout << std::setw(10) << head->name << ":\t";
        DistanceNode *curD = head->column;
        while (curD) {
            std::cout << curD->distance << "\t";
            curD = curD->nextInColumn;
        }
        std::cout << std::endl;
        head = head->next;
    }
}


