//
// Created by song on 8/8/18.
//

#include "include/gtest/gtest.h"
#include "../../impl/impl.h"
#include <string>
#include <iostream>
#include <sstream>
TEST(sequencesToKmerCountMatrix, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/bs674/seq.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );
    for( std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > >::iterator it = kmerCountMatrix.begin(); it!=kmerCountMatrix.end(); ++it ){
        std::cout << it->first << std::endl;
        for( std::map<std::string/*kmer*/, int16_t >::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2 ){
            std::cout << it2->first << "\t" << it2->second << std::endl;
        }
    }
    ASSERT_EQ(0, 0);
}

TEST(sequencesToKmerDistanceMatrix, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/bs674/seq.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );
    std::vector<std::string> seqNames;
    float ** distanceMatrix = new float * [sequences.size()];

    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }
    kmerCountMatrixToDistanceMatrix(  kmerCountMatrix, seqNames, distanceMatrix );

    for(size_t i = 0; i < seqNames.size(); ++i){
        std::cout << seqNames[i] ;
        for(size_t j = 0; j < seqNames.size(); ++j){
            std::cout << "\t" << distanceMatrix[i][j];
        }
        std::cout<< std::endl;
    }

    for(size_t i = 0; i < sequences.size(); ++i){
        delete[] distanceMatrix[i];
    }
    delete[] distanceMatrix;

    ASSERT_EQ(0, 0);
}


TEST(sequencesToKmeruUpgma, c1){ // just to make sure that every line has been analysed
    std::string fastaFilePath = "/home/bs674/seq.fa";
    int8_t kmer_size = 12;
    std::map<std::string, std::string> sequences;
    readFastaFile( fastaFilePath, sequences);
    std::map<std::string/*sequence name*/,  std::map<std::string/*kmer*/, int16_t > > kmerCountMatrix;
    sequencesToKmerCountMatrix( sequences, kmer_size, kmerCountMatrix );
    std::vector<std::string> seqNames;
    float ** distanceMatrix = new float * [sequences.size()];

    for(size_t i = 0; i < sequences.size(); ++i){
        distanceMatrix[i] = new float [sequences.size()];
    }
    kmerCountMatrixToDistanceMatrix(  kmerCountMatrix, seqNames, distanceMatrix );

    for(size_t i = 0; i < seqNames.size(); ++i){
        std::cout << seqNames[i] ;
        for(size_t j = 0; j < seqNames.size(); ++j){
            std::cout << "\t" << distanceMatrix[i][j];
        }
        std::cout<< std::endl;
    }


    using namespace std;


    // let's start with empty DynMatrix:
    ClusterNode* head = NULL;
    ClusterNode* tail = NULL;

    for (string seqName : seqNames){
        addCluster(head, tail, seqName);

    }
    ClusterNode *node = head;
    for (int i=0; i<seqNames.size(); i++) {
        DistanceNode *newDistance = node->row;
        for (int j=0; j<seqNames.size(); j++) {
            double d = distanceMatrix[i][j];
            newDistance->distance = d;
            newDistance = newDistance->nextInRow;
        }
        node  = node->next;
    }

    // Run UPGMA until there is only one node left
    int j = 1;
    while (head->next) {
        cout << "---------------Printing UPGMA Round " << j << " (by rows)---------------" << endl;
        printRowByRow(head);
        cout << "--------------------------------------------------------------" << endl;

        cout << "---------------Printing UPGMA Round " << j << " (by columns)------------" << endl;
        printColumnByColumn(head);
        cout << "--------------------------------------------------------------" << endl;

        UPGMA(head,tail);

        j++;
    }
    // Print out name of last node
    if (head) {
        cout  << head->name <<endl;

    }
    // BONUS (optional): print the tree in a nice way

    ASSERT_EQ(0, 0);
}
