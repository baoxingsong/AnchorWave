//
// Created by baoxing on 10/10/17.
//

#include "getSubsequence.h"

std::string getSubsequence( std::map<std::string, Fasta>& sequences,const  std::string& seqName, const int& _start, const int& _end){
//    std::cout << seqName << " : " <<_start << "-" << _end  << std::endl;
    if( sequences.find(seqName)!=sequences.end() ){
        size_t start = _start;
        size_t end = _end;
        if( start > end ){
            size_t temp = start;
            start=end;
            end=temp;
        }
        if( start < 1 ){
            start = 1;
        }
        if( start > sequences[seqName].getSequence().size() ){
            start = sequences[seqName].getSequence().size();
            end = start;
        }else if( end > sequences[seqName].getSequence().size() ){
            end = sequences[seqName].getSequence().size();
        }
        return sequences[seqName].getSequence().substr(start-1, end-start+1);
    }
    return "";
}

std::string getSubsequence( const std::string & sequence, const int& _start, const int& _end){
    size_t start = _start;
    size_t end = _end;
    if( start > end ){
        size_t temp = start;
        start=end;
        end=temp;
    }
    if( start < 1 ){
        start = 1;
    }
    if( start > sequence.size() ){
        start = sequence.size();
        end = start;
    }else if( end > sequence.size() ){
        end = sequence.size();
    }
    return sequence.substr(start-1, end-start+1);
}

std::string getSubsequence( std::map<std::string, Fasta>& sequences,const  std::string& seqName, const int& _start, const int& _end, const STRAND& strand){
    if( sequences.find(seqName)!=sequences.end() ){
        if( strand == POSITIVE ){
            return  getSubsequence( sequences,  seqName,  _start, _end);
        }else{
            std::string seq = getSubsequence( sequences,  seqName,  _start, _end);
            return getReverseComplementary(seq);
        }
    }
    return "";
}
