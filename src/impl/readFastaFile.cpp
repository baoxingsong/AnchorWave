//
// Created by baoxing on 10/10/17.
//

#include "readFastaFile.h"


void readFastaFile( const std::string& filePath, std::map<std::string, Fasta>& sequences ){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }
    std::string name="";
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if( name.size()>0 ){
                std::string sequence = sequencestream.str();
                transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
                Fasta f (name, sequence);
                sequences[name]=f;
            }
            name=line.substr(1, line.find(" ", 0)-1);
            if( name[0] == '>' ){
                name=name.substr(1, name.find("\t", 0)-1);
            }else{
                name=name.substr(0, name.find("\t", 0)-1);
            }
            sequencestream.str(std::string());
        }else{
            sequencestream << line;
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        transform(sequence.begin(), sequence.end(), sequence.begin(),::toupper);
        Fasta f (name, sequence);
        if( sequence.size()>0 ){
            sequences[name]=f;
        }
    }
}
