//
// Created by Baoxing song on 2019-01-02.
//

#include "SequenceCharToUInt8.h"

SequenceCharToUInt8::SequenceCharToUInt8( ){
    dna_acid_int_map['A'] = (int8_t)0;
    dna_acid_int_map['a'] = (int8_t)0;
    dna_acid_int_map['C'] = (int8_t)1;
    dna_acid_int_map['c'] = (int8_t)1;

    dna_acid_int_map['G'] = (int8_t)3;
    dna_acid_int_map['g'] = (int8_t)3;
    dna_acid_int_map['T'] = (int8_t)4;
    dna_acid_int_map['t'] = (int8_t)4;
    dna_acid_int_map['U'] = (int8_t)4;
    dna_acid_int_map['u'] = (int8_t)4;

    dna_int_acid_map[0]='A';
    dna_int_acid_map[1]='C';
    dna_int_acid_map[2]='N';
    dna_int_acid_map[3]='G';
    dna_int_acid_map[4]='T';
}

int8_t SequenceCharToUInt8::get_dna_acid_int_map( const char & c ){
    if( dna_acid_int_map.find(c) == dna_acid_int_map.end() ){
        return (int8_t)2;
    }else{
        return dna_acid_int_map[c];
    }
}

char SequenceCharToUInt8::get_dna_int_acid_map( const int8_t & c ){
    if( dna_int_acid_map.find(c) == dna_int_acid_map.end() ){
        return 'N';
    }else{
        return dna_int_acid_map[c];
    }
}

int8_t * SequenceCharToUInt8::seq_to_int8( std::string & seq ){
    int8_t * seq_int8 = new int8_t[seq.length()];
    for( int i=0; i < seq.length(); ++i ){
        seq_int8[i] = this->get_dna_acid_int_map(seq[i]);
    }
    return seq_int8;
}

std::string SequenceCharToUInt8::int8_to_seq( int8_t * c,  const uint32_t & length ){
    std::string seq;
    for( int i=0; i < length; ++i ){
        seq+=get_dna_int_acid_map(c[i]);
        //std::cout << get_dna_int_acid_map(c[i]) << std::endl;
    }
    return seq;
}
int8_t * SequenceCharToUInt8::rev_comp( const int8_t * c, const uint32_t & length ){
    int8_t * seq_int8 = new int8_t[length];
    for( int i=0; i < length; ++i ){
        seq_int8[i] = (int8_t)4 - c[length - i - 1];
    }
    return seq_int8;
}
