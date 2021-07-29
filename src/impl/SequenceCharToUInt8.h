//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_SEQUENCECHARTOUINT8_H
#define SONG_CNS_SEQUENCECHARTOUINT8_H

#include <string>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <iostream>
#include <sstream>
class SequenceCharToUInt8 {
    private:
        std::map<char, int8_t> dna_acid_int_map;
        std::map<int8_t, char > dna_int_acid_map;
    public:
        SequenceCharToUInt8();
        int8_t get_dna_acid_int_map( const char & c );
        char get_dna_int_acid_map( const int8_t & c );
        int8_t * seq_to_int8( std::string & seq );
        std::string int8_to_seq( int8_t * c,  const uint32_t & length );
        int8_t * rev_comp( const int8_t * c, const uint32_t & length );
};


#endif //SONG_CNS_SEQUENCECHARTOUINT8_H
