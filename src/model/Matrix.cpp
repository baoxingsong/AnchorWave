//
// Created by bs674 on 8/2/19.
//

#include "Matrix.h"

Matrix::Matrix(int32_t _length1, int32_t _length2){
    length1 = _length1;
    length2=_length2;
    T = new uint8_t *[length1];
    for (int32_t i = 0; i < (length1); ++i) {
        T[i] = new uint8_t[length2];
        std::fill_n(T[i], length2, 0);
    }
}

Matrix::~Matrix(){
    for ( int i=0; i<length1; ++i ){
        delete[] T[i];
    }
    delete[] T;
}

void Matrix::reset(int32_t a, int32_t b){
    for (int i = 0; i < a; ++i) {
        std::fill_n(T[i], b, 255);
    }
}

void Matrix::set(int32_t & a, int32_t & b, uint8_t & c){
    T[a][b] =c;
}

uint8_t & Matrix::get(int32_t & a, int32_t & b){
    return T[a][b];
}
