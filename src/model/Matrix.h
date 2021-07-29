//
// Created by bs674 on 8/2/19.
//

#ifndef AND_CNS_MATRIX_H
#define AND_CNS_MATRIX_H

#include <string>

class Matrix {
    private:
        uint8_t ** T;
        int32_t length1;
        int32_t length2;
    public:
        Matrix(int32_t _length1, int32_t _length2);
        ~Matrix();
        void reset(int32_t a, int32_t b);
        void set(int32_t & a, int32_t & b, uint8_t & c);
        uint8_t & get(int32_t & a, int32_t & b);
};

#endif //AND_CNS_MATRIX_H
