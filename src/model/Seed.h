//
// Created by bs674 on 8/9/19.
//

#ifndef AND_CNS_SEED_H
#define AND_CNS_SEED_H



#include <algorithm>


class Seed{
    private:
        int32_t start1;
        int32_t end1;
        int32_t start2;
        int32_t end2;
    public:
        Seed(int32_t & _start1, int32_t & _end1, int32_t & _start2, int32_t & _end2);
        int32_t getStart1() const;
        void setStart1(int32_t & start1);
        int32_t getEnd1() const;
        void setEnd1(int32_t & end1);
        int32_t getStart2() const ;
        void setStart2(int32_t & start2);
        int32_t getEnd2() const ;
        void setEnd2(int32_t & end2);
        bool operator< (const Seed & c) const {
            if( this->start1 == c.start1 && this->start2 == c.start2 ){
                return false;
            }
            if( this->end1 == c.end1 && this->end2 == c.end2 ){
                return false;
            }
            if ( this->start1 < c.start1 ){
                return true;
            }else if( this->start1 == c.start1 ){
                return this->start2 < c.start2;
            }
            return false;
        }
};

#endif //AND_CNS_SEED_H
