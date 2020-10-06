//
// Created by bs674 on 8/9/19.
//

#include "Seed.h"


Seed::Seed(int32_t & _start1, int32_t & _end1, int32_t & _start2, int32_t & _end2){
    start1 = _start1;
    end1 = _end1;
    start2 = _start2;
    end2 = _end2;
}
int32_t Seed::getStart1() const {
    return start1;
}
void Seed::setStart1(int32_t & start1) {
    Seed::start1 = start1;
}
int32_t Seed::getEnd1() const {
    return end1;
}
void Seed::setEnd1(int32_t & end1) {
    Seed::end1 = end1;
}
int32_t Seed::getStart2() const {
    return start2;
}
void Seed::setStart2(int32_t & start2) {
    Seed::start2 = start2;
}
int32_t Seed::getEnd2() const {
    return end2;
}
void Seed::setEnd2(int32_t & end2) {
    Seed::end2 = end2;
}
