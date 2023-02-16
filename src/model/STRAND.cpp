#include "STRAND.h"

std::ostream &operator<<(std::ostream &out, const STRAND value) {
    if(value == POSITIVE)
        return out << "POSITIVE";

    return out << "NEGATIVE";
}
