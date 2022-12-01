#include "STRAND.h"

std::ostream &operator<<(std::ostream &out, const STRAND value) {
    static std::map<STRAND, std::string> strings;
    if (strings.size() == 0) {
        strings[POSITIVE] = "POSITIVE";
        strings[NEGATIVE] = "NEGATIVE";
    }
    return out << strings[value];
}

