//
// Created by baoxing on 10/10/17.
//

#include "WriteFasta.h"
#include <iostream>

void writeFasta(  std::ostream& out, const std::string& seqname, const std::string& sequence, const int& linewidth) {
    out << ">" << seqname << std::endl;
    std::size_t pos = 0;
    while (pos < sequence.length()) {
        out << sequence.substr(pos, linewidth) << std::endl;
        pos += linewidth;
    }
}
