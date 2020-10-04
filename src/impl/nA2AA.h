//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_NA2AA_H
#define ANNOTATIONLIFTOVER_NA2AA_H

#include "../model/model.h"
#include "../util/util.h"
std::string nA2AA(std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
void na2aaLocal( std::string& seq, std::stringstream& ssSb, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);


#endif //ANNOTATIONLIFTOVER_NA2AA_H
