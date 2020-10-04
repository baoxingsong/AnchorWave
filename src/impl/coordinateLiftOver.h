//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_COORDINATELIFTOVER_H
#define ANNOTATIONLIFTOVER_COORDINATELIFTOVER_H

#include "../model/model.h"

int getChangedFromBasement(const std::string & chromosomeName, const int & basement, std::map<std::string, std::vector<Variant> >& variantsMap);
int getBasementFromChanged(const std::string &  chromsomeName, const int & changed, std::map<std::string, std::vector<Variant> >& variantsMap);

#endif //ANNOTATIONLIFTOVER_COORDINATELIFTOVER_H
