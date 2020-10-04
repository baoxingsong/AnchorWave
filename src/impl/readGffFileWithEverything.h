//
// Created by Baoxing song on 08.08.18.
//

#ifndef ANNOTATIONLIFTOVER_READGFFFILEWITHEVERYTHING_H
#define ANNOTATIONLIFTOVER_READGFFFILEWITHEVERYTHING_H

#include "../model/model.h"
#include "../util/util.h"

void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, std::vector<std::string> > & geneNameMap,
                                std::map<std::string, Gene > & geneHashMap,
                                std::map<std::string, Transcript> & transcriptHashMap);
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, std::vector<std::string> > & geneNameMap,
                                std::map<std::string, Gene > & geneHashMap,
                                std::map<std::string, Transcript> & transcriptHashMap,
                                std::set<std::string> & cdsParentRegex, std::set<std::string> & transcriptParentRegex,
                                std::set<std::string> & transcriptIdRegex, std::set<std::string> & geneIdRegex, std::set<std::string> & transcriptNickNames,
                                std::set<std::string> & geneNickNames, std::set<std::string> & ignoreTypes);



#endif //ANNOTATIONLIFTOVER_READGFFFILEWITHEVERYTHING_H
