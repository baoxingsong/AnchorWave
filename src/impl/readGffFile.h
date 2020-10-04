//
// Created by baoxing on 10/10/17.
//

#ifndef ANNOTATIONLIFTOVER_READGFFFILE_H
#define ANNOTATIONLIFTOVER_READGFFFILE_H

#include "../model/model.h"

void get_transcript_to_gene_map_from_gff (const std::string& filePath, std::map<std::string, std::string >& transcript_to_gene_map);
void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet, const std::string& cdsParentRegex);
void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet);

#endif //ANNOTATIONLIFTOVER_READGFFFILE_H
