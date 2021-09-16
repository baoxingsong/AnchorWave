//
// Created by bs674 on 12/19/20.
//

#ifndef PROALI_REFORMAT_H
#define PROALI_REFORMAT_H

#include <iostream>
#include <sstream>

#include "../util/util.h"
#include "../model/model.h"
#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include "../service/service.h"


void mafTovcf( std::string & mafFile,  std::string & refGenomeFilePath, std::string & outputovcffile, const bool & gvcf );
void samToMaf (const std::string& samFilePath, const std::string& refGenomeFile,  const std::string& queryGenomeFile, const std::string& output);
void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& vcfFix, std::map<std::string, std::string>& referenceGenome);
void sdiToMaf(std::string & fastaFilePath, std::string & targetFastaFilePath, std::string & sdiFile, std::string & outfile );
#endif //PROALI_REFORMAT_H
