//
// Created by bs674 on 12/19/20.
//

#ifndef PROALI_EVALUATION_H
#define PROALI_EVALUATION_H


#include <iostream>
#include <sstream>
#include "../model/model.h"
#include "../util/util.h"
#include "../impl/impl.h"
#include "../service/service.h"


void vcfToVariant( std::string & vcfFilePath,  std::map< std::string, std::vector<Variant>> & variants , std::string & chrTotest);
bool equalVariant( Variant & v1, Variant & v2,  std::map<std::string, std::string>  & referenceGenome );
bool variantWith( Variant & v, std::map< std::string, std::vector<Variant>> variants, std::map<std::string, std::string>  & referenceGenome );

void gffToVariant(std::string & fastaFilePath,  std::string & gffFile, std::string & chrTotest, std::vector<Variant> & sdiRecordsThisOne );

#endif //PROALI_EVALUATION_H
