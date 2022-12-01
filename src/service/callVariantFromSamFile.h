//
// Created by bs674 on 6/5/21.
//

#ifndef ANCHORWAVE_CALLVARIANTFROMSAMFILE_H
#define ANCHORWAVE_CALLVARIANTFROMSAMFILE_H

#include <map>
#include <regex>
#include <cstdlib>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>

#include "../util/util.h"
#include "../model/model.h"
#include "../impl/impl.h"
#include "../myImportandFunction/myImportantFunction.h"
#include "../service/service.h"

void samToVcf(const std::string &samFilePath, const std::string &refGenomeFile, const std::string &queryGenomeFile, const int32_t &range, const std::string &output);

#endif //ANCHORWAVE_CALLVARIANTFROMSAMFILE_H
