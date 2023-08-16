/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

#pragma once

#include "InputParser.h"
#include "./impl/getSequencesFromGff.h"
#include "./impl/deNovoGenomeVariantCalling.h"
#include "./model/AlignmentMatch.h"
#include "./myImportandFunction/alignSlidingWindow.h"
#include "./service/TransferGffWithNucmerResult.h"
#include "./version.h"

#include <iostream>
#include <sstream>
#include <vector>

int gff2seq(int argc, char **argv);

int genomeAlignment(int argc, char **argv);

int proportionalAlignment(int argc, char **argv);

int ali(int argc, char **argv);

int geno(int argc, char **argv);

int pro(int argc, char **argv);
