//
// Created by song on 8/4/18.
//

#include "callVariantFromSamFile.h"

void samToVcf(const std::string &samFilePath, const std::string &refGenomeFile, const std::string &queryGenomeFile, const int32_t &range, const std::string &output) {
    std::map<std::string, std::string> refGenome;
    readFastaFile(refGenomeFile, refGenome);
    std::map<std::string, std::string> queryGenome;
    readFastaFile(queryGenomeFile, queryGenome);
    std::string line;

    std::string databaseChr;
    int32_t databaseStart;
    int32_t databaseEnd;
    std::string queryChr;
    int32_t queryStart;
    int32_t queryEnd;

    std::map<std::string, std::string> lastChr;
    std::map<std::string, int32_t> lastPosition;

    std::ofstream ovcffile;
    ovcffile.open(output);

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon));
    if (ltm->tm_mday < 10) {
        filedate = filedate + "0" + std::to_string((ltm->tm_mday));
    } else {
        filedate = filedate + std::to_string((ltm->tm_mday));
    }

    std::string refFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refGenomeFile, delim, elems);
    refFileName = elems.back();

    delim = '\t';
    ovcffile << "##fileformat=VCFv4.3" << std::endl;
    ovcffile << "##fileDate=" << filedate << std::endl;
    ovcffile << "##source=proali" << std::endl;
    ovcffile << "##reference=" << refFileName << std::endl;
    ovcffile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
    ovcffile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    ovcffile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
    ovcffile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
    ovcffile << "##FILTER=<ID=q30,Description=\"Quality below 30\">" << std::endl;
    std::string accession = "query";
    accession = songStrReplaceAll(accession, ".fasta", "");
    accession = songStrReplaceAll(accession, ".fa", "");
    accession.erase(std::remove(accession.begin(), accession.end(), ' '), accession.end());
    ovcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << accession << std::endl;


    std::map<std::string, std::map<std::string, std::vector<double>>> geneScores; //fist key is gene name, second key is chr value is a vector of similarity
    std::ifstream infile(samFilePath);
    while (std::getline(infile, line)) { // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if (line[0] != '@') { //ignore the header
            elems.clear();
            split(line, delim, elems);
            queryStart = 1;
//            if(databaseStart>0 && transcriptHashMap.find(elems[0])!=transcriptHashMap.end() ){ // ignore those none mapping records
            queryChr = elems[0];
            if (queryChr.compare("*") != 0) { // ignore those none mapping records
                //std::cout << "begain to analysis " << line << std::endl;
                databaseChr = elems[2];

                databaseStart = std::stoi(elems[3]);
                databaseEnd = databaseStart;
                queryEnd = queryStart;
                std::vector<std::string> cigarElems;
                splitCIGAR(elems[5], cigarElems);
                int headClipping = 0;
                int tailClipping = 0;
                for (int i = 0; i < cigarElems.size(); ++i) {
                    std::string cVal = cigarElems[i];
                    char cLetter = cVal[cVal.length() - 1];
                    int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                    if (i == cigarElems.size() - 1 && (cLetter == 'H' || cLetter == 'S')) { // ignore the last soft/hard clipping
                        tailClipping += cLen;
                        continue;
                    }
                    switch (cLetter) {
                        case 'H':
                            headClipping += cLen;
                            queryStart += cLen;
                            queryEnd += cLen;
                            break;
                        case 'S':
                            headClipping += cLen;
                            queryStart += cLen;
                            queryEnd += cLen;
                            break;
                        case 'M':
                            databaseEnd += cLen;
                            queryEnd += cLen;
                            break;
                        case '=':
                            databaseEnd += cLen;
                            queryEnd += cLen;
                            break;
                        case 'X':
                            databaseEnd += cLen;
                            queryEnd += cLen;
                            break;
                        case 'I':
                            queryEnd += cLen;
                            break;
                        case 'D':
                            databaseEnd += cLen;
                            break;
                        case 'N':
                            databaseEnd += cLen;
                            break;
                        case 'P':
                            std::cout << "current we could not deal with cigar letter P" << std::endl;
                            break;
                        default:
                            std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                            break;
                    }
                }
                --queryEnd;
                --databaseEnd;


                int samFlag = stoi(elems[1]);

                if (0 == samFlag % 32) {
                    int32_t databasePosition = databaseStart;
                    int32_t queryPosition = queryStart;
                    for (int i = 0; i < cigarElems.size(); ++i) {
                        std::string cVal = cigarElems[i];
                        char cLetter = cVal[cVal.length() - 1];
                        int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                        if (i == cigarElems.size() - 1 && (cLetter == 'H' || cLetter == 'S')) { // ignore the last soft/hard clipping
                            tailClipping += cLen;
                            continue;
                        }
                        std::string refSeq;
                        std::string querySeq;
                        std::string preQuerySeq;
                        std::string postQuerySeq;
                        std::string preRefSeq;
                        std::string postRefSeq;

                        switch (cLetter) {
                            case 'M':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);
                                for (int32_t ii = 0; ii < refSeq.size(); ++ii) {
                                    if (refSeq[ii] != querySeq[ii]) {
                                        ovcffile << databaseChr << "\t" << std::to_string(databasePosition + ii) << "\t" << databaseChr << "_" << std::to_string(databasePosition + ii) << "\t" << refSeq[ii] << "\t" << querySeq[ii] << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1"
                                                 << std::endl;
                                    }
                                }
                                databasePosition += cLen;
                                queryPosition += cLen;
                                break;
                            case '=':
                                databasePosition += cLen;
                                queryPosition += cLen;
                                break;
                            case 'X':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);
                                for (int32_t ii = 0; ii < refSeq.size(); ++ii) {
                                    if (refSeq[ii] != querySeq[ii]) {
                                        ovcffile << databaseChr << "\t" << std::to_string(databasePosition + ii) << "\t" << databaseChr << "_" << std::to_string(databasePosition + ii) << "\t" << refSeq[ii] << "\t" << querySeq[ii] << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1"
                                                 << std::endl;
                                    }
                                }
                                databasePosition += cLen;
                                queryPosition += cLen;
                                break;
                            case 'I':
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);

                                preRefSeq = getSubsequence(refGenome, databaseChr, databasePosition - range, databasePosition - 1);
                                postRefSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + range - 1);

                                preQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition - range, queryPosition - 1);
                                postQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition + cLen, queryPosition + cLen + range - 1);
                                ovcffile << databaseChr << "\t" << std::to_string(databasePosition) << "\t" << databaseChr << "_" << std::to_string(databasePosition) << "\t" << "-" << "\t" << querySeq << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" /*<< "\t" <<
                                            preRefSeq << "\t" << postRefSeq << "\t" << preQuerySeq << "\t" << postQuerySeq*/ << std::endl;
                                queryPosition += cLen;
                                break;
                            case 'D':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                preRefSeq = getSubsequence(refGenome, databaseChr, databasePosition - range, databasePosition - 1);
                                postRefSeq = getSubsequence(refGenome, databaseChr, databasePosition + cLen, databasePosition + cLen + range - 1);

                                preQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition - range, queryPosition - 1);
                                postQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + range - 1);

                                ovcffile << databaseChr << "\t" << std::to_string(databasePosition) << "\t" << databaseChr << "_" << std::to_string(databasePosition) << "\t" << refSeq << "\t" << "-" << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" /*<< "\t" <<
                                           preRefSeq << "\t" << postRefSeq << "\t" << preQuerySeq << "\t" << postQuerySeq */<< std::endl;
                                databasePosition += cLen;
                                break;
                            case 'N': // intron looks like deletion
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                preRefSeq = getSubsequence(refGenome, databaseChr, databasePosition - range, databasePosition - 1);
                                postRefSeq = getSubsequence(refGenome, databaseChr, databasePosition + cLen, databasePosition + cLen + range - 1);

                                preQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition - range, queryPosition - 1);
                                postQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + range - 1);

                                ovcffile << databaseChr << "\t" << std::to_string(databasePosition) << "\t" << databaseChr << "_" << std::to_string(databasePosition) << "\t" << refSeq << "\t" << "-" << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" /*<< "\t" <<
                                            preRefSeq << "\t" << postRefSeq << "\t" << preQuerySeq << "\t" << postQuerySeq */<< std::endl;
                                databasePosition += cLen;
                                break;
                            case 'H':
                                break;
                            case 'S':
                                break;
                            case 'P':
                                std::cout << "current we could not deal with cigar letter P" << std::endl;
                                break;
                            default:
                                std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                                break;
                        }
                    }
                } else {
                    int32_t databasePosition = databaseStart;
                    int32_t queryPosition = queryEnd;
                    for (int i = 0; i < cigarElems.size(); ++i) {
                        std::string cVal = cigarElems[i];
                        char cLetter = cVal[cVal.length() - 1];
                        int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                        if (i == cigarElems.size() - 1 && (cLetter == 'H' || cLetter == 'S')) { // ignore the last soft/hard clipping
                            tailClipping += cLen;
                            continue;
                        }
                        std::string refSeq;
                        std::string querySeq;

                        std::string preQuerySeq;
                        std::string postQuerySeq;
                        std::string preRefSeq;
                        std::string postRefSeq;

                        switch (cLetter) {
                            case 'M':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);
                                for (int32_t ii = 0; ii < refSeq.size(); ++ii) {
                                    if (refSeq[ii] != querySeq[ii]) {
                                        ovcffile << databaseChr << "\t" << std::to_string(databasePosition + ii) << "\t" << databaseChr << "_" << std::to_string(databasePosition + ii) << "\t" << refSeq[ii] << "\t" << querySeq[ii] << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1"
                                                 << std::endl;
                                    }
                                }
                                databasePosition += cLen;
                                queryPosition -= cLen;
                                break;
                            case '=':
                                databasePosition += cLen;
                                queryPosition -= cLen;
                                break;
                            case 'X':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);
                                for (int32_t ii = 0; ii < refSeq.size(); ++ii) {
                                    if (refSeq[ii] != querySeq[ii]) {
                                        ovcffile << databaseChr << "\t" << std::to_string(databasePosition + ii) << "\t" << databaseChr << "_" << std::to_string(databasePosition + ii) << "\t" << refSeq[ii] << "\t" << querySeq[ii] << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1"
                                                 << std::endl;
                                    }
                                }
                                databasePosition += cLen;
                                queryPosition -= cLen;
                                break;
                            case 'I':
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);

                                preRefSeq = getSubsequence(refGenome, databaseChr, databasePosition - range, databasePosition - 1);
                                postRefSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + range - 1);

                                preQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition + 1, queryPosition + range, NEGATIVE);
                                postQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition - cLen - range + 1, queryPosition - cLen, NEGATIVE);

                                ovcffile << databaseChr << "\t" << std::to_string(databasePosition) << "\t" << databaseChr << "_" << std::to_string(databasePosition) << "\t" << "-" << "\t" << querySeq << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" << "\t" <<
                                         preRefSeq << "\t" << postRefSeq << "\t" << preQuerySeq << "\t" << postQuerySeq << std::endl;
                                queryPosition -= cLen;
                                break;
                            case 'D':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);

                                preRefSeq = getSubsequence(refGenome, databaseChr, databasePosition - range, databasePosition - 1);
                                postRefSeq = getSubsequence(refGenome, databaseChr, databasePosition + cLen, databasePosition + cLen + range - 1);

                                preQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition + 1, queryPosition + range - 1, NEGATIVE);
                                postQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition - range + 1, queryPosition, NEGATIVE);

                                ovcffile << databaseChr << "\t" << std::to_string(databasePosition) << "\t" << databaseChr << "_" << std::to_string(databasePosition) << "\t" << refSeq << "\t" << "-" << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" << "\t" <<
                                         preRefSeq << "\t" << postRefSeq << "\t" << preQuerySeq << "\t" << postQuerySeq << std::endl;
                                databasePosition += cLen;
                                break;
                            case 'N':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                preRefSeq = getSubsequence(refGenome, databaseChr, databasePosition - range, databasePosition - 1);
                                postRefSeq = getSubsequence(refGenome, databaseChr, databasePosition + cLen, databasePosition + cLen + range - 1);

                                preQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition + 1, queryPosition + range - 1, NEGATIVE);
                                postQuerySeq = getSubsequence(queryGenome, queryChr, queryPosition - range + 1, queryPosition, NEGATIVE);
                                ovcffile << databaseChr << "\t" << std::to_string(databasePosition) << "\t" << databaseChr << "_" << std::to_string(databasePosition) << "\t" << refSeq << "\t" << "-" << "\t50\tPASS\tDP=1\tGT:GQ:DP\t1|1:50:1" << "\t" <<
                                         preRefSeq << "\t" << postRefSeq << "\t" << preQuerySeq << "\t" << postQuerySeq << std::endl;
                                databasePosition += cLen;
                                break;
                            case 'H':
                                break;
                            case 'S':
                                break;
                            case 'P':
                                std::cout << "current we could not deal with cigar letter P" << std::endl;
                                break;
                            default:
                                std::cout << "current we could not deal with cigar letter " << cLetter << std::endl;
                                break;
                        }
                    }
                }
            }
        }
    }
    ovcffile.close();
}
