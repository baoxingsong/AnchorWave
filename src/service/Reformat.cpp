//
// Created by bs674 on 12/19/20.
//

#include "Reformat.h"



void mafTovcf( std::string & mafFile,  std::string & refGenomeFilePath, std::string & outputovcffile, const bool & gvcf ){
    std::ofstream ovcffile;
    ovcffile.open(outputovcffile);

    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string filedate = std::to_string((1900 + ltm->tm_year)) + std::to_string((1 + ltm->tm_mon));
    if( ltm->tm_mday < 10 ){
        filedate = filedate + "0" + std::to_string((ltm->tm_mday));
    }else{
        filedate = filedate + std::to_string((ltm->tm_mday));
    }

    std::string refFileName;
    std::vector<std::string> elems;
    char delim = '/';
    split(refGenomeFilePath, delim, elems);
    refFileName = elems.back();

    if( gvcf ){
        ovcffile << "##fileformat=VCFv4.2" << std::endl;
        ovcffile << "##fileDate=" << filedate << std::endl;
        ovcffile << "##source=proali" << std::endl;
        ovcffile <<"##reference=" << refFileName << std::endl;
        ovcffile << "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
        ovcffile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">" << std::endl;
        ovcffile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
        ovcffile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
        ovcffile << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" << std::endl;
        ovcffile << "##INFO=<ID=AF,Number=3,Type=Integer,Description=\"Allele Frequency\">" << std::endl;
        ovcffile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
        ovcffile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">" << std::endl;
        ovcffile << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
    }else{
        ovcffile << "##fileformat=VCFv4.3" << std::endl;
        ovcffile << "##fileDate=" << filedate << std::endl;
        ovcffile << "##source=AnchorWave" << std::endl;
        ovcffile <<"##reference=" << refFileName << std::endl;
        ovcffile <<"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
        ovcffile <<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
        ovcffile <<"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
        ovcffile <<"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
        ovcffile <<"##FILTER=<ID=q30,Description=\"Quality below 30\">" << std::endl;
    }
    std::string accession = "query";
    accession = songStrReplaceAll(accession, ".fasta", "");
    accession = songStrReplaceAll(accession, ".fa", "");
    accession.erase(std::remove(accession.begin(), accession.end(), ' '), accession.end());

    ovcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << accession << std::endl;

    std::map <std::string, std::string> refSequences;
    readFastaFileWorkWithIUPACcode(refGenomeFilePath, refSequences);

    std::ifstream infile(mafFile);
    if( ! infile.good()){
        std::cerr << "error in opening variants file " << mafFile << std::endl;
        exit (1);
    }
    std::string line;
    while (std::getline(infile, line)){
        if( line.compare(0, 1, "#")==0 ){
            continue;
        }
        songStrReplaceAll(line, "\\s+", "\t");
        std::vector<std::string> splits;
        char sep = '\t';
        split(line, sep, splits);
        if( splits.size() ==7 ){
            std::string refchr = splits[1];
            int refstart = std::stoi(splits[2]);
            int reflength = std::stoi(splits[3]);
            std::string refstrand = splits[4];
            int refchrLength = std::stoi(splits[5]);
            std::string refali = splits[6];

            std::string line2;
            std::getline(infile, line2);
            songStrReplaceAll(line2, "\\s+", "\t");
            std::vector<std::string> splits2;
            split(line2, sep, splits2);

            std::string querychr = splits2[1];
            int querystart = std::stoi(splits2[2]);
            int querylength = std::stoi(splits2[3]);
            std::string querystrand = splits2[4];
            int queryChrLength = std::stoi(splits2[5]);
            std::string queryali = splits2[6];
            songStrReplaceAll(refchr, "Zea_mays.AGPv4.dna.toplevel.fa.", "");
            songStrReplaceAll(refchr, "col.", "");
            songStrReplaceAll(refchr, "hg38.fa.", "");
            songStrReplaceAll(refchr, "B73.ref.fa.", "");
            refchr = songStrReplaceAll(refchr, "\\s", "");
            alignmentToVcf(queryali, refali, ovcffile, refchr, refSequences, refstart, gvcf);
        }
    }
    ovcffile.close();
}










void samToMaf (const std::string& samFilePath, const std::string& refGenomeFile,  const std::string& queryGenomeFile, const std::string& output){
    std::map<std::string, std::string> refGenome;
    readFastaFile(refGenomeFile, refGenome);
    std::map<std::string, std::string> queryGenome;
    readFastaFile(queryGenomeFile, queryGenome);
    std::string line;
    char delim = '\t';
    std::vector<std::string> elems;

    std::string databaseChr;
    int32_t databaseStart;
    int32_t databaseEnd;
    std::string queryChr;
    int32_t queryStart;
    int32_t queryEnd;

    std::map<std::string, std::string> lastChr;
    std::map<std::string, int32_t > lastPosition;


    std::ofstream maffile;
    maffile.open(output);
    maffile << "##maf version=1" << std::endl;

    std::map<std::string, std::map<std::string, std::vector<double>>> geneScores; //fist key is gene name, second key is chr value is a vector of similarity
    std::ifstream infile(samFilePath);
    while (std::getline(infile, line)){ // no matter the transcript in on forward strand or reverse strand, it should do not matter
        if( line[0] != '@' ){ //ignore the header
            elems.clear();
            split(line, delim, elems);
            queryStart= 1;
//            if(databaseStart>0 && transcriptHashMap.find(elems[0])!=transcriptHashMap.end() ){ // ignore those none mapping records
            queryChr=elems[0];
            if( queryChr.compare("*") != 0){ // ignore those none mapping records
                //std::cout << "begain to analysis " << line << std::endl;
                databaseChr=elems[2];

                databaseStart = std::stoi(elems[3]);
                databaseEnd = databaseStart;
                queryEnd=queryStart;
                std::vector<std::string> cigarElems;
                splitCIGAR(elems[5], cigarElems);
                int headClipping = 0;
                int tailClipping = 0;
                for(int i=0; i<cigarElems.size(); ++i) {
                    std::string cVal = cigarElems[i];
                    char cLetter = cVal[cVal.length() - 1];
                    int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                    if( i == cigarElems.size()-1 && (cLetter == 'H' || cLetter == 'S' ) ){ // ignore the last soft/hard clipping
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

                std::stringstream refAlign;
                std::stringstream queryAlign;

                if ( 0 == samFlag % 32){
                    int32_t databasePosition = databaseStart;
                    int32_t queryPosition = queryStart;
                    for(int i=0; i<cigarElems.size(); ++i) {
                        std::string cVal = cigarElems[i];
                        char cLetter = cVal[cVal.length() - 1];
                        int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                        if( i == cigarElems.size()-1 && (cLetter == 'H' || cLetter == 'S' ) ){ // ignore the last soft/hard clipping
                            tailClipping += cLen;
                            continue;
                        }
                        std::string refSeq;
                        std::string querySeq;
                        switch (cLetter) {
                            case 'M':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);
                                databasePosition += cLen;
                                queryPosition += cLen;
                                refAlign << refSeq;
                                queryAlign << querySeq;
                                break;
                            case '=':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);
                                databasePosition += cLen;
                                queryPosition += cLen;
                                refAlign << refSeq;
                                queryAlign << querySeq;
                                break;
                            case 'X':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);
                                databasePosition += cLen;
                                queryPosition += cLen;
                                refAlign << refSeq;
                                queryAlign << querySeq;
                                break;
                            case 'I':
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition + cLen - 1);
                                queryPosition += cLen;
                                queryAlign << querySeq;
                                for (int repeatI = 0; repeatI < querySeq.length(); ++repeatI) {
                                    refAlign << "-";
                                }
                                break;
                            case 'D':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                databasePosition += cLen;
                                refAlign << refSeq;
                                for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                                    queryAlign << "-";
                                }

                                break;
                            case 'N':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                databasePosition += cLen;
                                refAlign << refSeq;
                                for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                                    queryAlign << "-";
                                }
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
                    std::string temp1 = refAlign.str();
                    std::string temp2 = queryAlign.str();
                    temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
                    temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());
                    int chrWidth = 15;
                    databaseStart = databaseStart - 1; // maf file is 0-based
                    queryStart = queryStart -1; // maf file is 0-based
                    maffile << "a\tscore=" << 0 << std::endl;
                    maffile << "s\t" << std::left << std::setw(chrWidth) <<   "col." + databaseChr << "\t" << std::right
                            << std::setw(9) << databaseStart << "\t" << std::setw(9)
                            << temp1.size() << "\t+\t" << refGenome[databaseChr].size() << "\t"
                            << refAlign.str() << std::endl;

                    maffile << "s\t" << std::left << std::setw(chrWidth) << "query." + queryChr << "\t" << std::right
                            << std::setw(9) << queryStart << "\t" << std::setw(9)
                            << temp2.size() << "\t+" << "\t"<< queryGenome[queryChr].size() << "\t"
                            << queryAlign.str() << std::endl;
                    maffile << std::endl;
                }else{
                    int32_t databasePosition = databaseStart;
                    int32_t queryPosition = queryEnd;
                    for(int i=0; i<cigarElems.size(); ++i) {
                        std::string cVal = cigarElems[i];
                        char cLetter = cVal[cVal.length() - 1];
                        int cLen = stoi(cVal.substr(0, cVal.length() - 1));
                        if( i == cigarElems.size()-1 && (cLetter == 'H' || cLetter == 'S' ) ){ // ignore the last soft/hard clipping
                            tailClipping += cLen;
                            continue;
                        }
                        std::string refSeq;
                        std::string querySeq;
                        switch (cLetter) {
                            case 'M':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);
                                databasePosition += cLen;
                                queryPosition -= cLen;
                                refAlign << refSeq;
                                queryAlign << querySeq;
                                break;
                            case '=':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);
                                databasePosition += cLen;
                                queryPosition -= cLen;
                                refAlign << refSeq;
                                queryAlign << querySeq;
                                break;
                            case 'X':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);
                                databasePosition += cLen;
                                queryPosition -= cLen;
                                refAlign << refSeq;
                                queryAlign << querySeq;
                                break;
                            case 'I':
                                querySeq = getSubsequence(queryGenome, queryChr, queryPosition, queryPosition - cLen + 1, NEGATIVE);
                                queryPosition -= cLen;
                                queryAlign << querySeq;
                                for (int repeatI = 0; repeatI < querySeq.length(); ++repeatI) {
                                    refAlign << "-";
                                }
                                break;
                            case 'D':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                databasePosition += cLen;
                                refAlign << refSeq;
                                for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                                    queryAlign << "-";
                                }

                                break;
                            case 'N':
                                refSeq = getSubsequence(refGenome, databaseChr, databasePosition, databasePosition + cLen - 1);
                                databasePosition += cLen;
                                refAlign << refSeq;
                                for (int repeatI = 0; repeatI < refSeq.length(); ++repeatI) {
                                    queryAlign << "-";
                                }
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
                    std::string temp1 = refAlign.str();
                    std::string temp2 = queryAlign.str();
                    temp1.erase(std::remove(temp1.begin(), temp1.end(), '-'), temp1.end());
                    temp2.erase(std::remove(temp2.begin(), temp2.end(), '-'), temp2.end());
                    int chrWidth = 15;

                    databaseStart = databaseStart - 1; // maf file is 0-based
                    queryStart = queryStart -1; // maf file is 0-based

                    maffile << "a\tscore=" << 0 << std::endl;
                    maffile << "s\t" << std::left << std::setw(chrWidth) <<   "col." + databaseChr << "\t" << std::right
                            << std::setw(9) << databaseStart << "\t" << std::setw(9)
                            << temp1.size() << "\t+\t" << refGenome[databaseChr].size() << "\t"
                            << refAlign.str() << std::endl;

                    maffile << "s\t" << std::left << std::setw(chrWidth) << "query." + queryChr << "\t" << std::right
                            << std::setw(9) << queryStart << "\t" << std::setw(9)
                            << temp2.size() << "\t-" << "\t"<< queryGenome[queryChr].size() << "\t"
                            << queryAlign.str() << std::endl;
                    maffile << std::endl;
                }
            }
        }
    }
    maffile.close();
}



void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& vcfFix, std::map<std::string, std::string>& referenceGenome){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening variants file " << filePath << std::endl;
        exit (1);
    }
//    std::regex reg("^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");
    std::string line;
//    std::cout << "143" << std::endl;
    while (std::getline(infile, line)){
        if( line.compare(0, 1, "#")==0 ){
            //std::cout << line << std::endl;
            continue;
        }

        std::vector<std::string> splits;
        char sep = '\t';
        split(line, sep, splits);
        if( splits.size() >=5 ){
            std::string thisChromosome;
            if( vcfFix.length()>0 ){
                thisChromosome = vcfFix + splits[0];
            }else{
                thisChromosome = splits[0];
            }
            int position = std::stoi(splits[1]);
            std::string reference = splits[3];
            std::string alternative = splits[4];
            transform(reference.begin(), reference.end(), reference.begin(),::toupper);
            transform(alternative.begin(), alternative.end(), alternative.begin(),::toupper);

            if (reference.compare("-") != 0) {
                if( referenceGenome.find(thisChromosome) == referenceGenome.end() ){
                } else if ( referenceGenome[thisChromosome].size() < (position+reference.size()-1) ){
                    std::cerr << "the variant position is not in the reference genome range: " << thisChromosome << " " << (position+reference.size()-1) << ". The variant record will be ignored" << std::endl;
                }else if (referenceGenome[thisChromosome].substr(position - 1, reference.size()).compare(reference) !=0) {
                    std::cerr << "the record does not confirm with reference genome sequence at: " << thisChromosome << "\t"
                              << position << std::endl<<
                              "reference genome: " << referenceGenome[thisChromosome].substr(position - 1, reference.size()) << std::endl <<
                              "sdi: " << reference << ". The variant record will be ignored" << std::endl;
                }else{
                    Variant variant(thisChromosome, position, reference, alternative);
                    if( variantsMap.find(thisChromosome) == variantsMap.end() ){
                        variantsMap[thisChromosome]=std::vector<Variant>();
                    }
                    variantsMap[thisChromosome].push_back(variant);
                }
            }else{
                Variant variant(thisChromosome, position, reference, alternative);
                if( variantsMap.find(thisChromosome) == variantsMap.end() ){
                    variantsMap[thisChromosome]=std::vector<Variant>();
                }
                variantsMap[thisChromosome].push_back(variant);
            }
        }
    }
    infile.close();

    for(std::map<std::string, std::vector<Variant> >::iterator it=variantsMap.begin(); it!=variantsMap.end(); ++it) {
        std::sort(it->second.begin(), it->second.end());
    }

    for(std::map<std::string, std::vector<Variant> >::iterator it=variantsMap.begin(); it!=variantsMap.end(); ++it){
        int lastTotalChanged=0;
        for( std::vector<Variant>::iterator it2=(*it).second.begin(); it2!=(*it).second.end(); it2++ ){
            (*it2).setLastTotalChanged(lastTotalChanged);

            int changedPoint= (*it2).getPosition() + lastTotalChanged;
            (*it2).setChangedPoint(changedPoint);
            lastTotalChanged +=(*it2).getChanginglength();
        }
    }
}



void sdiToMaf(std::string & fastaFilePath, std::string & targetFastaFilePath, std::string & sdiFile, std::string & outfile ){

    std::map<std::string, std::string> referenceSequences;
    readFastaFileWorkWithIUPACcode( fastaFilePath, referenceSequences);

    std::map<std::string, std::string> targetSequences;
    readFastaFileWorkWithIUPACcode( targetFastaFilePath, targetSequences);

    std::string  vcfFix = "";
    std::map<std::string, std::vector<Variant> > variantsMap;
    readSdiFile(sdiFile, variantsMap, vcfFix, referenceSequences);

    std::ofstream omaffile;
    omaffile.open(outfile);
    omaffile << "##maf version=1" << std::endl;

    for( std::map<std::string, std::string>::iterator iterFasta=referenceSequences.begin(); iterFasta!=referenceSequences.end(); iterFasta++ ) {
        std::string chromosome = iterFasta->first;
        std::string referenceSequence = iterFasta->second;
        int totalSize = referenceSequence.size();
        std::string reserveString;
        reserveString.reserve(totalSize * 2.5);
        std::stringstream alternativeSequencestream(reserveString);
        std::stringstream referenceSequencestream(reserveString);
        int currentPosition = 1;
        for (std::vector<Variant>::iterator iterSdi = variantsMap[chromosome].begin(); iterSdi != variantsMap[chromosome].end(); iterSdi++) {
            int position = (*iterSdi).getPosition();
            std::string ref = (*iterSdi).getReference();
            std::string alter = (*iterSdi).getAlternative();
//            std::cout << position << "\t" << ref << "\t" << alter << std::endl;
            int changingLength = (*iterSdi).getChanginglength();
            if (currentPosition < position) {
                referenceSequencestream << referenceSequence.substr(currentPosition - 1, position - currentPosition);
                alternativeSequencestream << referenceSequence.substr(currentPosition - 1, position - currentPosition);
            } else if( currentPosition == position+1 && changingLength>0 && ref.compare("-")==0 && (*(iterSdi-1)).getChanginglength()==0 && (*(iterSdi-1)).getReference().size()==1 ){
                currentPosition = position+1;
                alternativeSequencestream.seekp(-1, alternativeSequencestream.cur);
                alternativeSequencestream << alter << (*(iterSdi-1)).getAlternative();

                referenceSequencestream.seekp(-1, referenceSequencestream.cur);
                for ( int i=0; i<alter.size(); ++i ){
                    referenceSequencestream << "-";
                }
                referenceSequencestream << (*(iterSdi-1)).getReference();
                continue;
            } else if( currentPosition == position+1 && changingLength>0 && ref.compare("-")==0 && (*(iterSdi-1)).getChanginglength() < 0 && (*(iterSdi-1)).getAlternative().compare("-")==0 ){ //todo might have problem
                currentPosition = position+1;
                alternativeSequencestream.seekp(-1, alternativeSequencestream.cur);
                alternativeSequencestream << alter << (*(iterSdi-1)).getAlternative();


                referenceSequencestream.seekp(-1, referenceSequencestream.cur);
                for ( int i=0; i<alter.size(); ++i ){
                    referenceSequencestream << "-";
                }
                referenceSequencestream << (*(iterSdi-1)).getReference();
                continue;
            } else if (currentPosition > position) {
                std::cerr << "the sdi file is not well sorted, it should be sorted with coordinate " <<std::endl
                          << chromosome << ": currentPosition:"<< currentPosition << " position:" << position << std::endl;
                exit(1);
            }
            if (changingLength == 0) {
                currentPosition = position + ref.size();
                alternativeSequencestream << alter;
                referenceSequencestream << ref;
            }
            if (changingLength > 0) {
                if ( ref.compare("-")==0 ) {
                    alternativeSequencestream << alter;
                    currentPosition = position;
                    for ( int i=0; i<alter.size(); ++i ){
                        referenceSequencestream << "-";
                    }
                } else {
                    alternativeSequencestream << alter;
                    currentPosition = position + ref.size();
                    referenceSequencestream << ref;
                    if( alter.size() > ref.size() ){
                        for ( int i=0; i<alter.size() - ref.size(); ++i ){
                            referenceSequencestream << "-";
                        }
                    }else{
                        for ( int i=0; i<ref.size() - alter.size(); ++i ){
                            alternativeSequencestream << "-";
                        }
                    }
                }
            }
            if (changingLength < 0) {
                if (alter.compare("-")==0 ) {
                    currentPosition = position - changingLength;
                    referenceSequencestream << ref;
                    for ( int i=0; i<ref.size(); ++i ){
                        alternativeSequencestream << "-";
                    }
                } else {
                    alternativeSequencestream << alter;
                    currentPosition = position + ref.size();
                    referenceSequencestream << ref;
                    if ( alter.size() > ref.size() ){
                        for ( int i=0; i<alter.size() - ref.size(); ++i ){
                            referenceSequencestream << "-";
                        }
                    }else{
                        for ( int i=0; i<ref.size() - alter.size(); ++i ){
                            alternativeSequencestream << "-";
                        }
                    }
                }
            }
        }

        alternativeSequencestream << referenceSequence.substr(currentPosition - 1, totalSize - currentPosition + 1);
        referenceSequencestream << referenceSequence.substr(currentPosition - 1, totalSize - currentPosition + 1);

        omaffile << "a\tscore=" << 0 << std::endl;
        omaffile << "s\t" << std::left << std::setw(20) <<  "col." + chromosome << "\t" << std::right
                 << std::setw(9) << 0 << "\t" << std::setw(9)
                 << referenceSequences[chromosome].size() << "\t+\t" << referenceSequences[chromosome].size() << "\t"
                 << referenceSequencestream.str() << std::endl;
        omaffile << "s\t" << std::left << std::setw(20) << "query." + chromosome << "\t" << std::right
                 << std::setw(9) << 0 << "\t" << std::setw(9)
                 << targetSequences[chromosome].size() << "\t+\t" << targetSequences[chromosome].size() << "\t"
                 << alternativeSequencestream.str() << std::endl;
        omaffile << std::endl;
    }
}
