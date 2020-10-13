/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/

#include "controlLayer.h"

std::string softwareName = "proali";


int gff2seq(int argc, char** argv, std::map<std::string, std::string>& parameters ){
    std::stringstream usage;
    usage << "Usage: "<<softwareName<<" gff2seq -i inputGffFile -r inputGenome -o outputCdsSequences " << std::endl<<
          "Options" << std::endl <<
          " -h        produce help message" << std::endl <<
          " -i FILE   reference genome in GFF/GTF format" << std::endl <<
          " -r FILE   genome sequence in fasta format" << std::endl <<
          " -o FILE   output file of the longest CDS for each gene" << std::endl <<
          " -m INT    minimum exon to output (default: 20)" << std::endl <<std::endl;

    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-o") ){
        std::string inputGffFile = inputParser.getCmdOption("-i");
        std::string genome = inputParser.getCmdOption("-r");
        std::string outputCdsSequences = inputParser.getCmdOption("-o");

        int minExon;
        if( inputParser.cmdOptionExists("-m") ){
            minExon = std::stoi( inputParser.getCmdOption("-m") );
        }else{
            minExon=20;
        }

        getSequences( inputGffFile, genome, outputCdsSequences, parameters, minExon);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 0;
}


int genomeAlignment(int argc, char** argv, std::map<std::string, std::string>& parameters ) {
    std::stringstream usage;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    double inversion_PENALTY = -1;
    double MIN_ALIGNMENT_SCORE = 3;
    bool considerInversion = false;

    usage << "Usage: " << softwareName
          << " genoAli -i refGffFile -r refGenome  -t targetGff -s targetGenome -o output GFF/GTF file " << std::endl <<
          "Options" << std::endl <<
          " -h          produce help message" << std::endl <<
          " -i  FILE    reference GFF/GTF file" << std::endl <<
          " -r  FILE    reference genome sequence" << std::endl <<
          " -a  FILE    sam file" << std::endl <<
          " -s  FILE    target genome sequence" << std::endl <<
          " -n  FILE    output anchors" << std::endl <<
          " -o  FILE    output file in maf format" << std::endl <<
          " -f  FILE    output sequence alignment for each block in maf format (any block with size larger than window width would be ignored)" << std::endl <<
          " -v  FILE    output variant calling in vcf format (conflict with -IV)" << std::endl <<
          " -m  INT     minimum exon to use (default: 20, should be identical with the setting of gff2seq function)"
          " -A  INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B  INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -IV         whether call inversions (default: false)"  << std::endl<<
          " -IC DOUBLE  penlty to have a non-linear match in inversion region (default: " << inversion_PENALTY << ")"  << std::endl<<
          " -I  DOUBLE  minimum score to keep an inversion (default: " << MIN_ALIGNMENT_SCORE << ")"  << std::endl<<
          " -w  INT     sequence alignment window width (default: 10000)" << std::endl << std::endl;

    InputParser inputParser(argc, argv);
    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help")) {
        std::cerr << usage.str();
    } else if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
               inputParser.cmdOptionExists("-a") && inputParser.cmdOptionExists("-s") &&
             (inputParser.cmdOptionExists("-n") || inputParser.cmdOptionExists("-f") || inputParser.cmdOptionExists("-v")  || inputParser.cmdOptionExists("-o")  ) ) {

        std::string refGffFilePath = inputParser.getCmdOption("-i");
        std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
        std::string samFilePath=inputParser.getCmdOption("-a");
        std::string targetGenomeSequence = inputParser.getCmdOption("-s");

        std::string outPutMafFile="";
        if (inputParser.cmdOptionExists("-o") ){
            outPutMafFile = inputParser.getCmdOption("-o");
        }

        std::string outPutFragedFile;
        if ( inputParser.cmdOptionExists("-f") ){
            outPutFragedFile = inputParser.getCmdOption("-f");
        }

        std::string outPutVcfFile = "";
        if (inputParser.cmdOptionExists("-v") ){
            outPutVcfFile = inputParser.getCmdOption("-v");
        }

        if( inputParser.cmdOptionExists("-A") ){
            matchingScore = std::stoi( inputParser.getCmdOption("-A") );
        }
        if( inputParser.cmdOptionExists("-B") ){
            mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
        }
        if( inputParser.cmdOptionExists("-O1") ){
            openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
        }
        if( inputParser.cmdOptionExists("-E1") ){
            extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
        }

        size_t widownWidth = 10000;
        if( inputParser.cmdOptionExists("-w") ){
            widownWidth = std::stoi( inputParser.getCmdOption("-w") );
        }

        if( inputParser.cmdOptionExists("-I")){
            MIN_ALIGNMENT_SCORE = std::stod(inputParser.getCmdOption("-I"));
        }
        if( MIN_ALIGNMENT_SCORE <= 1 ){
            std::cout << "minimum score to keep an inversion is not larger than 1, maybe output weired inversions. please change it" << std::endl;
            return 1;
        }
        if( inputParser.cmdOptionExists("-IC")){
            inversion_PENALTY = std::stod(inputParser.getCmdOption("-IC"));
        }
        if( inputParser.cmdOptionExists("-IV")){
            considerInversion = true;
        }

        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMap;
        if (considerInversion){
            if (inputParser.cmdOptionExists("-v") ){
                std::cout << "-v and -IV are conflict with each other" << std::endl;
                return 1;
            }
        }
        int minExon;
        if( inputParser.cmdOptionExists("-m") ){
            minExon = std::stoi( inputParser.getCmdOption("-m") );
        }else{
            minExon=20;
        }

        setupAnchorsWithSpliceAlignmentResult( refGffFilePath, samFilePath, alignmentMatchsMap, inversion_PENALTY,  MIN_ALIGNMENT_SCORE, considerInversion, minExon);

        if( inputParser.cmdOptionExists("-n") ) {
            std::ofstream ofile;
            ofile.open(inputParser.getCmdOption("-n"));
            int blockIndex = 0;
            ofile << "refChr" << "\t"
                  << "referenceStart" << "\t"
                  << "referenceEnd" << "\t"
                  << "queryChr" << "\t"
                  << "queryStart" << "\t"
                  << "queryEnd" << "\t"
                  << "strand" << "\t"
                  << "gene" << "\t"
                  << "blockIndex" << std::endl;

            for (std::map<std::string, std::vector<AlignmentMatch>>::iterator it = alignmentMatchsMap.begin(); it != alignmentMatchsMap.end(); ++it ) {
                ofile << "#block begin" << std::endl;
                blockIndex++;
                for (int rangeIndex = 0; rangeIndex <  it->second.size(); ++rangeIndex) {
                    if (rangeIndex > 0) {
                        if (it->second[rangeIndex].getStrand() == POSITIVE && it->second[rangeIndex-1].getStrand()== POSITIVE ) {
                            ofile << it->second[rangeIndex].getRefChr() << "\t"
                                   << it->second[rangeIndex - 1].getRefEndPos() + 1 << "\t"
                                   << it->second[rangeIndex].getRefStartPos() - 1 << "\t"
                                   << it->second[rangeIndex].getQueryChr() << "\t"
                                   << it->second[rangeIndex - 1].getQueryEndPos() + 1 << "\t"
                                   << it->second[rangeIndex].getQueryStartPos() - 1 << "\t"
                                    << "+" << "\t"
                                    << "intergenetic" << "\t" <<
                                                                        blockIndex << std::endl;
                        } else if (it->second[rangeIndex].getStrand() == NEGATIVE && it->second[rangeIndex-1].getStrand()== NEGATIVE ) {
                            ofile << it->second[rangeIndex].getRefChr() << "\t"
                                   << it->second[rangeIndex - 1].getRefEndPos() + 1 << "\t"
                                   << it->second[rangeIndex].getRefStartPos() - 1 << "\t"
                                   << it->second[rangeIndex].getQueryChr() << "\t"
                                   << it->second[rangeIndex].getQueryEndPos() + 1 << "\t"
                                   << it->second[rangeIndex - 1].getQueryStartPos() - 1 << "\t"
                                    << "-" << "\t"
                                    << "intergenetic" << "\t" <<
                                    blockIndex << std::endl;
                        }
                    }
                    std::string thisStrand = "+";
                    if( it->second[rangeIndex].getStrand() == NEGATIVE ){
                        thisStrand = "-";
                    }
                    ofile << it->second[rangeIndex].getRefChr() << "\t"
                          << it->second[rangeIndex].getRefStartPos() << "\t"
                          << it->second[rangeIndex].getRefEndPos() << "\t"
                          << it->second[rangeIndex].getQueryChr() << "\t"
                          << it->second[rangeIndex].getQueryStartPos() << "\t"
                          << it->second[rangeIndex].getQueryEndPos() << "\t"
                          << thisStrand << "\t"
                          << it->second[rangeIndex].getReferenceGeneName() << "\t" <<
                          blockIndex << std::endl;
                }
                ofile << "#block end" << std::endl;
            }
            ofile.close();
        }

        genomeAlignmentAndVariantCalling(alignmentMatchsMap, referenceGenomeSequence, targetGenomeSequence, widownWidth,
                                         outPutMafFile, outPutVcfFile, outPutFragedFile,
                                         matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1);

        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 0;
}

int proportationalAlignment(int argc, char** argv, std::map<std::string, std::string>& parameters ) {
    std::stringstream usage;
    int windownWidth = 15000;
    double calculateIndelDistance = 3;
    double GAP_OPEN_PENALTY=-0.03;
    double INDEL_SCORE=-0.01;
    double MIN_ALIGNMENT_SCORE = 4;
    int MAX_DIST_BETWEEN_MATCHES=25;  // between maize and sorghum set it as 27
    int refMaximumTimes=1;
    int queryMaximumTimes=2;
    bool outPutAlignmentForEachInterval = false;
    bool localAlignment = false;

    int32_t matchingScore = 2;
    int32_t mismatchingPenalty = -3;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t seed_window_size = 38;
    int32_t mini_cns_score = 40;
    int32_t step_size = 8;
    int32_t matrix_boundary_distance = 0;

    int32_t scoreThreshold = 54;

    int32_t w = 10;  //this is the band width for band sequence alignments
    int32_t xDrop = 20;

    usage << "Usage: " << softwareName
          << " genoAli -i proali -r refGenome -a samFile -s targetGenome -o output mafFile " << std::endl <<
          "Options" << std::endl <<
          " -h          produce help message" << std::endl <<
          " -i  FILE    reference GFF/GTF file" << std::endl <<
          " -r  FILE    reference genome sequence" << std::endl <<
          " -a  FILE    sam file" << std::endl <<
          " -s  FILE    target genome sequence" << std::endl <<
          " -o  STRING  output file prefix" << std::endl <<
          " -w  INT     sequence alignment window width (default: " << windownWidth << ")" << std::endl <<
          " -R  INT     reference coverage(default: " << refMaximumTimes << ")" << std::endl<<
          " -Q  INT     query coverage (default: " << queryMaximumTimes << ")"  << std::endl<<
          " -A  INT     Matching score (default: " << matchingScore << ")" << std::endl <<
          " -B  INT     Mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1 INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1 INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -m  INT     minimum exon to use (default: 20, should be identical with the setting of gff2seq function)"
          " advanced options" << std::endl<<
          " -y  INT     minimum score to report a local sequence alignment (default: "<<scoreThreshold<<")" << std::endl <<
          " -sw  INT    the windows size used to run the smith-waterman algorithm to get the alignment seed (default: "<<seed_window_size<<")" << std::endl <<
          " -c  INT     minimum seeds score to trigger a local alignment extension (default: " << mini_cns_score << ")" << std::endl <<
          " -st  INT    step size for sliding the smith-waterman seeds alignment window (default: " << step_size << ")" << std::endl <<
          " -x  INT     x-drop for local alignment (default: " << xDrop << ")" << std::endl <<
          " -u  INT     xextend alignment band width (default: " << w << ")" << std::endl <<
          " -d  DOUBLE  calculateIndelDistance (default: " << calculateIndelDistance << ")"  << std::endl<<
          " -OC DOUBLE  chain open gap penalty (default: " << GAP_OPEN_PENALTY << ")"  << std::endl<<
          " -EC DOUBLE  chain extend gap penalty (default: " << INDEL_SCORE << ")"  << std::endl<<
          " -I  DOUBLE  minimum chain score (default: " << MIN_ALIGNMENT_SCORE << ")"  << std::endl<<
          " -D  INT     maximum gap size for chain (default: " << MAX_DIST_BETWEEN_MATCHES << ")"  << std::endl<<
          " -f          output alignment for each interval (default: false)"  << std::endl<<
          " -l          perform local alignment for each interval (default: false)"  << std::endl<<
          std::endl;

    InputParser inputParser(argc, argv);
    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help")) {
        std::cerr << usage.str();
    } else if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
               inputParser.cmdOptionExists("-a") && inputParser.cmdOptionExists("-s") &&
               inputParser.cmdOptionExists("-o")) {
        std::string refGffFilePath = inputParser.getCmdOption("-i");
        std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
        std::string targetGenomeSequence = inputParser.getCmdOption("-s");
        std::string outPutFilePath = inputParser.getCmdOption("-o");
        std::string samFilePath=inputParser.getCmdOption("-a");

        if( inputParser.cmdOptionExists("-w")){
            windownWidth = std::stoi(inputParser.getCmdOption("-w"));
        }
        if( inputParser.cmdOptionExists("-d")){
            calculateIndelDistance = std::stod(inputParser.getCmdOption("-d"));
        }
        if( inputParser.cmdOptionExists("-OC")){
            GAP_OPEN_PENALTY = std::stod(inputParser.getCmdOption("-OC"));
        }
        if( inputParser.cmdOptionExists("-EC")){
            INDEL_SCORE = std::stod(inputParser.getCmdOption("-EC"));
        }
        if( inputParser.cmdOptionExists("-I")){
            MIN_ALIGNMENT_SCORE = std::stod(inputParser.getCmdOption("-I"));
        }
        if( inputParser.cmdOptionExists("-D")){
            MAX_DIST_BETWEEN_MATCHES = std::stoi(inputParser.getCmdOption("-D"));
        }
        if( inputParser.cmdOptionExists("-R")){
            refMaximumTimes = std::stoi(inputParser.getCmdOption("-R"));
        }
        if( inputParser.cmdOptionExists("-Q")){
            queryMaximumTimes = std::stoi(inputParser.getCmdOption("-Q"));
        }
        if( inputParser.cmdOptionExists("-f")){
            outPutAlignmentForEachInterval = true;
        }
        if( inputParser.cmdOptionExists("-l")){
            localAlignment = true;
        }

        if( inputParser.cmdOptionExists("-A") ){
            matchingScore = std::stoi( inputParser.getCmdOption("-A") );
        }
        if( inputParser.cmdOptionExists("-B") ){
            mismatchingPenalty = std::stoi( inputParser.getCmdOption("-B") );
        }
        if( inputParser.cmdOptionExists("-O1") ){
            openGapPenalty1 = std::stoi( inputParser.getCmdOption("-O1") );
        }
        if( inputParser.cmdOptionExists("-E1") ){
            extendGapPenalty1 = std::stoi( inputParser.getCmdOption("-E1") );
        }

        if( inputParser.cmdOptionExists("-sw") ){
            seed_window_size = std::stoi( inputParser.getCmdOption("-sw") );
        }
        if( inputParser.cmdOptionExists("-c") ){
            mini_cns_score = std::stoi( inputParser.getCmdOption("-c") );
        }
        if( inputParser.cmdOptionExists("-st") ){
            step_size = std::stoi( inputParser.getCmdOption("-st") );
        }
        if( inputParser.cmdOptionExists("-y") ){
            scoreThreshold = std::stoi( inputParser.getCmdOption("-y") );
        }
        if( inputParser.cmdOptionExists("-u") ){
            w = std::stoi( inputParser.getCmdOption("-u") );
        }
        if( inputParser.cmdOptionExists("-x") ){
            xDrop = std::stoi( inputParser.getCmdOption("-x") );
        }

        if( localAlignment && outPutAlignmentForEachInterval ){
            std::cout << "please do not perform local alignment and global alignment for each interval at the same time" << std::endl;
            return 1;
        }
        int minExon;
        if( inputParser.cmdOptionExists("-m") ){
            minExon = std::stoi( inputParser.getCmdOption("-m") );
        }else{
            minExon=20;
        }

        std::vector<std::vector<AlignmentMatch>> alignmentMatchsMap;

        setupAnchorsWithSpliceAlignmentResultQuota( refGffFilePath, samFilePath, alignmentMatchsMap, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE,
                                                    MAX_DIST_BETWEEN_MATCHES, refMaximumTimes, queryMaximumTimes,
                                                    calculateIndelDistance, minExon);

        std::ofstream ofile;
        ofile.open(outPutFilePath + ".forPlotQuota");

        ofile << "refChr" << "\t"
              << "referenceStart" << "\t"
              << "referenceEnd" << "\t"
              << "queryChr" << "\t"
              << "queryStart" << "\t"
              << "queryEnd" << "\t"
              << "strand" << "\t"
              << "gene" << "\t"
              << "blockIndex" << std::endl;

        size_t  totalAnchors = 0;
        int blockIndex = 0;
        for ( std::vector<AlignmentMatch> alignmentMatchs : alignmentMatchsMap ){
            ofile << "#block begin" << std::endl;
            totalAnchors += alignmentMatchs.size();
            blockIndex++;
            for( int rangeIndex = 0; rangeIndex <alignmentMatchs.size();  ++rangeIndex){

                std::string thisStrand = "+";
                if( alignmentMatchs[rangeIndex].getStrand() == NEGATIVE ){
                    thisStrand = "-";
                }
                if (  rangeIndex >0 ){
                    if( alignmentMatchs[rangeIndex].getStrand() == POSITIVE &&  alignmentMatchs[rangeIndex-1].getStrand() == POSITIVE ){
                        ofile << alignmentMatchs[rangeIndex].getRefChr() << "\t"
                               << alignmentMatchs[rangeIndex-1].getRefEndPos()+1 << "\t"
                               << alignmentMatchs[rangeIndex].getRefStartPos() - 1 << "\t"
                               << alignmentMatchs[rangeIndex].getQueryChr() << "\t"
                               << alignmentMatchs[rangeIndex-1].getQueryEndPos()+1 << "\t"
                               << alignmentMatchs[rangeIndex].getQueryStartPos()-1 << "\t"
                                << "+" << "\t" <<
                                "intergenetic" << "\t" <<
                                blockIndex << std::endl;
                    }else if( alignmentMatchs[rangeIndex].getStrand() == NEGATIVE &&  alignmentMatchs[rangeIndex-1].getStrand() == NEGATIVE ) {
                        ofile << alignmentMatchs[rangeIndex].getRefChr() << "\t"
                               << alignmentMatchs[rangeIndex-1].getRefEndPos()+1 << "\t"
                               << alignmentMatchs[rangeIndex].getRefStartPos() - 1 << "\t"
                               << alignmentMatchs[rangeIndex].getQueryChr() << "\t"
                               << alignmentMatchs[rangeIndex].getQueryEndPos() +1 << "\t"
                               << alignmentMatchs[rangeIndex-1].getQueryStartPos()-1 << "\t"
                                << "-" << "\t" <<
                                "intergenetic" << "\t" <<
                                blockIndex << std::endl;
                    }
                }
                ofile << alignmentMatchs[rangeIndex].getRefChr() << "\t"
                      << alignmentMatchs[rangeIndex].getRefStartPos() << "\t"
                      << alignmentMatchs[rangeIndex].getRefEndPos() << "\t"
                      << alignmentMatchs[rangeIndex].getQueryChr() << "\t"
                      << alignmentMatchs[rangeIndex].getQueryStartPos() << "\t"
                      << alignmentMatchs[rangeIndex].getQueryEndPos() << "\t"
                      << thisStrand << "\t"
                      << alignmentMatchs[rangeIndex].getReferenceGeneName() << "\t" <<
                      blockIndex << std::endl;
            }
            ofile << "#block end" << std::endl;
        }
        ofile.close();

        genomeAlignment(alignmentMatchsMap, referenceGenomeSequence,  targetGenomeSequence, windownWidth, outPutFilePath, outPutAlignmentForEachInterval,localAlignment,
                        matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, seed_window_size, mini_cns_score, step_size, matrix_boundary_distance, scoreThreshold, w, xDrop);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 0;
}
