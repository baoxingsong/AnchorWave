/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  09/25/2020 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

#include "controlLayer.h"

int gff2seq(int argc, char **argv) {
    std::stringstream usage;
    usage << "Usage: " << PROGRAMNAME << " gff2seq -i inputGffFile -r inputGenome -o outputSequences " << std::endl <<
          "Options" << std::endl <<
          " -h        produce help message" << std::endl <<
          " -i FILE   reference genome annotation in GFF/GTF format" << std::endl <<
          " -r FILE   reference genome sequence in fasta format" << std::endl <<
          " -o FILE   output file of the longest CDS/exon for each gene" << std::endl <<
          " -x        use exon records instead of CDS from the GFF file" << std::endl <<
          " -m INT    minimum exon length to output (default: 20)" << std::endl << std::endl;

    InputParser inputParser(argc, argv);
    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help") || !inputParser.cmdOptionExists("-i") || !inputParser.cmdOptionExists("-r") || !inputParser.cmdOptionExists("-o")) {
        std::cerr << usage.str();
        return 1;
    } else if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-o") ) {

        std::string inputGffFile = inputParser.getCmdOption("-i");
        std::string genome = inputParser.getCmdOption("-r");
        std::string outputCdsSequences = inputParser.getCmdOption("-o");

        int minExon;
        if (inputParser.cmdOptionExists("-m")) {
            minExon = std::stoi(inputParser.getCmdOption("-m"));
        } else {
            minExon = 20;
        }

        bool exonModel = false;
        if (inputParser.cmdOptionExists("-x")) {
            exonModel = true;
        }
        getSequences(inputGffFile, genome, outputCdsSequences, minExon, exonModel);

        return 0;
    }else{
        std::cerr << usage.str();
        return 1;
    }
}

// generate alignmentMatchsMap for genomeAlignment from anchors file.
bool genGenoMapFromFile(std::string path_anchors, std::string wholeCommand, std::map<std::string, std::vector<AlignmentMatch>> & map_v_am) {
    std::ifstream infile_anchors(path_anchors);

    std::string line;
    bool same_command = false;
    bool has_begin = false;
    bool has_end = false;
    int count_begin = 0;
    int count_end = 0;
    int line_no = 0;

    while (std::getline(infile_anchors, line)) {
        if(line_no++ == 0) {
            if(line == std::string("#" + std::string(PROGRAMNAME) + " " + wholeCommand)) {
                same_command = true;
            }
            else {
                break;
            }
        }

        if (line[0] != '#') {
            continue;
        }

        if(line.find("#block begin")) {
            has_begin = true;
            count_begin++;
        }

        if(line.find("#block end")) {
            has_end = true;
            count_end++;
        }
    }

    // has generated anchors file before.
    if(same_command && has_begin && has_end && (count_begin == count_end)) {
        // get alignmentMatchsMap from anchors file.
        std::ifstream infile(path_anchors);
        if (!infile.good()) {
            std::cerr << "error in opening anchors file " << path_anchors << std::endl;
        }

        std::vector <AlignmentMatch> v_am;
        std::string line;
        bool flag_begin = false;
        int line_std_n_no = 0;
        int line_std_p_no = 0;
        STRAND last_strand = POSITIVE;

        while (std::getline(infile, line)) {
            // --line like :    5B	387734189	387739576	GWHBJBH00000005	293115001	293120759	+	transcript:TraesCS5B02G214600.2	1	0.99872
            // refChr	referenceStart	referenceEnd	queryChr	queryStart	queryEnd	strand	gene	score
            if (!flag_begin && line.find("#block begin") != std::string::npos) {
                flag_begin = true;
                line_std_n_no = 0;
                line_std_p_no = 0;
                continue;
            }

            if (!flag_begin) {
                continue;
            }

            if (line.find("#block end") != std::string::npos) {
                map_v_am[v_am[0].getRefChr()] = v_am;
                v_am.clear();
                flag_begin = false;
            }

            int size = line.size();

            char refChr[size];
            uint32_t referenceStart, referenceEnd;
            uint32_t last_referenceEnd;
            char queryChr[size];
            uint32_t queryStart, queryEnd;
            uint32_t last_queryStart;
            char strand_[size];
            char c_gene[size];
            char c_score[size];

            int ret = sscanf(line.c_str(), "%s%d%d%s%d%d%s%s%*s%s", refChr, &referenceStart, &referenceEnd, queryChr, &queryStart, &queryEnd, strand_, c_gene, c_score);

            if (ret == 9) {
                STRAND strand;
                if (strand_[0] == '-') {
                    strand = NEGATIVE;
                } else {
                    strand = POSITIVE;
                }

                if (last_strand != strand) {
                    if (strand == POSITIVE) {
                        line_std_p_no = 0;
                        line_std_p_no++;
                    } else {
                        line_std_n_no = 0;
                        line_std_n_no++;
                    }
                    last_strand = strand;
                } else {
                    if (strand == POSITIVE) {
                        int r = line_std_p_no++ % 2;
                        if (r != 0) {
                            last_referenceEnd = referenceEnd;
                            last_queryStart = queryStart;
                            continue;
                        }
                    } else {
                        bool b_2 = (last_referenceEnd < referenceStart) && (last_queryStart > queryEnd);
                        if (b_2) {
                            if (line_std_n_no++ % 2 == 1) {
                                last_referenceEnd = referenceEnd;
                                last_queryStart = queryStart;
                                continue;
                            }
                        }
                    }
                }

                double score = 1;
                if (std::string(c_score) != "1") {
                    score = 0;
                }

                last_referenceEnd = referenceEnd;
                last_queryStart = queryStart;

                AlignmentMatch am = AlignmentMatch(std::string(refChr), std::string(queryChr), referenceStart, referenceEnd, queryStart, queryEnd, score, strand, std::string(c_gene), "");
                v_am.push_back(am);
            }
        }

        infile.close();
        return true;
    }

    return false;
}

int genomeAlignment(int argc, char **argv) {
    std::stringstream usage;

    int32_t matchingScore = 0;
    int32_t mismatchingPenalty = -6;
    int32_t openGapPenalty1 = -8;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -75;
    int32_t extendGapPenalty2 = -1;

    int minExon = 20;
    double minimumSimilarity = 0.95;
    double minimumSimilarity2 = 0.2;

    double inversion_PENALTY = -1;
    double MIN_ALIGNMENT_SCORE = 2;
    bool considerInversion = false;

    int32_t wfaSize3 = 100000; // if the inter-anchor length is shorter than this value, stop trying to find new anchors
    int64_t windowWidth = 100000;
    int expectedCopies = 1;
    double maximumSimilarity = 0.6; // the maximum simalarity between secondary hist the primary hit. If the second hit is too similary with primary hit, that is unwanted duplications

    bool searchForNewAnchors = true;

    int threads = 1;
    bool exonModel = false;

    usage << "Usage: " << PROGRAMNAME
          << " genoAli -i refGffFile -r refGenome -a cds.sam -as cds.fa -ar ref.sam -s targetGenome -n outputAnchorFile -o output.maf -f output.fragmentation.maf " << std::endl <<
          "Options" << std::endl <<
          " -h           produce help message" << std::endl <<
          " -i   FILE    reference GFF/GTF file" << std::endl <<
          " -r   FILE    reference genome sequence file in fasta format" << std::endl <<
          " -as  FILE    anchor sequence file. (output from the gff2seq command)" << std::endl <<
          " -a   FILE    sam file generated by mapping conserved sequence to query genome" << std::endl <<
          " -s   FILE    target genome sequence file in fasta format" << std::endl <<
          "              Those sequences with the same name in the reference genome and query genome file would be aligned" << std::endl <<
          " -n   FILE    output anchors file" << std::endl <<
          " -o   FILE    output file in maf format" << std::endl <<
          " -f   FILE    output sequence alignment for each anchor/inter-anchor region in maf format" << std::endl <<
          " -t   INT     number of threads (default: " << threads << ")" << std::endl <<
          " -m   INT     minimum exon length to use (default: " << minExon << ", should be identical with the setting of gff2seq function)" << std::endl <<
          " -mi  DOUBLE  minimum full-length CDS anchor hit similarity to use (default:" << minimumSimilarity << ")" << std::endl <<
          " -mi2 DOUBLE  minimum novel anchor hit similarity to use (default:" << minimumSimilarity2 << ")" << std::endl <<
          " -ar  FILE    sam file generated by mapping conserved sequence to reference genome" << std::endl <<
          " -w   INT     sequence alignment window width (default: " << windowWidth << ")" << std::endl <<
          " -fa3 INT     if the inter-anchor length is shorter than this value, stop trying to find new anchors (default: " << wfaSize3 << ")" << std::endl <<
          " -B   INT     mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1  INT     gap open penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1  INT     gap extension penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -O2  INT     gap open penalty 2 (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2  INT     gap extension penalty 2 (default: " << extendGapPenalty2 << ")" << std::endl <<
          " -IV          whether to call inversions (default: false)" << std::endl <<
          " -IC  DOUBLE  penalty for having a non-linear match in inversion region (default: " << inversion_PENALTY << ")" << std::endl <<
          "              We use IC * alignment_similarity as the penalty in the inversion block" << std::endl <<
          " -I   DOUBLE  minimum score to keep an inversion (default: " << MIN_ALIGNMENT_SCORE << ")" << std::endl <<
          " -e   INT     maximum expected copy number of each gene on each chromosome (default: " << expectedCopies << ")" << std::endl << // this is used to duplicated anchors from the sam file
          "              This prevents using tandem duplicated genes to identify collinear block" << std::endl <<
          " -y   DOUBLE  minimal ratio of e+1 similarity to 1 similarity to drop an anchor (default: " << maximumSimilarity << ")" << std::endl <<
          " -ns          do not search for new anchors (default: false)" << std::endl <<
          " -x           use exon records instead of CDS from the GFF file (should be identical with the setting of gff2seq function)" << std::endl
          << std::endl;

    InputParser inputParser(argc, argv);

    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help")) {
        std::cerr << usage.str();
        return 1;
    }
    else if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
               inputParser.cmdOptionExists("-a") && inputParser.cmdOptionExists("-s") && inputParser.cmdOptionExists("-as") &&
               (inputParser.cmdOptionExists("-n") || inputParser.cmdOptionExists("-f") || inputParser.cmdOptionExists("-o") || inputParser.cmdOptionExists("-l"))) {

        std::string refGffFilePath = inputParser.getCmdOption("-i");
        std::string path_ref_GenomeSequence = inputParser.getCmdOption("-r");
        std::string cdsSequenceFile = inputParser.getCmdOption("-as");
        std::string samFilePath = inputParser.getCmdOption("-a");
        std::string path_target_GenomeSequence = inputParser.getCmdOption("-s");
        std::string outPutMafFile;

        if (inputParser.cmdOptionExists("-o")) {
            outPutMafFile = inputParser.getCmdOption("-o");
        }

        std::string outPutFragedFile;
        if (inputParser.cmdOptionExists("-f")) {
            outPutFragedFile = inputParser.getCmdOption("-f");
        }

        if (inputParser.cmdOptionExists("-t")) {
            threads = std::stoi(inputParser.getCmdOption("-t"));
        }

        if (inputParser.cmdOptionExists("-m")) {
            minExon = std::stoi(inputParser.getCmdOption("-m"));
        }
        if (inputParser.cmdOptionExists("-mi")) {
            minimumSimilarity = std::stod(inputParser.getCmdOption("-mi"));
        }
        if (inputParser.cmdOptionExists("-mi2")) {
            minimumSimilarity2 = std::stod(inputParser.getCmdOption("-mi2"));
        }

        std::string referenceSamFilePath;
        if (inputParser.cmdOptionExists("-ar")) {
            referenceSamFilePath = inputParser.getCmdOption("-ar");
        }

        if (inputParser.cmdOptionExists("-fa3")) {
            wfaSize3 = std::stoi(inputParser.getCmdOption("-fa3"));
        }

        if (inputParser.cmdOptionExists("-B")) {
            mismatchingPenalty = std::stoi(inputParser.getCmdOption("-B"));
            if (mismatchingPenalty >= 0) {
                std::cout << "parameter of B should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-O1")) {
            openGapPenalty1 = std::stoi(inputParser.getCmdOption("-O1"));
            if (openGapPenalty1 >= 0) {
                std::cout << "parameter of O1 should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-E1")) {
            extendGapPenalty1 = std::stoi(inputParser.getCmdOption("-E1"));
            if (extendGapPenalty1 >= 0) {
                std::cout << "parameter of E1 should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-O2")) {
            openGapPenalty2 = std::stoi(inputParser.getCmdOption("-O2"));
            if (openGapPenalty2 >= 0) {
                std::cout << "parameter of O1 should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-E2")) {
            extendGapPenalty2 = std::stoi(inputParser.getCmdOption("-E2"));
            if (extendGapPenalty2 > 0) {
                std::cout << "parameter of E1 should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-IV")) {
            considerInversion = true;
        }
        if (inputParser.cmdOptionExists("-IC")) {
            inversion_PENALTY = std::stod(inputParser.getCmdOption("-IC"));
            if (inversion_PENALTY >= 0) {
                std::cout << "parameter of IC should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-I")) {
            MIN_ALIGNMENT_SCORE = std::stod(inputParser.getCmdOption("-I"));
        }

        if (MIN_ALIGNMENT_SCORE < 1) {
            std::cout << "minimum score to keep an inversion is not larger than 1, maybe output weired inversions. please change it" << std::endl;
            return 1;
        }
        if (inputParser.cmdOptionExists("-e")) {
            expectedCopies = std::stoi(inputParser.getCmdOption("-e"));
        }

        if (inputParser.cmdOptionExists("-y")) {
            maximumSimilarity = std::stod(inputParser.getCmdOption("-y"));
        }
        if (inputParser.cmdOptionExists("-w")) {
            windowWidth = std::stoi(inputParser.getCmdOption("-w"));
        }

        if (inputParser.cmdOptionExists("-ns")) {
            searchForNewAnchors = false;
        }

        if (inputParser.cmdOptionExists("-x")) {
            exonModel = true;
        }

        std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
        readFastaFile(path_ref_GenomeSequence, map_ref);

        std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
        readFastaFile(path_target_GenomeSequence, map_qry);

        std::cout << "setupAnchorsWithSpliceAlignmentResult begin " << std::endl;
        std::map<std::string, std::vector<AlignmentMatch>> map_v_am;

        std::string path_anchors = inputParser.getCmdOption("-n");

        bool flag_CanGen = false;
        std::ifstream infile_anchors(path_anchors);
        if (infile_anchors.good()) {
            std::string wholeCommand = argv[0];
            for (int i = 1; i < argc; ++i) {
                wholeCommand = wholeCommand + " " + argv[i];
            }

            flag_CanGen = genGenoMapFromFile(path_anchors, wholeCommand, map_v_am);
            std::cout << "map_v_am generated!" << std::endl;
        }

        if(!flag_CanGen) {
            setupAnchorsWithSpliceAlignmentResult(refGffFilePath, cdsSequenceFile, samFilePath, map_v_am,
                                                  inversion_PENALTY, MIN_ALIGNMENT_SCORE, considerInversion, minExon, windowWidth, minimumSimilarity, minimumSimilarity2,
                                                  map_ref, map_qry,
                                                  expectedCopies, maximumSimilarity, referenceSamFilePath, wfaSize3, searchForNewAnchors, exonModel);

            std::cout << "setupAnchorsWithSpliceAlignmentResult done!" << std::endl;

            for (std::map < std::string, std::vector < AlignmentMatch >> ::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
                myAlignmentMatchSort(it->second, inversion_PENALTY, MIN_ALIGNMENT_SCORE, false, false);
            }

            std::cout << "myAlignmentMatchSort done!" << std::endl;

            if (inputParser.cmdOptionExists("-n")) {
                std::string wholeCommand = argv[0];
                for (int i = 1; i < argc; ++i) {
                    wholeCommand = wholeCommand + " " + argv[i];
                }

                std::ofstream ofile;
                ofile.open(inputParser.getCmdOption("-n"));
                ofile << "#" << PROGRAMNAME << " " << wholeCommand << std::endl;
                int blockIndex = 0;
                ofile << "refChr" << "\t"
                      << "referenceStart" << "\t"
                      << "referenceEnd" << "\t"
                      << "queryChr" << "\t"
                      << "queryStart" << "\t"
                      << "queryEnd" << "\t"
                      << "strand" << "\t"
                      << "gene" << "\t"
                      << "blockIndex" << "\tscore" << std::endl;

                for (std::map < std::string, std::vector < AlignmentMatch >> ::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
                    ofile << "#block begin" << std::endl;
                    std::vector < AlignmentMatch > v_am = it->second;
                    blockIndex++;
                    bool hasInversion = false;
                    for (size_t i = 0; i < v_am.size(); ++i) {
                        if (i > 0) {
                            if (v_am[i].getStrand() == POSITIVE && v_am[i - 1].getStrand() == POSITIVE) {
                                ofile << v_am[i].getRefChr() << "\t"
                                      << v_am[i - 1].getRefEndPos() + 1 << "\t"
                                      << v_am[i].getRefStartPos() - 1 << "\t"
                                      << v_am[i].getQueryChr() << "\t"
                                      << v_am[i - 1].getQueryEndPos() + 1 << "\t"
                                      << v_am[i].getQueryStartPos() - 1 << "\t"
                                      << "+" << "\t"
                                      << "interanchor" << "\t"
                                      << blockIndex << "\tNA" << std::endl;
                            } else if (v_am[i].getStrand() == NEGATIVE && v_am[i - 1].getStrand() == NEGATIVE
                                       && v_am[i - 1].getRefEndPos() < v_am[i].getRefStartPos()
                                       && v_am[i - 1].getQueryStartPos() > v_am[i].getQueryEndPos()) {

                                ofile << v_am[i].getRefChr() << "\t"
                                      << v_am[i - 1].getRefEndPos() + 1 << "\t"
                                      << v_am[i].getRefStartPos() - 1 << "\t"
                                      << v_am[i].getQueryChr() << "\t"
                                      << v_am[i].getQueryEndPos() + 1 << "\t"
                                      << v_am[i - 1].getQueryStartPos() - 1 << "\t"
                                      << "-" << "\t"
                                      << "interanchor" << "\t"
                                      << blockIndex << "\tNA" << std::endl;
                            }
                        }

                        std::string thisStrand = "+";
                        if (v_am[i].getStrand() == NEGATIVE) {
                            thisStrand = "-";
                            hasInversion = true;
                        }

                        ofile <<  v_am[i].getRefChr() << "\t"
                              << v_am[i].getRefStartPos() << "\t"
                              << v_am[i].getRefEndPos() << "\t"
                              << v_am[i].getQueryChr() << "\t"
                              << v_am[i].getQueryStartPos() << "\t"
                              << v_am[i].getQueryEndPos() << "\t"
                              << thisStrand << "\t"
                              << v_am[i].getReferenceGeneName() << "\t";

                        if (v_am[i].getReferenceGeneName().find("localAlignment") == std::string::npos) {
                            ofile << blockIndex << "\t" << v_am[i].getScore() << std::endl;
                        } else {
                            ofile << blockIndex << "\t" << "NA" << std::endl;
                        }
                    }

                    int i = v_am.size() - 1;
                    if (!hasInversion) {
                        size_t size_sr2 = getSequenceSizeFromPath2(map_ref[v_am[i].getRefChr()]);
                        size_t size_sq2 = getSequenceSizeFromPath2(map_qry[v_am[i].getQueryChr()]);

                        ofile << v_am[i].getRefChr() << "\t"
                              << v_am[i].getRefEndPos() + 1 << "\t"
                              << size_sr2 << "\t"
                              << v_am[i].getQueryChr() << "\t"
                              << v_am[i].getQueryEndPos() + 1 << "\t"
                              << size_sq2 << "\t"
                              << "+" << "\t"
                              << "interanchor" << "\t"
                              << blockIndex << "\tNA" << std::endl;
                    }

                    ofile << "#block end" << std::endl;
                }
                ofile.close();
            }

            std::cout << "anchors generate done!" << std::endl;
        }

        if (inputParser.cmdOptionExists("-f") || inputParser.cmdOptionExists("-o") || inputParser.cmdOptionExists("-l")) {
            genomeAlignmentAndVariantCalling(map_v_am, path_ref_GenomeSequence, path_target_GenomeSequence,
                                             windowWidth,
                                             outPutMafFile, outPutFragedFile,
                                             matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                                             openGapPenalty2, extendGapPenalty2,
                                             threads);

            std::cout << "AnchorWave done!" << std::endl;
        }
    }
    else {
        std::cerr << usage.str();
        return 1;
    }

    return 0;
}

// generate alignmentMatchsMap for proportionalAlignment from anchors file.
bool genProVectorFromFile(std::string path_anchors, std::string wholeCommand, std::vector<std::vector<AlignmentMatch>> & v_v_am) {
    std::ifstream infile_anchors(path_anchors);

    std::string line;
    bool same_command = false;
    bool has_begin = false;
    bool has_end = false;
    int count_begin = 0;
    int count_end = 0;
    int line_no = 0;

    while (std::getline(infile_anchors, line)) {
        if(line_no++ == 0) {
            if(line == std::string("#" + std::string(PROGRAMNAME) + " " + wholeCommand)) {
                same_command = true;
            }
            else {
                break;
            }
        }

        if (line[0] != '#') {
            continue;
        }

        if(line.find("#block begin")) {
            has_begin = true;
            count_begin++;
        }

        if(line.find("#block end")) {
            has_end = true;
            count_end++;
        }
    }

    // has generated anchors file before.
    if(same_command && has_begin && has_end && (count_begin == count_end)) {
        // get alignmentMatchsMap from anchors file.
        std::ifstream infile(path_anchors);
        if (!infile.good()) {
            std::cerr << "error in opening anchors file " << path_anchors << std::endl;
        }

        std::vector <AlignmentMatch> v_am;
        std::string line;
        bool flag_begin = false;
        int line_std_n_no = 0;
        int line_std_p_no = 0;
        STRAND last_strand = POSITIVE;

        while (std::getline(infile, line)) {
            // --line like :    5B	387734189	387739576	GWHBJBH00000005	293115001	293120759	+	transcript:TraesCS5B02G214600.2	1	0.99872
            // refChr	referenceStart	referenceEnd	queryChr	queryStart	queryEnd	strand	gene	score
            if (!flag_begin && line.find("#block begin") != std::string::npos) {
                flag_begin = true;
                line_std_n_no = 0;
                line_std_p_no = 0;
                continue;
            }

            if (!flag_begin) {
                continue;
            }

            if (line.find("#block end") != std::string::npos) {
                v_v_am.push_back(v_am);
                v_am.clear();
                flag_begin = false;
            }

            int size = line.size();

            char refChr[size];
            uint32_t referenceStart, referenceEnd;
            char queryChr[size];
            uint32_t queryStart, queryEnd;
            char strand_[size];
            char c_gene[size];
            char c_score[size];

            int ret = sscanf(line.c_str(), "%s%d%d%s%d%d%s%s%*s%s", refChr, &referenceStart, &referenceEnd, queryChr, &queryStart, &queryEnd, strand_, c_gene, c_score);

            if (ret == 9) {
                STRAND strand;
                if (strand_[0] == '-') {
                    strand = NEGATIVE;
                } else {
                    strand = POSITIVE;
                }

                if (last_strand != strand) {
                    if (strand == POSITIVE) {
                        line_std_p_no = 0;
                        line_std_p_no++;
                    } else {
                        line_std_n_no = 0;
                        line_std_n_no++;
                    }
                    last_strand = strand;
                } else {
                    if (strand == POSITIVE) {
                        int r = line_std_p_no++ % 2;
                        if (r != 0) {
                            continue;
                        }
                    } else {
                        if (line_std_n_no++ % 2 == 1) {
                            continue;
                        }
                    }
                }

                double score = 1;
                if (std::string(c_score) != "1") {
                    score = 0;
                }

                AlignmentMatch am = AlignmentMatch(std::string(refChr), std::string(queryChr), referenceStart, referenceEnd, queryStart, queryEnd, score, strand, std::string(c_gene), "");
                v_am.push_back(am);
            }
        }

        infile.close();
        return true;
    }

    return false;
}

int proportionalAlignment(int argc, char **argv) {

    int32_t matchingScore = 0;
    int32_t mismatchingPenalty = -4;
    int32_t openGapPenalty1 = -4;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -80;
    int32_t extendGapPenalty2 = -1;

    double MIN_ALIGNMENT_SCORE = 2;

    int minExon = 20;
    double minimumSimilarity = 0;
    double minimumSimilarity2 = 0;

    int32_t wfaSize3 = 100000; // if the inter-anchor length is shorter than this value, stop trying to find new anchors
    int64_t windowWidth = 100000;
    int expectedCopies = 1;
    double maximumSimilarity = 0.6;

    bool searchForNewAnchors = true;

    double calculateIndelDistance = 3;
    double GAP_OPEN_PENALTY = -0.03;
    double INDEL_SCORE = -0.01;

    int MAX_DIST_BETWEEN_MATCHES = 25;  // between maize and sorghum set it as 25*3000
    int refMaximumTimes = 1;
    int queryMaximumTimes = 2;

    int threads = 1;
    bool exonModel = false;
    std::stringstream usage;

    usage << "Usage: " << PROGRAMNAME
          << " proali -i refGffFile -r refGenome -a cds.sam -as cds.fa -ar ref.sam -s targetGenome -n outputAnchorFile -o output.maf -f output.fragmentation.maf -R 1 -Q 1" << std::endl <<
          "Options" << std::endl <<
          " -h           produce help message" << std::endl <<
          " -i   FILE    reference GFF/GTF file" << std::endl <<
          " -r   FILE    reference genome sequence" << std::endl <<
          " -as  FILE    anchor sequence file. (output from the gff2seq command)" << std::endl <<
          " -a   FILE    sam file by mapping conserved sequence to query genome" << std::endl <<
          " -s   FILE    target genome sequence" << std::endl <<
          " -n   FILE    output anchors file" << std::endl <<
          " -o   FILE    output file in maf format" << std::endl <<
          " -f   FILE    output sequence alignment for each anchor/inter-anchor region in maf format" << std::endl <<
          " -t   INT     number of threads (default: " << threads << ")" << std::endl <<
          " -fa3 INT     if the inter-anchor length is shorter than this value, stop trying to find new anchors (default: " << wfaSize3 << ")" << std::endl <<
          " -w   INT     sequence alignment window width (default: " << windowWidth << ")" << std::endl <<
          " -R   INT     reference genome maximum alignment coverage " << std::endl <<
          " -Q   INT     query genome maximum alignment coverage " << std::endl <<
          " -B   INT     mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1  INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1  INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -O2  INT     open gap penalty 2 (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2  INT     extend gap penalty 2 (default: " << extendGapPenalty2 << ")" << std::endl <<
          " -m   INT     minimum exon length to use (default: " << minExon << ", should be identical with the setting of gff2seq function)" << std::endl <<
          " -mi  DOUBLE  minimum full-length CDS anchor hit similarity to use (default:" << minimumSimilarity << ")" << std::endl <<
          " -mi2 DOUBLE  minimum novel anchor hit similarity to use (default:" << minimumSimilarity2 << ")" << std::endl <<
          " -e   INT     maximum expected copy number of each gene on each chromosome (default: " << expectedCopies << ")" << std::endl << // this is used to duplicated anchors from the sam file
          "              This prevents using tandem duplicated genes to identify collinear block" << std::endl <<
          " -y   DOUBLE  minimal ratio of e+1 similarity to 1 similarity to drop an anchor (default: " << maximumSimilarity << ")" << std::endl <<
          " -ar  FILE    sam file by mapping conserved sequence to reference genome" << std::endl <<
          "              this is used to improve the accuracy of anchors mapping" << std::endl <<
          " Following parameters are to identify collinear blocks" << std::endl <<
          " -d   DOUBLE  calculate IndelDistance (default: " << calculateIndelDistance << ")" << std::endl <<
          " -O   DOUBLE  chain open gap penalty (default: " << GAP_OPEN_PENALTY << ")" << std::endl <<
          " -E   DOUBLE  chain extend gap penalty (default: " << INDEL_SCORE << ")" << std::endl <<
          " -I   DOUBLE  minimum chain score (default: " << MIN_ALIGNMENT_SCORE << ")" << std::endl <<
          " -D   INT     maximum gap size for chain (default: " << MAX_DIST_BETWEEN_MATCHES << ")" << std::endl <<
          " -ns          do not search for new anchors (default: false)" << std::endl <<
          " -x           use exon records instead of CDS from the GFF file (should be identical with the setting of gff2seq function)" << std::endl <<
          std::endl;

    InputParser inputParser(argc, argv);
    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help")) {
        std::cerr << usage.str();
        return 1;
    } else if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
               inputParser.cmdOptionExists("-a") && inputParser.cmdOptionExists("-s") && inputParser.cmdOptionExists("-as")
               && (inputParser.cmdOptionExists("-n") || inputParser.cmdOptionExists("-f") || inputParser.cmdOptionExists("-o") ) ) {

        std::string refGffFilePath = inputParser.getCmdOption("-i");
        std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
        std::string cdsSequenceFile = inputParser.getCmdOption("-as");
        std::string samFilePath = inputParser.getCmdOption("-a");
        std::string targetGenomeSequence = inputParser.getCmdOption("-s");

        std::string outPutMafFile;
        if (inputParser.cmdOptionExists("-o")) {
            outPutMafFile = inputParser.getCmdOption("-o");
        }

        std::string outPutFragedFile;
        if (inputParser.cmdOptionExists("-f")) {
            outPutFragedFile = inputParser.getCmdOption("-f");
        }

        if (inputParser.cmdOptionExists("-t")) {
            threads = std::stoi(inputParser.getCmdOption("-t"));
        }

        if (inputParser.cmdOptionExists("-fa3")) {
            wfaSize3 = std::stoi(inputParser.getCmdOption("-fa3"));
        }

        if (inputParser.cmdOptionExists("-w")) {
            windowWidth = std::stoi(inputParser.getCmdOption("-w"));
        }

        if (inputParser.cmdOptionExists("-R")) {
            refMaximumTimes = std::stoi(inputParser.getCmdOption("-R"));
        } else {
            std::cerr << "parameter -R is required" << std::endl;
            std::cerr << usage.str();
            return 1;
        }

        if (inputParser.cmdOptionExists("-Q")) {
            queryMaximumTimes = std::stoi(inputParser.getCmdOption("-Q"));
        } else {
            std::cerr << "parameter -Q is required" << std::endl;
            std::cerr << usage.str();
            return 1;
        }

        if (inputParser.cmdOptionExists("-B")) {
            mismatchingPenalty = std::stoi(inputParser.getCmdOption("-B"));
            if (mismatchingPenalty >= 0) {
                std::cout << "parameter of B should be a negative value" << std::endl;
                return 1;
            }
        }
        if (inputParser.cmdOptionExists("-O1")) {
            openGapPenalty1 = std::stoi(inputParser.getCmdOption("-O1"));
            if (openGapPenalty1 >= 0) {
                std::cout << "parameter of O1 should be a negative value" << std::endl;
                return 1;
            }
        }
        if (inputParser.cmdOptionExists("-E1")) {
            extendGapPenalty1 = std::stoi(inputParser.getCmdOption("-E1"));
            if (extendGapPenalty1 >= 0) {
                std::cout << "parameter of E1 should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-O2")) {
            openGapPenalty2 = std::stoi(inputParser.getCmdOption("-O2"));
            if (openGapPenalty2 >= 0) {
                std::cout << "parameter of O1 should be a negative value" << std::endl;
                return 1;
            }
        }
        if (inputParser.cmdOptionExists("-E2")) {
            extendGapPenalty2 = std::stoi(inputParser.getCmdOption("-E2"));
            if (extendGapPenalty2 > 0) {
                std::cout << "parameter of E1 should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-m")) {
            minExon = std::stoi(inputParser.getCmdOption("-m"));
        }
        if (inputParser.cmdOptionExists("-mi")) {
            minimumSimilarity = std::stod(inputParser.getCmdOption("-mi"));
        }
        if (inputParser.cmdOptionExists("-mi2")) {
            minimumSimilarity2 = std::stod(inputParser.getCmdOption("-mi2"));
        }

        if (inputParser.cmdOptionExists("-e")) {
            expectedCopies = std::stoi(inputParser.getCmdOption("-e"));
        }

        if (inputParser.cmdOptionExists("-y")) {
            maximumSimilarity = std::stod(inputParser.getCmdOption("-y"));
        }

        std::string referenceSamFilePath;
        if (inputParser.cmdOptionExists("-ar")) {
            referenceSamFilePath = inputParser.getCmdOption("-ar");
        }

        if (inputParser.cmdOptionExists("-d")) {
            calculateIndelDistance = std::stod(inputParser.getCmdOption("-d"));
        }
        if (inputParser.cmdOptionExists("-O")) {
            GAP_OPEN_PENALTY = std::stod(inputParser.getCmdOption("-O"));
        }
        if (inputParser.cmdOptionExists("-E")) {
            INDEL_SCORE = std::stod(inputParser.getCmdOption("-E"));
        }
        if (inputParser.cmdOptionExists("-I")) {
            MIN_ALIGNMENT_SCORE = std::stod(inputParser.getCmdOption("-I"));
        }
        if (inputParser.cmdOptionExists("-D")) {
            MAX_DIST_BETWEEN_MATCHES = std::stoi(inputParser.getCmdOption("-D"));
        }

        if (inputParser.cmdOptionExists("-ns")) {
            searchForNewAnchors = false;
        }

        if (inputParser.cmdOptionExists("-ua")) {
            GAP_OPEN_PENALTY = -4;
            INDEL_SCORE = -2;
            MAX_DIST_BETWEEN_MATCHES = 25;
        }

        std::vector<std::vector<AlignmentMatch>> v_v_am;

        std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
        readFastaFile(referenceGenomeSequence, map_ref);

        std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
        readFastaFile(targetGenomeSequence, map_qry);

        if (inputParser.cmdOptionExists("-x")) {
            exonModel = true;
        }

        std::string path_anchors = inputParser.getCmdOption("-n");

        bool flag_CanGen = false;
        std::ifstream infile_anchors(path_anchors);
        if (infile_anchors.good()) {
            std::string wholeCommand = argv[0];
            for (int i = 1; i < argc; ++i) {
                wholeCommand = wholeCommand + " " + argv[i];
            }
            flag_CanGen = genProVectorFromFile(path_anchors, wholeCommand, v_v_am);
        }

        if(!flag_CanGen) {
            std::cout << "setupAnchorsWithSpliceAlignmentResultQuota begin!" << std::endl;
            // generate alignmentMatchsMap.
            setupAnchorsWithSpliceAlignmentResultQuota(refGffFilePath, samFilePath, cdsSequenceFile, v_v_am, INDEL_SCORE, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE,
                                                       MAX_DIST_BETWEEN_MATCHES, refMaximumTimes, queryMaximumTimes,
                                                       calculateIndelDistance, minExon, windowWidth, minimumSimilarity, minimumSimilarity2,
                                                       map_ref, map_qry,
                                                       expectedCopies, wfaSize3, maximumSimilarity, referenceSamFilePath,
                                                       searchForNewAnchors, exonModel);

            // generated anchors file.
            if (inputParser.cmdOptionExists("-n")) {
                std::string wholeCommand = argv[0];
                for (int i = 1; i < argc; ++i) {
                    wholeCommand = wholeCommand + " " + argv[i];
                }

                std::ofstream ofile;
                ofile.open(inputParser.getCmdOption("-n"));
                ofile << "#" << PROGRAMNAME << " " << wholeCommand << std::endl;
                ofile << "refChr" << "\t"
                      << "referenceStart" << "\t"
                      << "referenceEnd" << "\t"
                      << "queryChr" << "\t"
                      << "queryStart" << "\t"
                      << "queryEnd" << "\t"
                      << "strand" << "\t"
                      << "gene" << "\t"
                      << "blockIndex" << "\t"
                      << "score" << std::endl;

                size_t totalAnchors = 0;
                int blockIndex = 0;
                for (std::vector <AlignmentMatch> v_am: v_v_am) {
                    ofile << "#block begin" << std::endl;
                    blockIndex++;

                    for (size_t i = 0; i < v_am.size(); i++) {
                        std::string thisStrand = "+";
                        if (v_am[i].getStrand() == NEGATIVE) {
                            thisStrand = "-";
                        }

                        if (i > 0) {
                            if (v_am[i].getStrand() == POSITIVE && v_am[i - 1].getStrand() == POSITIVE) {
                                ofile << v_am[i].getRefChr() << "\t"
                                      << v_am[i - 1].getRefEndPos() + 1 << "\t"
                                      << v_am[i].getRefStartPos() - 1 << "\t"
                                      << v_am[i].getQueryChr() << "\t"
                                      << v_am[i - 1].getQueryEndPos() + 1 << "\t"
                                      << v_am[i].getQueryStartPos() - 1 << "\t"
                                      << "+" << "\t"
                                      << "interanchor" << "\t"
                                      << blockIndex << "\tNA"
                                      << std::endl;
                            } else if (v_am[i].getStrand() == NEGATIVE && v_am[i - 1].getStrand() == NEGATIVE) {
                                ofile << v_am[i].getRefChr() << "\t"
                                      << v_am[i - 1].getRefEndPos() + 1 << "\t"
                                      << v_am[i].getRefStartPos() - 1 << "\t"
                                      << v_am[i].getQueryChr() << "\t"
                                      << v_am[i].getQueryEndPos() + 1 << "\t"
                                      << v_am[i - 1].getQueryStartPos() - 1 << "\t"
                                      << "-" << "\t"
                                      << "interanchor" << "\t"
                                      << blockIndex << "\tNA"
                                      << std::endl;
                            }
                        }

                        ofile << v_am[i].getRefChr() << "\t"
                              << v_am[i].getRefStartPos() << "\t"
                              << v_am[i].getRefEndPos() << "\t"
                              << v_am[i].getQueryChr() << "\t"
                              << v_am[i].getQueryStartPos() << "\t"
                              << v_am[i].getQueryEndPos() << "\t"
                              << thisStrand << "\t"
                              << v_am[i].getReferenceGeneName() << "\t";

                        if (v_am[i].getReferenceGeneName().find("localAlignment") == std::string::npos) {
                            totalAnchors++;
                            ofile << blockIndex << "\t" << v_am[i].getScore() << std::endl;
                        } else {
                            ofile << blockIndex << "\t" << "NA" << std::endl;
                        }
                    }
                    ofile << "#block end" << std::endl;
                }

                ofile.close();
                std::cout << "totalAnchors:" << totalAnchors << std::endl;
            }

            std::cout << "anchors generate done!" << std::endl;
        }

        if (inputParser.cmdOptionExists("-f") || inputParser.cmdOptionExists("-o") || inputParser.cmdOptionExists("-l")) {
            genomeAlignment(v_v_am, referenceGenomeSequence, targetGenomeSequence, windowWidth,
                            outPutMafFile, outPutFragedFile, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1,
                            openGapPenalty2, extendGapPenalty2, threads);
            std::cout << "AnchorWave done!" << std::endl;
        }
    }
    else {
        std::cerr << usage.str();
        return 1;
    }

    return 0;
}

int ali(int argc, char **argv) {

    int32_t matchingScore = 0;
    int32_t mismatchingPenalty = -6;
    int32_t openGapPenalty1 = -8;
    int32_t extendGapPenalty1 = -2;

    int32_t openGapPenalty2 = -75;
    int32_t extendGapPenalty2 = -1;
    int64_t windowWidth = 100000;

    std::stringstream usage;
    usage << "Usage: " << PROGRAMNAME <<
          " ali -r refSeq.fa -s querySeq.fa" << std::endl <<
          "Options" << std::endl <<
          " -h           produce help message" << std::endl <<
          " -r   FILE    reference sequence (single sequence in FASTA format)" << std::endl <<
          " -s   FILE    target sequence (single sequence in FASTA format)" << std::endl <<
          " -w   INT     sequence alignment window width (default: " << windowWidth << ")" << std::endl <<
          " -B   INT     mismatching penalty (default: " << mismatchingPenalty << ")" << std::endl <<
          " -O1  INT     open gap penalty (default: " << openGapPenalty1 << ")" << std::endl <<
          " -E1  INT     extend gap penalty (default: " << extendGapPenalty1 << ")" << std::endl <<
          " -O2  INT     open gap penalty 2 (default: " << openGapPenalty2 << ")" << std::endl <<
          " -E2  INT     extend gap penalty 2 (default: " << extendGapPenalty2 << ")" << std::endl << std::endl;

    InputParser inputParser(argc, argv);
    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help") || !inputParser.cmdOptionExists("-r") || !inputParser.cmdOptionExists("-s")) {
        std::cerr << usage.str();
        return 1;
    }

    std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
    std::string targetGenomeSequence = inputParser.getCmdOption("-s");

    if (inputParser.cmdOptionExists("-w")) {
        windowWidth = std::stoi(inputParser.getCmdOption("-w"));
    }

    if (inputParser.cmdOptionExists("-B")) {
        mismatchingPenalty = std::stoi(inputParser.getCmdOption("-B"));
        if (mismatchingPenalty >= 0) {
            std::cout << "parameter of B should be a negative value" << std::endl;
            return 1;
        }
    }

    if (inputParser.cmdOptionExists("-O1")) {
        openGapPenalty1 = std::stoi(inputParser.getCmdOption("-O1"));
        if (openGapPenalty1 >= 0) {
            std::cout << "parameter of O1 should be a negative value" << std::endl;
            return 1;
        }
    }
    if (inputParser.cmdOptionExists("-E1")) {
        extendGapPenalty1 = std::stoi(inputParser.getCmdOption("-E1"));
        if (extendGapPenalty1 >= 0) {
            std::cout << "parameter of E1 should be a negative value" << std::endl;
            return 1;
        }
    }

    if (inputParser.cmdOptionExists("-O2")) {
        openGapPenalty2 = std::stoi(inputParser.getCmdOption("-O2"));
        if (openGapPenalty2 >= 0) {
            std::cout << "parameter of O1 should be a negative value" << std::endl;
            return 1;
        }
    }
    if (inputParser.cmdOptionExists("-E2")) {
        extendGapPenalty2 = std::stoi(inputParser.getCmdOption("-E2"));
        if (extendGapPenalty2 > 0) {
            std::cout << "parameter of E1 should be a negative value" << std::endl;
            return 1;
        }
    }

    std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
    readFastaFile(referenceGenomeSequence, map_ref);

    if (map_ref.size() != 1) {
        std::cerr << "There should be one and only one sequence in the reference FASTA file" << std::endl;
    }

    std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
    readFastaFile(targetGenomeSequence, map_qry);

    if (map_qry.size() != 1) {
        std::cerr << "There should be one and only one sequence in the query FASTA file" << std::endl;
    }

    std::string _alignment_q;
    std::string _alignment_d;

    std::string refSeqStr = getSubsequence2(map_ref, map_ref.begin()->first);
    std::string querySeqStr = getSubsequence2(map_qry, map_qry.begin()->first);

    alignSlidingWindow(querySeqStr, refSeqStr, _alignment_q, _alignment_d, windowWidth, matchingScore, mismatchingPenalty, openGapPenalty1, extendGapPenalty1, openGapPenalty2, extendGapPenalty2);

    std::cout << ">" << map_ref.begin()->first << std::endl;
    std::cout << _alignment_d << std::endl;
    std::cout << ">" << map_qry.begin()->first << std::endl;
    std::cout << _alignment_q << std::endl;

    return 0;
}




int geno(int argc, char **argv) {
/*
    std::stringstream usage;

    int minExon = 20;
    double minimumSimilarity = 0.95;
    double minimumSimilarity2 = 0.2;

    double inversion_PENALTY = -1;
    double MIN_ALIGNMENT_SCORE = 2;
    bool considerInversion = false;

    int32_t wfaSize3 = 100000; // if the inter-anchor length is shorter than this value, stop trying to find new anchors
    int64_t windowWidth = 100000;
    int expectedCopies = 1;
    double maximumSimilarity = 0.6; // the maximum simalarity between secondary hist the primary hit. If the second hit is too similary with primary hit, that is unwanted duplications

    bool searchForNewAnchors = true;

    int threads = 1;
    bool exonModel = false;
    std::string blastpResultFile = "";

    usage << "Usage: " << PROGRAMNAME
          << " genoAli -i refGffFile -r refGenome -a cds.sam -as cds.fa -ar ref.sam -s targetGenome -n outputAnchorFile -o output.maf -f output.fragmentation.maf " << std::endl <<
          "Options" << std::endl <<
          " -h           produce help message" << std::endl <<
          " -i   FILE    reference GFF/GTF file" << std::endl <<
          " -r   FILE    reference genome sequence file in fasta format" << std::endl <<
          " -as  FILE    anchor sequence file. (output from the gff2seq command)" << std::endl <<
          " -a   FILE    sam file generated by mapping conserved sequence to query genome" << std::endl <<
          " -s   FILE    target genome sequence file in fasta format" << std::endl <<
          "              Those sequences with the same name in the reference genome and query genome file would be aligned" << std::endl <<
          " -n   FILE    output anchors file" << std::endl <<
          " -o   FILE    output file in maf format" << std::endl <<
          " -f   FILE    output sequence alignment for each anchor/inter-anchor region in maf format" << std::endl <<
          " -t   INT     number of threads (default: " << threads << ")" << std::endl <<
          " -m   INT     minimum exon length to use (default: " << minExon << ", should be identical with the setting of gff2seq function)" << std::endl <<
          " -mi  DOUBLE  minimum full-length CDS anchor hit similarity to use (default:" << minimumSimilarity << ")" << std::endl <<
          " -mi2 DOUBLE  minimum novel anchor hit similarity to use (default:" << minimumSimilarity2 << ")" << std::endl <<
          " -ar  FILE    sam file generated by mapping conserved sequence to reference genome" << std::endl <<
          " -w   INT     sequence alignment window width (default: " << windowWidth << ")" << std::endl <<
          " -IV          whether to call inversions (default: false)" << std::endl <<
          " -IC  DOUBLE  penalty for having a non-linear match in inversion region (default: " << inversion_PENALTY << ")" << std::endl <<
          "              We use IC * alignment_similarity as the penalty in the inversion block" << std::endl <<
          " -I   DOUBLE  minimum score to keep an inversion (default: " << MIN_ALIGNMENT_SCORE << ")" << std::endl <<
          " -e   INT     maximum expected copy number of each gene on each chromosome (default: " << expectedCopies << ")" << std::endl << // this is used to duplicated anchors from the sam file
          "              This prevents using tandem duplicated genes to identify collinear block" << std::endl <<
          " -y   DOUBLE  minimal ratio of e+1 similarity to 1 similarity to drop an anchor (default: " << maximumSimilarity << ")" << std::endl <<
          " -ns          do not search for new anchors (default: false)" << std::endl <<
          " -x           use exon records instead of CDS from the GFF file (should be identical with the setting of gff2seq function)" << std::endl
          << std::endl;

    InputParser inputParser(argc, argv);

    if (inputParser.cmdOptionExists("-h") || inputParser.cmdOptionExists("--help")) {
        std::cerr << usage.str();
    }
    else if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") &&
             inputParser.cmdOptionExists("-a") && inputParser.cmdOptionExists("-s") && inputParser.cmdOptionExists("-as") &&
             (inputParser.cmdOptionExists("-n") || inputParser.cmdOptionExists("-f") || inputParser.cmdOptionExists("-o") || inputParser.cmdOptionExists("-l"))) {

        std::string refGffFilePath = inputParser.getCmdOption("-i");
        std::string path_ref_GenomeSequence = inputParser.getCmdOption("-r");
        std::string cdsSequenceFile = inputParser.getCmdOption("-as");
        std::string samFilePath = inputParser.getCmdOption("-a");
        std::string path_target_GenomeSequence = inputParser.getCmdOption("-s");
        std::string outPutMafFile;

        if (inputParser.cmdOptionExists("-o")) {
            outPutMafFile = inputParser.getCmdOption("-o");
        }

        std::string outPutFragedFile;
        if (inputParser.cmdOptionExists("-f")) {
            outPutFragedFile = inputParser.getCmdOption("-f");
        }

        if (inputParser.cmdOptionExists("-t")) {
            threads = std::stoi(inputParser.getCmdOption("-t"));
        }

        if (inputParser.cmdOptionExists("-m")) {
            minExon = std::stoi(inputParser.getCmdOption("-m"));
        }
        if (inputParser.cmdOptionExists("-mi")) {
            minimumSimilarity = std::stod(inputParser.getCmdOption("-mi"));
        }
        if (inputParser.cmdOptionExists("-mi2")) {
            minimumSimilarity2 = std::stod(inputParser.getCmdOption("-mi2"));
        }

        std::string referenceSamFilePath;
        if (inputParser.cmdOptionExists("-ar")) {
            referenceSamFilePath = inputParser.getCmdOption("-ar");
        }



        if (inputParser.cmdOptionExists("-IV")) {
            considerInversion = true;
        }
        if (inputParser.cmdOptionExists("-IC")) {
            inversion_PENALTY = std::stod(inputParser.getCmdOption("-IC"));
            if (inversion_PENALTY >= 0) {
                std::cout << "parameter of IC should be a negative value" << std::endl;
                return 1;
            }
        }

        if (inputParser.cmdOptionExists("-I")) {
            MIN_ALIGNMENT_SCORE = std::stod(inputParser.getCmdOption("-I"));
        }

        if (MIN_ALIGNMENT_SCORE < 1) {
            std::cout << "minimum score to keep an inversion is not larger than 1, maybe output weired inversions. please change it" << std::endl;
            return 1;
        }
        if (inputParser.cmdOptionExists("-e")) {
            expectedCopies = std::stoi(inputParser.getCmdOption("-e"));
        }

        if (inputParser.cmdOptionExists("-y")) {
            maximumSimilarity = std::stod(inputParser.getCmdOption("-y"));
        }
        if (inputParser.cmdOptionExists("-w")) {
            windowWidth = std::stoi(inputParser.getCmdOption("-w"));
        }

        if (inputParser.cmdOptionExists("-ns")) {
            searchForNewAnchors = false;
        }

        if (inputParser.cmdOptionExists("-x")) {
            exonModel = true;
        }

        std::map<std::string, std::vector<AlignmentMatch>> alignmentMatchsMapT;

        std::ifstream infile(blastpResultFile);
        if (!infile.good()) {
            std::cerr << "error in opening fasta file " << blastpResultFile << std::endl;
            exit(1);
        }
        std::string line;

        while (std::getline(infile, line)) {
            std::vector<std::string> elems;
            char seperator = '\t';
            split(line, seperator, elems);
            std::string refGene = elems[0];
            std::string refChr = elems[1];
            int32_t refStart = stoi(elems[2]);
            int32_t refEnd = stoi(elems[3]);
            std::string refStrand = elems[4];

            std::string queryGene = elems[5];
            std::string queryChr = elems[6];
            int32_t queryStart = stoi(elems[7]);
            int32_t queryEnd = stoi(elems[8]);
            std::string queryStrand = elems[9];

            double alignmentScore = stoi(elems[10]);
            AlignmentMatch orthologPair(refChr, queryChr, startRef + r->rs,
                                        startRef + newAnchorRefEnd, startQuery + r->qs,
                                        startQuery + newAnchorQueryEnd, thisScore, POSITIVE,
                                        alignmentName, alignmentName);


        }

        infile.close();

        if (temp.size() > 0) {
            myAlignmentMatchSort(temp, inversion_PENALTY, MIN_ALIGNMENT_SCORE, keepTandemDuplication, considerInversion);

            std::vector<AlignmentMatch> sortedAlignmentMatchs;

            if (temp.size() > 1) {
                longestPath(temp, sortedAlignmentMatchs, keepTandemDuplication, MIN_ALIGNMENT_SCORE);
            } else {
                sortedAlignmentMatchs = temp;
            }

            for (unsigned long i = 0; i < sortedAlignmentMatchs.size(); ++i) {
                map_v_am[it->first].push_back(sortedAlignmentMatchs[i]);
            }
        }




        std::map<std::string, std::tuple<std::string, long, long, int> > map_ref;
        readFastaFile(path_ref_GenomeSequence, map_ref);

        std::map<std::string, std::tuple<std::string, long, long, int> > map_qry;
        readFastaFile(path_target_GenomeSequence, map_qry);

        std::cout << "setupAnchorsWithSpliceAlignmentResult begin " << std::endl;
        std::map<std::string, std::vector<AlignmentMatch>> map_v_am;

        std::string path_anchors = inputParser.getCmdOption("-n");

        bool flag_CanGen = false;
        std::ifstream infile_anchors(path_anchors);
        if (infile_anchors.good()) {
            std::string wholeCommand = argv[0];
            for (int i = 1; i < argc; ++i) {
                wholeCommand = wholeCommand + " " + argv[i];
            }

            flag_CanGen = genGenoMapFromFile(path_anchors, wholeCommand, map_v_am);
            std::cout << "map_v_am generated!" << std::endl;
        }

        if(!flag_CanGen) {
            setupAnchorsWithSpliceAlignmentResult(refGffFilePath, cdsSequenceFile, samFilePath, map_v_am,
                                                  inversion_PENALTY, MIN_ALIGNMENT_SCORE, considerInversion, minExon, windowWidth, minimumSimilarity, minimumSimilarity2,
                                                  map_ref, map_qry,
                                                  expectedCopies, maximumSimilarity, referenceSamFilePath, wfaSize3, searchForNewAnchors, exonModel);

            std::cout << "setupAnchorsWithSpliceAlignmentResult done!" << std::endl;

            for (std::map < std::string, std::vector < AlignmentMatch >> ::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
                myAlignmentMatchSort(it->second, inversion_PENALTY, MIN_ALIGNMENT_SCORE, false, false);
            }

            std::cout << "myAlignmentMatchSort done!" << std::endl;

            if (inputParser.cmdOptionExists("-n")) {
                std::string wholeCommand = argv[0];
                for (int i = 1; i < argc; ++i) {
                    wholeCommand = wholeCommand + " " + argv[i];
                }

                std::ofstream ofile;
                ofile.open(inputParser.getCmdOption("-n"));
                ofile << "#" << PROGRAMNAME << " " << wholeCommand << std::endl;
                int blockIndex = 0;
                ofile << "refChr" << "\t"
                      << "referenceStart" << "\t"
                      << "referenceEnd" << "\t"
                      << "queryChr" << "\t"
                      << "queryStart" << "\t"
                      << "queryEnd" << "\t"
                      << "strand" << "\t"
                      << "gene" << "\t"
                      << "blockIndex" << "\tscore" << std::endl;

                for (std::map < std::string, std::vector < AlignmentMatch >> ::iterator it = map_v_am.begin(); it != map_v_am.end(); ++it) {
                    ofile << "#block begin" << std::endl;
                    std::vector < AlignmentMatch > v_am = it->second;
                    blockIndex++;
                    bool hasInversion = false;
                    for (size_t i = 0; i < v_am.size(); ++i) {
                        if (i > 0) {
                            if (v_am[i].getStrand() == POSITIVE && v_am[i - 1].getStrand() == POSITIVE) {
                                ofile << v_am[i].getRefChr() << "\t"
                                      << v_am[i - 1].getRefEndPos() + 1 << "\t"
                                      << v_am[i].getRefStartPos() - 1 << "\t"
                                      << v_am[i].getQueryChr() << "\t"
                                      << v_am[i - 1].getQueryEndPos() + 1 << "\t"
                                      << v_am[i].getQueryStartPos() - 1 << "\t"
                                      << "+" << "\t"
                                      << "interanchor" << "\t"
                                      << blockIndex << "\tNA" << std::endl;
                            } else if (v_am[i].getStrand() == NEGATIVE && v_am[i - 1].getStrand() == NEGATIVE
                                       && v_am[i - 1].getRefEndPos() < v_am[i].getRefStartPos()
                                       && v_am[i - 1].getQueryStartPos() > v_am[i].getQueryEndPos()) {

                                ofile << v_am[i].getRefChr() << "\t"
                                      << v_am[i - 1].getRefEndPos() + 1 << "\t"
                                      << v_am[i].getRefStartPos() - 1 << "\t"
                                      << v_am[i].getQueryChr() << "\t"
                                      << v_am[i].getQueryEndPos() + 1 << "\t"
                                      << v_am[i - 1].getQueryStartPos() - 1 << "\t"
                                      << "-" << "\t"
                                      << "interanchor" << "\t"
                                      << blockIndex << "\tNA" << std::endl;
                            }
                        }

                        std::string thisStrand = "+";
                        if (v_am[i].getStrand() == NEGATIVE) {
                            thisStrand = "-";
                            hasInversion = true;
                        }

                        ofile <<  v_am[i].getRefChr() << "\t"
                              << v_am[i].getRefStartPos() << "\t"
                              << v_am[i].getRefEndPos() << "\t"
                              << v_am[i].getQueryChr() << "\t"
                              << v_am[i].getQueryStartPos() << "\t"
                              << v_am[i].getQueryEndPos() << "\t"
                              << thisStrand << "\t"
                              << v_am[i].getReferenceGeneName() << "\t";

                        if (v_am[i].getReferenceGeneName().find("localAlignment") == std::string::npos) {
                            ofile << blockIndex << "\t" << v_am[i].getScore() << std::endl;
                        } else {
                            ofile << blockIndex << "\t" << "NA" << std::endl;
                        }
                    }

                    int i = v_am.size() - 1;
                    if (!hasInversion) {
                        size_t size_sr2 = getSequenceSizeFromPath2(map_ref[v_am[i].getRefChr()]);
                        size_t size_sq2 = getSequenceSizeFromPath2(map_qry[v_am[i].getQueryChr()]);

                        ofile << v_am[i].getRefChr() << "\t"
                              << v_am[i].getRefEndPos() + 1 << "\t"
                              << size_sr2 << "\t"
                              << v_am[i].getQueryChr() << "\t"
                              << v_am[i].getQueryEndPos() + 1 << "\t"
                              << size_sq2 << "\t"
                              << "+" << "\t"
                              << "interanchor" << "\t"
                              << blockIndex << "\tNA" << std::endl;
                    }

                    ofile << "#block end" << std::endl;
                }
                ofile.close();
            }

            std::cout << "anchors generate done!" << std::endl;
        }

    }
    else {
        std::cerr << usage.str();
    }
 */
    return 0;
}




int pro(int argc, char **argv) {
    double MIN_ALIGNMENT_SCORE = 2;
    double GAP_OPEN_PENALTY = -0.02;
    double GAP_EXTENSION_PENALTY = -0.005;
    int MAX_DIST_BETWEEN_MATCHES = 25;
    int refMaximumTimes = 1;
    int queryMaximumTimes = 1;
    int OVER_LAP_WINDOW = 5;
    int delete_tandem = 0;
    int get_all_collinear_gene_pair = 0;
    int count_style = 1;

    std::stringstream usage;

    usage << "Usage: " << PROGRAMNAME
          << " pro -i inputFile -n outputAnchorFile -R refMaximumTimes -Q queryMaximumTimes" << std::endl <<
          "Options" << std::endl <<
          " -h           produce help message" << std::endl <<
          " -i   FILE    input file" << std::endl <<
          " -o   FILE    output anchors file" << std::endl <<
          " -R   INT     reference genome maximum alignment coverage " << std::endl <<
          " -Q   INT     query genome maximum alignment coverage " << std::endl <<
          " -O   DOUBLE  chain open gap penalty (default: " << GAP_OPEN_PENALTY << ")" << std::endl <<
          " -E   DOUBLE  chain gap extend penalty (default: " << GAP_EXTENSION_PENALTY << ")" << std::endl <<
          " -I   DOUBLE  minimum chain score (default: " << MIN_ALIGNMENT_SCORE << ")" << std::endl <<
          " -D   INT     maximum gap size for chain (default: " << MAX_DIST_BETWEEN_MATCHES << ")" << std::endl <<
          " -m   INT     Specify whether to delete tandem gene pairs in the input file (default: " << delete_tandem << ")" << std::endl <<
          "              Options: 0 to retain tandem gene pairs; 1 or any other integer to delete them." << std::endl <<
          " -W   INT     When -m is set to delete tandem gene pairs(1 or any other integer), " << std::endl <<
          "              specify the maximum distance allowed between two homologous gene pairs before they are considered for deletion" << "(default: " << OVER_LAP_WINDOW << ")" << std::endl <<
          "              This parameter is ignored if -m is not set to delete tandem gene pairs(set 0)." << std::endl <<
          " -a   INT     Enable this flag to get all collinear results and disable the R and Q parameters (default: " << get_all_collinear_gene_pair << ")" << std::endl <<
          "              Options: 0:enable R Q parameter; 1 or other integer: all collinear result" << std::endl <<
          " -c   INT     For R Q count, 0: gene number count; 1 or other integer: block count (default: " << count_style << ")"
          << std::endl;

    InputParser inputParser(argc, argv);
    if (inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-o")) {
        std::string fakeBlastpResultFile = inputParser.getCmdOption("-i");
        std::string path_anchors = inputParser.getCmdOption("-o");
        std::ofstream outputFile_test(path_anchors);          // avoid output path exists spelling error
        if (!outputFile_test.good()) {
            std::cerr << "error in creating outputResultFile file, please check your output path :" << path_anchors
                      << std::endl;
            exit(1);
        }
        outputFile_test.close();

        if (inputParser.cmdOptionExists("-R")) {
            refMaximumTimes = std::stoi(inputParser.getCmdOption("-R"));
        }
        if (inputParser.cmdOptionExists("-Q")) {
            queryMaximumTimes = std::stoi(inputParser.getCmdOption("-Q"));
        }
        if (inputParser.cmdOptionExists("-a")) {
            get_all_collinear_gene_pair = std::stoi(inputParser.getCmdOption("-a"));
        }

        if (inputParser.cmdOptionExists("-O")) {
            GAP_OPEN_PENALTY = std::stod(inputParser.getCmdOption("-O"));
        }
        if (inputParser.cmdOptionExists("-E")) {
            GAP_EXTENSION_PENALTY = std::stod(inputParser.getCmdOption("-E"));
        }
        if (inputParser.cmdOptionExists("-I")) {
            MIN_ALIGNMENT_SCORE = std::stod(inputParser.getCmdOption("-I"));
        }
        if (inputParser.cmdOptionExists("-D")) {
            MAX_DIST_BETWEEN_MATCHES = std::stoi(inputParser.getCmdOption("-D"));
        }
        if (inputParser.cmdOptionExists("-W")) {
            OVER_LAP_WINDOW = std::stoi(inputParser.getCmdOption("-W"));
        }
        if (inputParser.cmdOptionExists("-m")) {
            delete_tandem = std::stoi(inputParser.getCmdOption("-m"));
        }


        std::vector<std::vector<AlignmentMatch>> alignmentMatchsMap;
        std::vector<AlignmentMatch> alignmentMatchsMapT;
        std::ifstream infile(fakeBlastpResultFile);
        if (!infile.good()) {
            std::cerr << "error in opening fakeBlastpResultFile file:" << fakeBlastpResultFile << std::endl;
            exit(1);
        }
        std::string line;
        // prepare data in RAM begin
        while (std::getline(infile, line)) {
            std::vector<std::string> elems;
            char separator = '\t';
            split(line, separator, elems);
            std::string refGene = elems[0];
            std::string refChr = elems[1];
            int refId = std::stoi(elems[2]);
            int32_t refStart = stoi(elems[3]);
            int32_t refEnd = stoi(elems[4]);
            std::string refStrand = elems[5];
            std::string queryGene = elems[6];
            std::string queryChr = elems[7];
            int queryId = std::stoi(elems[8]);
            int32_t queryStart = stoi(elems[9]);
            int32_t queryEnd = stoi(elems[10]);
            std::string queryStrand = elems[11];
            double alignmentScore = stod(elems[12])/100.0;
            STRAND thisStrand = POSITIVE;
            if (refStrand != queryStrand) {
                thisStrand = NEGATIVE;
            }

            if (refGene == queryGene) {
                continue;
            }

            AlignmentMatch orthologPair(refChr, queryChr, refStart,
                                        refEnd, queryStart,
                                        queryEnd, alignmentScore, thisStrand,
                                        refGene, queryGene);
            orthologPair.setRefId(refId);
            orthologPair.setQueryId(queryId);
            alignmentMatchsMapT.push_back(orthologPair);
        }
        infile.close();
        // prepare data in RAM end
        if (delete_tandem !=0 ) {
            // sort by querychr refchr querystart(query gene name) and identity/100,respectively.
            orthologPairSortQuery(alignmentMatchsMapT);
            // filter ref tandem gene by threshold five. This thought comes from MCScanX.
            std::vector<AlignmentMatch>::const_iterator it1, prev_pair1;
            std::vector<AlignmentMatch> alignmentMatchsMapT_cpy1;
            std::vector<AlignmentMatch> match_bin1;
            prev_pair1 = it1 = alignmentMatchsMapT.begin();
            it1++;
            // insert first pair
            match_bin1.push_back(*prev_pair1);
            // match_bin1[0] has a maximum identity value in the match_bin by orthologPairSortQuery function.
//            std::cout << "\t" << "line1606" << "\t" << alignmentMatchsMapT.size() << std::endl;
            for (; it1 != alignmentMatchsMapT.end(); it1++) {
                if (it1->getRefChr() == prev_pair1->getRefChr() && it1->getQueryChr() == prev_pair1->getQueryChr()) {
                    if ((it1->getQueryId() != prev_pair1->getQueryId()) ||
                        (std::abs(it1->getRefId() - prev_pair1->getRefId()) > OVER_LAP_WINDOW)) {
                        alignmentMatchsMapT_cpy1.push_back(match_bin1[0]);
                        match_bin1.clear();
                    }
                    match_bin1.push_back(*it1);
                } else if (it1->getRefChr() != prev_pair1->getRefChr() ||
                           it1->getQueryChr() != prev_pair1->getQueryChr()) {
                    alignmentMatchsMapT_cpy1.push_back(match_bin1[0]);
                    match_bin1.clear();
                    prev_pair1 = it1;
                    match_bin1.push_back(*prev_pair1);
                    continue;
                }
                prev_pair1 = it1;
            }
            alignmentMatchsMapT_cpy1.push_back(match_bin1[0]);  //the last bin
            alignmentMatchsMapT.clear();
            alignmentMatchsMapT = alignmentMatchsMapT_cpy1;
//            std::cout << "\t" << "line1628" << "\t" << alignmentMatchsMapT.size() << std::endl;
            alignmentMatchsMapT_cpy1.clear();

            // sort by querychr refchr refstart(ref gene name) and identity/100,respectively.
            orthologPairSortReference(alignmentMatchsMapT);
            // filter query tandem gene by threshold five. This thought comes from MCScanX.
            std::vector<AlignmentMatch>::const_iterator it2, prev_pair2;
            std::vector<AlignmentMatch> alignmentMatchsMapT_cpy2;
            std::vector<AlignmentMatch> match_bin2;
            prev_pair2 = it2 = alignmentMatchsMapT.begin();
            it2++;
            // insert first pair
            match_bin2.push_back(*prev_pair2);
            // match_bin2[0] has a maximum identity value in the match_bin by orthologPairSortReference.
//            std::cout << "\t" << "line1642" << "\t" << alignmentMatchsMapT.size() << std::endl;
            for (; it2 != alignmentMatchsMapT.end(); it2++) {
                if (it2->getRefChr() == prev_pair2->getRefChr() && it2->getQueryChr() == prev_pair2->getQueryChr()) {
                    if ((it2->getRefId() != prev_pair2->getRefId()) ||
                        (std::abs(it2->getQueryId() - prev_pair2->getQueryId()) > OVER_LAP_WINDOW)) {
                        alignmentMatchsMapT_cpy2.push_back(match_bin2[0]);
                        match_bin2.clear();
                    }
                    match_bin2.push_back(*it2);
                } else if (it1->getRefChr() != prev_pair1->getRefChr() ||
                           it1->getQueryChr() != prev_pair1->getQueryChr()) {
                    alignmentMatchsMapT_cpy2.push_back(match_bin2[0]);
                    match_bin2.clear();
                    prev_pair2 = it2;
                    match_bin2.push_back(*prev_pair2);
                    continue;
                }
                prev_pair2 = it2;
            }
            alignmentMatchsMapT_cpy2.push_back(match_bin2[0]);
            alignmentMatchsMapT.clear();
            alignmentMatchsMapT = alignmentMatchsMapT_cpy2;
//            std::cout << "\t" << "line1664" << "\t" << alignmentMatchsMapT.size() << std::endl;
            alignmentMatchsMapT_cpy2.clear();
        }

        // begin setting index, they are necessary in the longest path approach
        std::map<std::string, std::map<int, std::string>> queryIndexMap; // chr : index : queryGeneName
        for (const auto &ii: alignmentMatchsMapT) {
            if (queryIndexMap.find(ii.getQueryChr()) == queryIndexMap.end()) {
                queryIndexMap[ii.getQueryChr()] = std::map<int, std::string>();
            } else {
                if (queryIndexMap[ii.getQueryChr()].find(ii.getQueryId())
                    == queryIndexMap[ii.getQueryChr()].end()) {
                    queryIndexMap[ii.getQueryChr()][ii.getQueryId()] =
                            ii.getQueryGeneName();
                }
            }
        }

        std::map<std::string, std::map<int, std::string>> refIndexMap;  // chr : index : refGeneName
        for (const auto &ii: alignmentMatchsMapT) {
            if (refIndexMap.find(ii.getRefChr()) == refIndexMap.end()) {
                refIndexMap[ii.getRefChr()] = std::map<int, std::string>();
            } else {
                if (refIndexMap[ii.getRefChr()].find(ii.getRefId())
                    == refIndexMap[ii.getRefChr()].end()) {
                    refIndexMap[ii.getRefChr()][ii.getRefId()] =
                            ii.getReferenceGeneName();
                }
            }
        }

        orthologPairSortPosition(alignmentMatchsMapT);
//        std::cout << "\t" << "line1739" << "\t" << alignmentMatchsMapT.size() << std::endl;
        std::vector<double> block_score;
//        for (AlignmentMatch &ortholog_pair: alignmentMatchsMapT) {
//            std::cout << "line1746" << "\t" << ortholog_pair.getQueryChr()  << "\t" << ortholog_pair.getQueryGeneName() << "\t"  << ortholog_pair.getQueryId() << "\t" << ortholog_pair.getRefChr() << "\t" << ortholog_pair.getReferenceGeneName() << "\t" << ortholog_pair.getRefId() << "\t" << ortholog_pair.getScore() << std::endl;
//        }
//        exit(1);
        longestPathQuotaGene(alignmentMatchsMapT, alignmentMatchsMap, refIndexMap, queryIndexMap,
                             GAP_EXTENSION_PENALTY, GAP_OPEN_PENALTY, MIN_ALIGNMENT_SCORE, MAX_DIST_BETWEEN_MATCHES,
                             refMaximumTimes, queryMaximumTimes, block_score, count_style, get_all_collinear_gene_pair);
        // output anchors file.
        if (inputParser.cmdOptionExists("-o")) {
            std::string wholeCommand = argv[0];
            for (int i = 1; i < argc; ++i) {
                wholeCommand = wholeCommand.append(" ");
                wholeCommand = wholeCommand.append(argv[i]);
            }
            std::ofstream offile;
            offile.open(inputParser.getCmdOption("-o"));
            offile << "#" << PROGRAMNAME << " " << wholeCommand << std::endl;
            offile << "refGene" << "\t"
                   << "refChr" << "\t"
                   << "refId" << "\t"
                   << "referenceStart" << "\t"
                   << "referenceEnd" << "\t"
                   << "queryGene" << "\t"
                   << "queryChr" << "\t"
                   << "queryId" << "\t"
                   << "queryStart" << "\t"
                   << "queryEnd" << "\t"
                   << "strand" << "\t"
                   << "score"
                   << std::endl;

            size_t totalAnchors = 0;
            int blockIndex = 1;
            for (std::vector<AlignmentMatch> v_am: alignmentMatchsMap) {
//                if (v_am.size() >= 4) {
                    std::string plus_minus = "NEGATIVE";
                    if (v_am[0].getStrand() == POSITIVE) {
                        plus_minus = "POSITIVE";
                    }
                    offile << "##Alignment" << "\t" << blockIndex << "\t" << "N=" << v_am.size() << "\t"
                           << "score=" << block_score[blockIndex - 1] << "\t" << v_am[0].getRefChr() << "&" <<
                           v_am[0].getQueryChr() << "\t" << plus_minus << std::endl;
                    blockIndex++;

                    for (const auto &i: v_am) {
                        std::string thisStrand = "+";
                        if (i.getStrand() == NEGATIVE) {
                            thisStrand = "-";
                        }

                        offile << i.getReferenceGeneName() << "\t"
                               << i.getRefChr() << "\t"
                               << i.getRefId() << "\t"
                               << i.getRefStartPos() << "\t"
                               << i.getRefEndPos() << "\t"
                               << i.getQueryGeneName() << "\t"
                               << i.getQueryChr() << "\t"
                               << i.getQueryId() << "\t"
                               << i.getQueryStartPos() << "\t"
                               << i.getQueryEndPos() << "\t"
                               << thisStrand << "\t"
                               << i.getScore() << std::endl;
                        totalAnchors++;
                    }
//                }
            }
            offile.close();
            std::cout << "totalAnchors:" << totalAnchors << std::endl;
            std::cout << "anchors generate done!" << std::endl;
        }
    } else {
        std::cerr << usage.str();
    }
    return 0;
}
