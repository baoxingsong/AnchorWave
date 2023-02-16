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
    }

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
    int64_t windowWidth = 38000;

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
    readFastaFile(referenceGenomeSequence, map_qry);

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
