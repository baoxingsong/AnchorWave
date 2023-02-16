// =====================================================================================
//
//       Filename:  Transcript.h
//
//    Description:
//
//        Version:  1.0
//        Created:  04/12/2017 02:13:56 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//
// =====================================================================================
/*************************************************************************

 Transcript records in i.e. gff file

*************************************************************************/

#pragma once

#include "GenomeBasicFeature.h"
#include "STRAND.h"
#include <vector>
#include <set>
#include <algorithm>
#include <limits>

class Transcript {
private:
    std::string _name;
    std::set<GenomeBasicFeature> _cdsHashSet;
    std::vector<GenomeBasicFeature> _cdsVector;
    std::set<GenomeBasicFeature> _exonHashSet;
    std::vector<GenomeBasicFeature> _exonVector;
    std::vector<GenomeBasicFeature> _threePrimerUtrVector;
    std::vector<GenomeBasicFeature> _fivePrimerUtrVector;
    std::string _chromeSomeName;
    STRAND _strand;
    int _start;
    int _end;
    int _pStart;
    int _pEnd;
    std::string _metaInformation; //for ORF lift over
    std::string _geneomeSequence;
    std::string _cdsSequence;
    std::string _exonSequence;
    bool _ifOrfShift;
    std::string _source;
    std::string _score;
    std::string _lastColumnInformation; // the last column in the gff/gtf file
    std::string _type;
    bool _ifFivePrimerUtr = false;
    bool _ifThreePrimerUtr = false;
public:
    bool is_ifFivePrimerUtr() const;

    void set_ifFivePrimerUtr(bool _ifFivePrimerUtr);

    bool is_ifThreePrimerUtr() const;

    void set_ifThreePrimerUtr(bool _ifThreePrimerUtr);

    Transcript();

    Transcript(const std::string &name, const std::string &chromeSomeName, const STRAND &strand);

    const std::string &getName() const;

    void setName(const std::string &name);

    std::set<GenomeBasicFeature> &getCdsHashSet();

    std::vector<GenomeBasicFeature> &getCdsVector();

    std::set<GenomeBasicFeature> &getExonHashSet();

    std::vector<GenomeBasicFeature> &getExonVector();

    std::vector<GenomeBasicFeature> &getThreePrimerUtr();

    void addThreePrimerUtr(const GenomeBasicFeature &genomeBasicFeature);

    std::vector<GenomeBasicFeature> &getFivePrimerUtr();

    void addFivePrimerUtr(const GenomeBasicFeature &genomeBasicFeature);

    const STRAND &getStrand() const;

    const std::string &getChromeSomeName() const;

    void addCds(const GenomeBasicFeature &cds);

    void addExon(const GenomeBasicFeature &exon);

    void setStart(const int &start);

    void setEnd(const int &end);

    const int &getStart() const;

    const int &getEnd() const;

    const int &getPStart() const;

    const int &getPEnd() const;

    void setPStart(const int &pStart);

    void setPEnd(const int &pEnd);

    void updateInforCds();

    const std::string &getGeneomeSequence() const;

    void setGeneomeSequence(const std::string &genomeSequence);

    const std::string &getCdsSequence() const;

    void setCdsSequence(const std::string &cdsSequence);

    const std::string &getExonSequence() const;

    void setExonSequence(const std::string &exonSequence);

    void setMetaInformation(const std::string &metaInformation);

    const std::string &getMetaInformation() const;

    const bool &getIfOrfShift() const;

    bool ifOverLap(const Transcript &transcript) const;

    bool ifOverLapIgnorStrand(const Transcript &transcript) const;

    void setIfOrfShift(const bool &ifOrfShift);

    const std::string &getScore() const;

    void setScore(const std::string &score);

    const std::string &getLastColumnInformation() const;

    void setLastColumnInformation(const std::string &lastColumnInformation);

    const std::string &getType() const;

    void setType(const std::string &type);

    void cleanCdsVector();

    void cleanSaveRam() {
        _geneomeSequence.clear();
        _cdsSequence.clear();
    }

    const std::string &getSource() const;

    void setSource(const std::string source);

    bool operator==(const Transcript &transcript) const {
        if (_chromeSomeName.compare(transcript._chromeSomeName) != 0) {
            return false;
        } else {
            if (_start != transcript._start) {
                return false;
            } else if (_end != transcript._end) {
                return false;
            } else {
                return true;
            }
        }
    }

    bool operator!=(const Transcript &transcript) const {
        if (_chromeSomeName.compare(transcript._chromeSomeName) != 0) {
            return true;
        } else {
            if (_start != transcript._start) {
                return true;
            } else if (_end != transcript._end) {
                return true;
            } else {
                return false;
            }
        }
    }

    bool operator<(const Transcript &transcript) const {
        if (_chromeSomeName.compare(transcript._chromeSomeName) != 0) {
            return _chromeSomeName < transcript._chromeSomeName;
        } else {
            if (_start != transcript._start) {
                if (_start < transcript._start) {
                    return true;
                } else {
                    return false;
                }
            } else if (_end != transcript._end) {
                if (_end < transcript._end) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false; //not perfecrt, but it should work for std::set<Transcripts>
            }
        }
    }

    bool operator>(const Transcript &transcript) const {

        if (_chromeSomeName.compare(transcript._chromeSomeName) != 0) {
            return _chromeSomeName > transcript._chromeSomeName;
        } else {
            if (_start != transcript._start) {
                if (_start > transcript._start) {
                    return true;
                } else {
                    return false;
                }
            } else if (_end != transcript._end) {
                if (_end > transcript._end) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return false; //not perfecrt, but it should work for std::set<Transcripts>
            }
        }
    }
};
