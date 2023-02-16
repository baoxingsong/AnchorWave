/*
 * =====================================================================================
 *
 *       Filename:  GenomeBasicFeature.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2017 15:47:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *        Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
/*************************************************************************

The class if designed for genetic elements CDS, i.e. from gff file


 ************************************************************************/

#pragma once

#include <string>

class GenomeBasicFeature {
private:
    int _start;
    int _end;
    std::string _score;
    std::string _codonFrame;
    std::string _lastColumnInformation;
    std::string type;
public:
    const std::string &getType() const;

    void setType(const std::string &type);

public:
    GenomeBasicFeature(const int &start, const int &end);

    GenomeBasicFeature(const int &start, const int &end, const std::string &score, const std::string &codonFrame, const std::string &lastColumnInformation);

    GenomeBasicFeature();

    const int &getStart() const;

    void setStart(const int &start);

    const int &getEnd() const;

    void setEnd(const int &end);

    const std::string &getScore() const;

    void setScore(const std::string &score);

    const std::string &getCodonFrame() const;

    void setCodonFrame(const std::string &codonFrame);

    const std::string &getLastColumnInformation() const;

    void setLastColumnInformation(const std::string &lastColumnInformation);

    bool operator<(const GenomeBasicFeature &cds) const {
        if (_start < cds._start) {
            return true;
        } else if (_start == cds._start && _end < cds._end) {
            return true;
        }
        return false;
    }

    bool operator>(const GenomeBasicFeature &cds) const {
        if (_start > cds._start) {
            return true;
        } else if (_start == cds._start && _end > cds._end) {
            return true;
        }
        return false;
    }

    bool operator==(const GenomeBasicFeature &cds) const {
        if (_start == cds._start && _end == cds._start) {
            return true;
        }
        return false;
    }

    bool operator!=(const GenomeBasicFeature &cds) const {
        if (_start == cds._start && _end == cds._start) {
            return false;
        }
        return true;
    }
};
