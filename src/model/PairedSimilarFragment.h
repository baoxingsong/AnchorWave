//
// Created by Baoxing song on 2019-01-02.
//

#ifndef SONG_CNS_PAIREDSIMILARFRAGMENTS_H
#define SONG_CNS_PAIREDSIMILARFRAGMENTS_H


#include <cstdint>
#include <vector>
#include <iostream>
#include <sstream>

class PairedSimilarFragment {
    private:
        uint32_t start1;
        uint32_t end1;
        uint32_t start2;
        uint32_t end2;
        uint32_t score;
        std::vector<uint32_t>cigar;
        double pValue;
        double eValue;
        std::string alignment1;
        std::string alignment2;
public:
        const std::string &getAlignment1() const;
        void setAlignment1(const std::string &alignment1);
        const std::string &getAlignment2() const;
        void setAlignment2(const std::string &alignment2);

        PairedSimilarFragment();
        PairedSimilarFragment(uint32_t start1, uint32_t end1, uint32_t start2, uint32_t end2, uint32_t score,
                          const std::vector<uint32_t> &cigar, double & pValue, double & eValue);

        double getPValue() const;

        void setPValue(double pValue);

        double getEValue() const;

        void setEValue(double eValue);

        uint32_t getStart1() const;

        void setStart1(uint32_t start1);

        uint32_t getEnd1() const;

        void setEnd1(uint32_t end1);

        uint32_t getStart2() const;

        void setStart2(uint32_t start2);

        uint32_t getEnd2() const;

        void setEnd2(uint32_t end2);

        uint32_t getScore() const;
        uint32_t getLength() const;

        void setScore(uint32_t score);

        const std::vector<uint32_t> &getCigar() const;

        void setCigar(const std::vector<uint32_t> &cigar);
};

class PairedSimilarFragment2 {
    private:
        std::string species;
        std::string queryChr;
        uint32_t start1;
        uint32_t end1;
        uint32_t start2;
        uint32_t end2;
        std::string alignment1;
        std::string alignment2;
        int8_t strand; //1 positive, 0 negative

    public:
        PairedSimilarFragment2(const std::string &species, const std::string &queryChr, uint32_t start1, uint32_t end1,
                               uint32_t start2, uint32_t end2, const std::string &alignment1, const std::string &alignment2,
                               int8_t strand);
        const std::string &getSpecies() const;
        void setSpecies(const std::string &species);
        const std::string &getQueryChr() const;
        void setQueryChr(const std::string &queryChr);
        uint32_t getStart1() const;
        void setStart1(uint32_t start1);
        uint32_t getEnd1() const;
        void setEnd1(uint32_t end1);
        uint32_t getStart2() const;
        void setStart2(uint32_t start2);
        uint32_t getEnd2() const;
        void setEnd2(uint32_t end2);
        const std::string &getAlignment1() const;
        void setAlignment1(const std::string &alignment1);
        const std::string &getAlignment2() const;
        void setAlignment2(const std::string &alignment2);
        int8_t getStrand() const;
        void setStrand(int8_t strand);
};

#endif //SONG_CNS_PAIREDSIMILARFRAGMENTS_H
