/*
 * =====================================================================================
 *
 *       Filename:  Gene.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/05/2017 15:47:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************


The gene. i.e. gff file

**************************************************************************/


#ifndef ANNOTATIONLIFTOVER_GENE_H
#define ANNOTATIONLIFTOVER_GENE_H

#include "STRAND.h"
#include "Transcript.h"
#include <string>
#include <vector>

class Gene{
    private:
        std::string _name;
        STRAND _strand;
        int _start;
        int _end;
        std::string _chromeSomeName;
        std::vector<std::string> _transcriptVector;
        std::vector<Transcript> _transcripts;
        void updateStartEnd ( const Transcript& transcript);
        std::string _source;
        std::string _score;
        std::string _codonFrame;
        std::string _lastColumnInformation;
    public:
        Gene(const std::string& name, const STRAND& strand);
        Gene(const std::string& name, const std::string & chromeSomeName, const STRAND& strand, const int & start, const int & end);
        Gene(const std::string& name, const std::string & chromeSomeName, const STRAND& strand);
        Gene();
        const std::string & getName();
        const STRAND & getStrand();
        const int& getStart() const;
        void setStart( const int & start);
        const int&  getEnd() const;
        void setEnd( const int & end);
        const std::string& getChromeSomeName() const;
        std::vector<std::string>& getTranscriptVector();
        std::vector<Transcript>& getTranscripts();
        //bool checkOverLapAndAddTranscript(Transcript& transcript);
        void addTranscript(const Transcript& transcript);
        std::string& getSource();
        void setSource(const std::string& source);
        const std::string &getScore() const;
        void setScore(const std::string &score);
        const std::string &getCodonFrame() const;
        void setCodonFrame(const std::string &codonFrame);
        const std::string &getLastColumnInformation() const;
        void setLastColumnInformation(const std::string &_lastColumnInformation);

};
struct compare_gene {
    inline bool operator() ( Gene& gene1, Gene& gene2 )
    {
        std::string chr1 = gene1.getChromeSomeName();
        std::string chr2 = gene2.getChromeSomeName();
        if( gene1.getChromeSomeName().compare(gene2.getChromeSomeName()) ==0 ){
            return gene1.getStart() < gene2.getStart();
        }else{
            return (chr1 < chr1);
        }
    }
};

#endif //ANNOTATIONLIFTOVER_GENE_H
