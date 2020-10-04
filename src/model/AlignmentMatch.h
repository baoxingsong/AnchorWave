//
// Created by song on 8/4/18.
//
//
// the strand of database is always POSITIVE
//
//
#ifndef ANNOTATIONLIFTOVER_ALIGNMENTMATCH_H
#define ANNOTATIONLIFTOVER_ALIGNMENTMATCH_H

#include "./Range.h"
class AlignmentMatch {
private:
    Range query;
    Range database;
    size_t windowSize=0;
public:
    size_t getWindowSize() const;

    void setWindowSize(size_t windowSize);

    AlignmentMatch(const std::string & queryChr, const size_t & queryStart, const size_t & queryEnd, const STRAND & queryStrand,
            const std::string & databaseChr, const size_t & databaseStart, const size_t & databaseEnd/*, const STRAND & databaseStrand*/ );
    AlignmentMatch(const std::string & queryChr, const size_t & queryStart, const size_t & queryEnd, const STRAND & queryStrand,
                   const std::string & databaseChr, const size_t & databaseStart, const size_t & databaseEnd, const size_t & _windowSize);
    AlignmentMatch();
    AlignmentMatch( const AlignmentMatch & alignmentMatch);
    const Range &getQuery() const;
    void setQuery(const Range &query);
    const Range &getDatabase() const;
    void setDatabase(const Range &database);

    const std::string & getQueryChr();
    void setQueryChr(const std::string & queryChr);
    const size_t & getQueryStart() const;
    void setQueryStart(const size_t & queryStart);
    const size_t & getQueryEnd() const;
    void setQueryEnd( const size_t & queryEnd );
    const STRAND & getQueryStrand() const;
    void setQueryStrand(const STRAND & queryStrand);
    const std::string & getDatabaseChr() const;
    void setDatabaseChr(const std::string & databaseChr);
    const size_t & getDatabaseStart() const;
    void setDatabaseStart(const size_t & databaseStart);
    const size_t & getDatabaseEnd() const;
    void setDatabaseEnd(const size_t & databaseEnd);
    const size_t & getRefLength() const;
    const int & getQueryIndex() const;
    void setQueryIndex( const int & index  ) ;
    const int & getDatabaseIndex() const;
    void setDatabaseIndex( const int & index );

    const std::string &getDatabaseName() const;
    void setDatabaseName(const std::string &name);
    const std::string & getQueryName() const;
    void setQueryName(const std::string &name);
    /*
    const STRAND & getDatabaseStrand() const;
    void setDatabaseStrand(const STRAND & strand);*/


    bool operator<( const AlignmentMatch& alignmentMatch ) const{
        if( database.getStart() < alignmentMatch.getDatabaseStart()){
            return true;
        }else if(database.getStart() == alignmentMatch.getDatabaseStart() && query.getStart()<alignmentMatch.getQueryStart() ){
            return true;
        }
        return false;
    }
    bool operator>(const AlignmentMatch& alignmentMatch )const {
        if( database.getStart() > alignmentMatch.getDatabaseStart()){
            return true;
        }else if(database.getStart() == alignmentMatch.getDatabaseStart() && query.getStart()>alignmentMatch.getQueryStart() ){
            return true;
        }
        return false;
    }
    bool operator==(const AlignmentMatch& alignmentMatch ) const{
        if( database.getStart() == alignmentMatch.getDatabaseStart() &&  query.getStart() == alignmentMatch.getQueryStart()){
            return true;
        }
        return false;
    }
    bool operator!=(const AlignmentMatch& alignmentMatch ) const{
        if( database.getStart() == alignmentMatch.getDatabaseStart()  &&  query.getStart() == alignmentMatch.getQueryStart() ){
            return false;
        }
        return true;
    }
};

#endif //ANNOTATIONLIFTOVER_ALIGNMENTMATCH_H
