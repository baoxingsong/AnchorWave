//
// Created by song on 8/4/18.
//

#include "AlignmentMatch.h"

const Range & AlignmentMatch::getQuery() const {
    return this->query;
}

void AlignmentMatch::setQuery(const Range &query) {
    this->query = query;
}

const Range & AlignmentMatch::getDatabase() const {
    return this->database;
}

void AlignmentMatch::setDatabase(const Range &database) {
    this->database = database;
}

AlignmentMatch::AlignmentMatch(const std::string & queryChr, const size_t & queryStart, const size_t & queryEnd, const STRAND & queryStrand,
                             const std::string & databaseChr, const size_t & databaseStart, const size_t & databaseEnd/*, const STRAND & databaseStrand*/ ){
    if( queryStart > queryEnd ){
        query.setStart(queryEnd);
        query.setEnd(queryStart);
    }else{
        query.setStart(queryStart);
        query.setEnd(queryEnd);
    }

    query.setChr(queryChr);
    query.setStrand(queryStrand);

    database.setChr(databaseChr);
//    database.setStrand(databaseStrand);

    if( databaseStart > databaseEnd ){
        database.setStart(databaseEnd);
        database.setEnd(databaseStart);
    }else{
        database.setStart(databaseStart);
        database.setEnd(databaseEnd);
    }
}
AlignmentMatch::AlignmentMatch(const std::string & queryChr, const size_t & queryStart, const size_t & queryEnd, const STRAND & queryStrand,
               const std::string & databaseChr, const size_t & databaseStart, const size_t & databaseEnd, const size_t & _windowSize){
    if( queryStart > queryEnd ){
        query.setStart(queryEnd);
        query.setEnd(queryStart);
    }else{
        query.setStart(queryStart);
        query.setEnd(queryEnd);
    }

    query.setChr(queryChr);
    query.setStrand(queryStrand);

    database.setChr(databaseChr);
//    database.setStrand(databaseStrand);

    if( databaseStart > databaseEnd ){
        database.setStart(databaseEnd);
        database.setEnd(databaseStart);
    }else{
        database.setStart(databaseStart);
        database.setEnd(databaseEnd);
    }
    this->windowSize=_windowSize;
}

const std::string & AlignmentMatch::getQueryChr(){
    return this->query.getChr();
}
void AlignmentMatch::setQueryChr(const std::string & queryChr){
    this->query.setChr(queryChr);
}
const size_t & AlignmentMatch::getQueryStart() const{
    return this->query.getStart();
}
void AlignmentMatch::setQueryStart(const size_t & queryStart){
    if( queryStart<1 ){
        this->query.setStart(1);
        return;
    }
    this->query.setStart(queryStart);
}
const size_t & AlignmentMatch::getQueryEnd() const{
    return query.getEnd();
}
void AlignmentMatch::setQueryEnd( const size_t & queryEnd ){
    this->query.setEnd(queryEnd);
}
const STRAND & AlignmentMatch::getQueryStrand() const{
    return this->query.getStrand();
}
void AlignmentMatch::setQueryStrand(const STRAND & queryStrand){
    this->query.setStrand(queryStrand);
}
const std::string & AlignmentMatch::getDatabaseChr() const{
    return this->database.getChr();
}
void AlignmentMatch::setDatabaseChr(const std::string & databaseChr){
    this->database.setChr(databaseChr);
}
const size_t & AlignmentMatch::getDatabaseStart() const{
    return this->database.getStart();
}
void AlignmentMatch::setDatabaseStart(const size_t & databaseStart){
    if( databaseStart<1 ){
        this->database.setStart(1);
        return;
    }
    this->database.setStart(databaseStart);
}
const size_t & AlignmentMatch::getDatabaseEnd() const{
    return this->database.getEnd();
}
void AlignmentMatch::setDatabaseEnd(const size_t & databaseEnd){
    this->database.setEnd(databaseEnd);
}

size_t AlignmentMatch::getWindowSize() const {
    return windowSize;
}

void AlignmentMatch::setWindowSize(size_t windowSize) {
    AlignmentMatch::windowSize = windowSize;
}

const int & AlignmentMatch::getQueryIndex() const {
    return this->query.getIndex();
}

void AlignmentMatch::setQueryIndex( const int & index  ) {
    this->query.setIndex(index);
}
const int & AlignmentMatch::getDatabaseIndex() const{
    return this->database.getIndex();
}
void AlignmentMatch::setDatabaseIndex( const int & index  ){
    this->database.setIndex(index);
}
const size_t & AlignmentMatch::getRefLength() const{
    return  this->database.getEnd() - this->database.getStart() + 1;
}

const std::string & AlignmentMatch::getDatabaseName() const{
    return this->database.getName();
}
void AlignmentMatch::setDatabaseName(const std::string &name){
    this->database.setName(name);
}
const std::string & AlignmentMatch::getQueryName() const{
    return this->query.getName();
}
void AlignmentMatch::setQueryName(const std::string &name){
    this->query.setName(name);
}

AlignmentMatch::AlignmentMatch(){

}
AlignmentMatch::AlignmentMatch( const AlignmentMatch & alignmentMatch){
    query = alignmentMatch.getQuery();
    database = alignmentMatch.getDatabase();
    windowSize=alignmentMatch.getWindowSize();
}

/*const STRAND & AlignmentMatch::getDatabaseStrand() const{
    return this->database.getStrand();
}
void AlignmentMatch::setDatabaseStrand(const STRAND & strand){
    this->database.setStrand(strand);
}
*/
