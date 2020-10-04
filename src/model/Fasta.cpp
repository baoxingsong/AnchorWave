/*
 * =====================================================================================
 *
 *       Filename:  Fasta.cpp
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




 ************************************************************************/


#include "Fasta.h"
//Fasta class begin
Fasta::Fasta(const std::string& _name, const std::string& _sequence){
    this->name=_name;
    this->sequence=_sequence;
}
Fasta::Fasta(){

}
void Fasta::setName(const std::string& _name){
    name=_name;
}
void Fasta::setSequence(const std::string& _sequence){
    sequence=_sequence;
}
Fasta::Fasta( const Fasta& f ){
    this->name=f.name;
    this->sequence=f.sequence;
}

const std::string& Fasta::getSequence() const{
    return sequence;
}

const std::string& Fasta::getName() const{
    return name;
}
std::ostream& print(std::ostream& out, const Fasta& f){
    out << ">" << f.getName() << std::endl << f.getSequence() << std::endl;
    return out;
}
