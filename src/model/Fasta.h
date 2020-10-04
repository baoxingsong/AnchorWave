/*
 * =====================================================================================
 *
 *       Filename:  Fasta.h
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

A class for single fasta records.
With two variables:
 1) sequence name
 2) sequence


 ************************************************************************/


#ifndef ANNOTATIONLIFTOVER_FASTA_H
#define ANNOTATIONLIFTOVER_FASTA_H

#include <string>
#include <iostream>
#include "STRAND.h"


//Fasta class begin
class Fasta{
    private:
        std::string name;
        std::string sequence;
    public:
        Fasta();
        Fasta(const std::string& _name, const std::string& _sequence);//When you declare any other constructor, the compiler will not generate the default constructor for you

        void setName(const std::string& _name);
        void setSequence(const std::string& _sequence);
        Fasta(const Fasta& f);//copy constructure
        const std::string& getName() const;
        const std::string& getSequence()const ;
};
std::ostream &print(std::ostream&, const Fasta&);

#endif //ANNOTATIONLIFTOVER_FASTA_H
