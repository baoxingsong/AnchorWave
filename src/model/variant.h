// =====================================================================================
//
//       Filename:  Variant.h
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

 A genotype variant records.

 It could be from sdi or vcf files and could output to sdi file



*************************************************************************/


#ifndef ANNOTATIONLIFTOVER_VARIANT_H
#define ANNOTATIONLIFTOVER_VARIANT_H

#include <string>
#include <iostream>
#include <cassert>

// variant record class begin
class Variant{
    private:
        std::string chromosome;
        int32_t position;
        int32_t changingLength;
        std::string reference;
        std::string alternative;
        int32_t lastTotalChanged;
        int32_t changedPoint;
    public:
        Variant();
        Variant(const std::string& _chromosome, const int32_t& _position, const std::string& _reference, const std::string& _alternative);
        const std::string & getChromosome() const;
        const int32_t& getPosition() const;
        void setLastTotalChanged(const int32_t & _lastTotalChanged);
        const int32_t& getLastTotalChanged( ) const;

        void setChangedPoint ( const int32_t & _changedPoint);
        const int32_t & getChangedPoint ( ) const;

        const int32_t & getChanginglength() const;
        const std::string & getReference() const;
        const std::string & getAlternative() const;
        Variant(const Variant& variant);
        bool overlap(const int32_t & position) const;

        bool operator<( const Variant& variant ) const {
            if(this->position < variant.position ){
                return true;
            }else if(this->position == variant.position && this->changingLength > variant.changingLength){
                return true;
            }else{
                return false;
            }
        }
        bool operator>( const Variant& variant ) const {
            if(this->position > variant.position ){
                return true;
            }else if(this->position == variant.position && this->changingLength < variant.changingLength){
                return true;
            }else{
                return false;
            }
        }
        bool operator == ( const Variant& variant  ) const {
            if(this->position == variant.position && this->changingLength == variant.changingLength){
                return true;
            }else{
                return false;
            }
        }
        bool operator != ( const Variant& variant  ) const {
            if(this->position == variant.position && this->changingLength == variant.changingLength){
                return false;
            }else{
                return true;
            }
        }
};


std::ostream &printSdi(std::ostream&, const Variant&);
std::ostream &printSdiln(std::ostream&, const Variant&);
struct compare_sdi_record {
    inline bool operator() ( Variant& variant1, Variant& variant2 ) {
        std::string chr1 = variant1.getChromosome();
        std::string chr2 = variant2.getChromosome();
        if( chr1.compare(chr2) ==0 ){
            if( variant1.getPosition() == variant2.getPosition() ){
                return 0-(variant1.getChanginglength() - variant2.getChanginglength());
            }else{
                return variant1.getPosition() < variant2.getPosition();
            }
        }else{
            return (chr1 < chr2);
        }
    }
};

#endif //ANNOTATIONLIFTOVER_VARIANT_H
