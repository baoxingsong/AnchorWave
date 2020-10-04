// =====================================================================================
//
//       Filename:  Variant.cpp
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





*************************************************************************/


#include "variant.h"


//variant record class begin

Variant::Variant(){
    changingLength=0;
    chromosome="";
    position=0;
    reference="";
    alternative="";
}
Variant::Variant(const std::string& _chromosome, const int32_t & _position, const std::string& _reference,
                 const std::string& _alternative){
    int32_t alternativeSize = _alternative.size();
    if( alternativeSize == 0 ){
        alternative="-";
    }else{
        alternative=_alternative;
    }
    if( _alternative.compare("-") ==0  ){
        alternativeSize = 0;
    }else{
        std::size_t found = _alternative.find("-");
        if (found!=std::string::npos){
            std::cerr << "ALTERNATIVE the record contains strange - at: " << _chromosome << "\t" << _position << "\t" << _alternative << std::endl;
            exit(1);
        }
    }
    int32_t referenceSize = _reference.size();
    if( referenceSize == 0 ){
        reference="-";
    }else{
        reference=_reference;
    }
    if( _reference.compare("-") ==0  ){
        referenceSize=0;
    }else {
        std::size_t found = _reference.find("-");
        if (found!=std::string::npos){
            std::cerr << "REFERENCE the record contains strange - at: " << _chromosome << "\t" << _position << "\t" << _reference << std::endl;
            exit(1);
        }
    }
    changingLength = alternativeSize - referenceSize;
    chromosome=_chromosome;
    position=_position;
}
const std::string & Variant::getChromosome() const{
    return this->chromosome;
}
const int32_t & Variant::getPosition() const {
    return this->position;
}
void Variant::setLastTotalChanged( const int32_t & _lastTotalChanged){
    this->lastTotalChanged=_lastTotalChanged;
}
const int32_t & Variant::getLastTotalChanged() const{
    return this->lastTotalChanged;
}
void Variant::setChangedPoint ( const int32_t & _changedPoint){
    this->changedPoint=_changedPoint;
}
const int32_t & Variant::getChangedPoint ( ) const{
    return this->changedPoint;
}
const int32_t & Variant::getChanginglength() const {
    return this->changingLength;
}
const std::string & Variant::getReference() const {
    return this->reference;
}
const std::string & Variant::getAlternative() const {
    return this->alternative;
}
Variant::Variant (const Variant& variant){
    this->chromosome=variant.chromosome;
    this->position = variant.position;
    this->changingLength = variant.changingLength;
    this->reference = variant.reference;
    this->alternative = variant.alternative;

}
bool Variant::overlap(const int32_t & _position) const{
    if( this->changingLength>=0 ){
        if( this->position == _position ){
            return true;
        }
    }else{
        if( (this->position-this->changingLength)>=_position && this->position<=_position  ){
            return true;
        }
    }
    return false;
}

std::ostream &printSdi(std::ostream& out, const Variant& variant){
    out << variant.getChromosome() << "\t" << variant.getPosition() << "\t" <<
        variant.getChanginglength() << "\t" << variant.getReference() << "\t" << variant.getAlternative();
    return out;
}
std::ostream &printSdiln(std::ostream& out, const Variant& variant){
    out << variant.getChromosome() << "\t" << variant.getPosition() << "\t" <<
        variant.getChanginglength() << "\t" << variant.getReference() << "\t" << variant.getAlternative() << std::endl;
    return out;
}

