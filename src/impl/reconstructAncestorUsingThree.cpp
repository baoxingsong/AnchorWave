//
// Created by bs674 on 6/16/21.
//

#include "reconstructAncestorUsingThree.h"

void reconstructAncestorUsingThree(const std::string & ref1, const std::string & ref2, const std::string & outGroup, std::string & ancestor ){
    assert ( ref1.length() == ref2.length() );
    assert ( ref1.length() == outGroup.length() );
    ancestor = "";
    for( int32_t i =0; i<ref1.length(); i++ ){
        if( ref1[i] == ref2[i] ){
            ancestor += ref1[i];
        }else if (ref2[i] == outGroup[i]){
            ancestor += ref2[i];
        }else{
            ancestor += ref1[i];
        }
    }
}
