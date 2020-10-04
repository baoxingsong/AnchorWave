//
// Created by baoxing on 10/10/17.
//

#include "coordinateLiftOver.h"

/**
 * if could not find the chrosomeName, the basement would be returned.
 * to accelerate, a binary search algorithm is used here
 * @param chromsomeName
 * @param basement the coordinate of reference genome
 * @param variantsMap
 * @return the coordinate of target genome
 */
int getChangedFromBasement(const std::string & chromosomeName, const int & basement, std::map<std::string, std::vector<Variant> >& variantsMap){
    if( variantsMap.find(chromosomeName)!=variantsMap.end() && variantsMap[chromosomeName].size()>0 ){
        int start = 0 ;
        int end = variantsMap[chromosomeName].size()-1;
        int lastStart = start;
        if(variantsMap[chromosomeName][start].getPosition() >= basement){
            return basement;
        }
        if(variantsMap[chromosomeName][end].getPosition() <= basement){
            return basement+variantsMap[chromosomeName][end].getChanginglength()+variantsMap[chromosomeName][end].getLastTotalChanged();//allChanged;
        }
        while(!((variantsMap[chromosomeName][start].getPosition() < basement) && (variantsMap[chromosomeName][start+1].getPosition() >= basement))){
            if((variantsMap[chromosomeName][start].getPosition() < basement)){
                lastStart = start;
                if(1 == (end - start)){
                    start = end;
                }else{
                    start = (start+end)/2;
                }
            }else{
                end = start;
                start = lastStart;
            }
        }
        if((variantsMap[chromosomeName][start].getChanginglength() < 0) &&  basement<=variantsMap[chromosomeName][start].getPosition()-variantsMap[chromosomeName][start].getChanginglength()){
            return variantsMap[chromosomeName][start].getPosition() + variantsMap[chromosomeName][start].getLastTotalChanged();//lastChangedPoint;
        }
        return basement+variantsMap[chromosomeName][start+1].getLastTotalChanged();//allChanged;
    }else{
        return basement;
    }
}

/**
 * if could not find the chrosomeName, the changed would be returned.
 * to accelerate, a binary search algorithm is used here
 * @param chromsomeName
 * @param changed the coordinate of target genome
 * @return the coordinate of reference genome
 */
int getBasementFromChanged(const std::string &  chromsomeName, const int & changed, std::map<std::string, std::vector<Variant> >& variantsMap){
    if(variantsMap.find(chromsomeName)!=variantsMap.end() && variantsMap[chromsomeName].size()>0){
        int start = 0 ;
        int end = variantsMap[chromsomeName].size()-1;
        int lastStart = start;

        if(variantsMap[chromsomeName][end].getChangedPoint() <= changed){
            Variant mEnd = variantsMap[chromsomeName][end];
            return changed - (mEnd.getChanginglength() + mEnd.getLastTotalChanged());
        }
        if(variantsMap[chromsomeName][start].getChangedPoint() >= changed){
            return changed;
        }
        while(!(variantsMap[chromsomeName][start].getChangedPoint()<changed && variantsMap[chromsomeName][start+1].getChangedPoint()>=changed)){
            if((variantsMap[chromsomeName][start].getChangedPoint() < changed)){
                lastStart = start;
                if(1 == (end - start)){
                    start = end;
                }else{
                    start = (start+end)/2;
                }
            }else{
                end = start;
                start = lastStart;
            }
        }
        //if( variantsMap[chromsomeName][start].getChanginglength() > 0 && changed <= (variantsMap[chromsomeName][start-1].getChangedPoint()+((variantsMap[chromsomeName][start]).getChanginglength())) ){
        //modified by Baoxing at 2018 SEP 03
        if( variantsMap[chromsomeName][start].getChanginglength() > 0 && changed <= (variantsMap[chromsomeName][start].getChangedPoint()+((variantsMap[chromsomeName][start]).getChanginglength())) ){
            return (variantsMap[chromsomeName][start]).getPosition();
        }
        return variantsMap[chromsomeName][start+1].getPosition() - (variantsMap[chromsomeName][start+1].getChangedPoint() - changed);
    }else{
        return changed;
    }
}
