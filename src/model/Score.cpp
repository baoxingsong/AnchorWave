//
// Created by Baoxing Song on 2019-04-04.
//

#include "Score.h"


//read score from files
Score::Score(const std::string & folder){
    std::ifstream infile(folder+"/matrixList");
    if( ! infile.good()){
        std::cerr << "error in opening score matrix list file " << folder+"/matrixList" << std::endl;
        exit (1);
    }
    std::string line="";
    while (std::getline(infile, line)){
        if( line.size()>0  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            int16_t index = std::stoi(line);
            std::ifstream infile2(folder+"/"+line);
            if( ! infile2.good()){
                std::cerr << "error in opening score matrix list file " << folder << "/" << line << std::endl;
                exit (1);
            }
            std::string line2="";
            int lineindex=0;
            int i;
            int32_t match = 0;
            int32_t misMatch = 0;
            while (std::getline(infile2, line2)){
                if( lineindex==0  ){
                    match = std::stoi(line2);
                }else if( lineindex==1  ){
                    misMatch = std::stoi(line2);
                }else if (lineindex==2){
                    this->openGapPenalty1[index]=std::stoi(line2);
                }else if (lineindex==3){
                    this->extendGapPenalty1[index]=std::stoi(line2);
                }else if (lineindex==4){
                    this->openGapPenalty2[index]=std::stoi(line2);
                }else if (lineindex==5){
                    this->extendGapPenalty2[index]=std::stoi(line2);
                }else if (lineindex==6){
                    this->zdrop[index]=std::stoi(line2);
                }

                ++lineindex;
            }

//            std::cout << "line 46 index: " << index << " match:" << match << " mis-match:" << misMatch << std::endl;
            int32_t ** score = new int32_t*[5];
            for (i = 0; i < 5; ++i) {
                score[i] = new int32_t[5];
                std::fill_n(score[i], 5, 0);
            }
            for ( i=0; i<5; ++i ){
                if  ( i!=2 ){
                    for ( int j=0; j<5; ++j ){
                        if( j!=2 ){
                            if ( i==j ){
                                score[i][j] = match;
                            }else{
                                score[i][j] = misMatch;
                            }
                        }
                    }
                }
            }
            this->m[index] = score;
//            std::cout << "line 65 index: " << index << " match:" << this->m[index][0][0] << " mis-match:" << this->m[index][0][1] << std::endl;
            infile2.close();
        }
    }
    infile.close();
}

Score::~Score() {
    for ( std::map<int16_t, int32_t**>::iterator it = m.begin(); it != m.end(); ++it ){
        for ( int i=0; i<5; ++i ){
            delete[] it->second[i];
        }
        delete[] it->second;
    }
}

int32_t ** Score::getM(const int16_t & category) {
    return m[category];
}
int32_t & Score::getOpenGapPenalty1(const int16_t & category) {
    return openGapPenalty1[category];
}
int32_t & Score::getExtendGapPenalty1(const int16_t & category){
    return extendGapPenalty1[category];
}
int32_t & Score::getOpenGapPenalty2(const int16_t & category) {
    return openGapPenalty2[category];
}
int32_t & Score::getExtendGapPenalty2(const int16_t & category) {
    return extendGapPenalty2[category];
}
int32_t & Score::getZdrop(const int16_t & category){
    return zdrop[category];
}


Scorei::Scorei( const int8_t & matchingScore, const int8_t & mismatchingPenalty ){
    int i, j;
    m = new int8_t*[5];
    for (i = 0; i < 5; ++i) {
        m[i] = new int8_t[5];
        std::fill_n(m[i], 5, 0);
    }
    for ( i=0; i<5; ++i ){
        if  ( i!=2 ){
            for ( j=0; j<5; ++j ){
                if( j!=2 ){
                    if ( i==j ){
                        m[i][j] = matchingScore;
                    }else{
                        m[i][j] = mismatchingPenalty;
                    }
                }
            }
        }
    }
}

Scorei::~Scorei(){
    for ( int i=0; i<5; ++i ){
        delete[] m[i];
    }
    delete[] m;
}
const int8_t & Scorei::getScore(const int8_t & a,const int8_t & b) const{
//    assert(a>=0);
//    assert(b>=0);
//    assert(a<5);
//    assert(b<5);
    return m[a][b];
}

int8_t ** Scorei::getScore() const{
    return m;
}
