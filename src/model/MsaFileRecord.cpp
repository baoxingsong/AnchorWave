// =====================================================================================
//
//       Filename:  MsaFileRecord.cpp
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




 ************************************************************************/


#include "MsaFileRecord.h"


MsaFileRecord::MsaFileRecord(const int & _start, const int & _end){
    this->start = _start;
    this->end = _end;
}
const int & MsaFileRecord::getStart(){
    return this->start;
}
const int & MsaFileRecord::getEnd(){
    return this->end;
}
void MsaFileRecord::setStart( const int & _start ){
    this->start = _start;
}
void MsaFileRecord::setEnd( const int & _end ){
    this->end = _end;
}
void MsaFileRecord::addMsaSingleRecord(std::string & lineName, MsaSingleRecord & msaSingleRecord){
    this->msaSingleRecordRecords[lineName] = msaSingleRecord;
}
std::map<std::string, MsaSingleRecord>& MsaFileRecord::getMsaSingleRecordRecords(){
    return this->msaSingleRecordRecords;
}


void msaFileRead(MsaFileRecord & msaFileRecord, std::string & fileLocation, std::map<std::string, std::string>& sdiFiles){
    std::set<std::string> accessionNames;
    for( std::map<std::string, std::string>::iterator itName=sdiFiles.begin(); itName!=sdiFiles.end(); ++itName ){
        accessionNames.insert(itName->first);
    }
    msaFileRead(msaFileRecord, fileLocation, accessionNames);
}
void msaFileRead(MsaFileRecord & msaFileRecord, std::string & fileLocation, std::set<std::string>& accessionNames){
    std::ifstream infile(fileLocation);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << fileLocation << std::endl;
        exit (1);
    }
    std::regex reg1("^>(\\S+)");
    std::string name="";
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        std::smatch match1;
        regex_search(line, match1, reg1);
        if( match1.empty()   ){
            sequencestream << line;
        }else{
//            std::cout << line << std::endl;
            if( name.size()>1 ){
//                std::string sequence = sequencestream.str();
                std::regex reg2("^(\\S+?)_(\\d+)_(\\d+)");
                std::smatch match2;
                regex_search(name, match2, reg2);
                if( ! match2.empty()){
                    std::string accessionName = match2[1];
                    if(accessionNames.find(accessionName)!=accessionNames.end() || accessionName.compare("ref")==0 ) {
                        std::string sequence = sequencestream.str();
                        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                        MsaSingleRecord MsaSingleRecord(std::stoi(match2[2]), std::stoi(match2[3]), match2[1], sequence);
                        msaFileRecord.addMsaSingleRecord(accessionName, MsaSingleRecord);
                    }
                }
            }
            name=match1[1];
            sequencestream.str(std::string());
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        std::regex reg2("^(\\S+?)_(\\d+)_(\\d+)");
        std::smatch match2;
        regex_search(name, match2, reg2);
        if( ! match2.empty()){
//            std::cout << name << std::endl;
            std::string accessionName = match2[1];
            if(accessionNames.find(accessionName)!=accessionNames.end() || accessionName.compare("ref")==0) {
//                std::cout << accessionName << std::endl;
                std::string sequence = sequencestream.str();
                transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
                MsaSingleRecord msaSingleRecord(std::stoi(match2[2]), std::stoi(match2[3]), match2[1], sequence);
                msaFileRecord.addMsaSingleRecord(accessionName, msaSingleRecord);
            }
        }
    }
}