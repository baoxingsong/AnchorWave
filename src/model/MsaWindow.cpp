//
// Created by baoxing on 10/10/17.
//

#include "MsaWindow.h"


// for MSA
Window::Window(const int & _start, const int & _end,  const bool & _preExtand, const bool & _postExtand){
    if (_start < _end) {
        this->start = _start;
        this->end = _end;
    } else {
        this->start = _end;
        this->end = _start;
    }
    transcript_name = "";
    this->preExtand = _preExtand;
    this->postExtand = _postExtand;
}
Window::Window(const int & _start, const int & _end, const std::string & _transcript, const REGION & _region,
               const bool & _preExtand, const bool & _postExtand){
    if (_start < _end) {
        this->start = _start;
        this->end = _end;
    } else {
        this->start = _end;
        this->end = _start;
    }
    this->transcript_name = _transcript;
    this->region = _region;
    this->preExtand = _preExtand;
    this->postExtand = _postExtand;
}
const int & Window::getStart() const{
    return this->start;
}
const int & Window::getEnd() const{
    return this->end;
}
const std::string & Window::get_transcript_name() const{
    return this->transcript_name;
}
const REGION & Window::get_region() const{
    return this->region;
}
const bool & Window::getPreExtand() const{
    return this->preExtand;
}
const bool & Window::getPostExtand() const{
    return this->postExtand;
}
