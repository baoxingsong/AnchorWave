/*
 * =====================================================================================
 *
 *       Filename:  FirstLastList.cpp
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




**************************************************************************/


#include "FirstLastList.h"



Data::Data(Variant& mapSingleRecord){
    this->mapSingleRecord = mapSingleRecord;
    this->prev = NULL;
    this->next = NULL;
}

const Variant& Data::getMapSingleRecord() const{
    return mapSingleRecord;
}
Data* Data::getNext() const{
    return this->next;
}
Data* Data::getPrev() const{
    return this->prev;
}
void Data::setNext( Data* data ){
    this->next = data;
}
void Data::setPrev( Data*  data){
    this->prev = data;
}
void Data::setMapSingleRecord( Variant& mapSingleRecord){
    this->mapSingleRecord=mapSingleRecord;
}


void FirstLastList::insertLast(Data* data){
    if(this->first == NULL){
        this->first = data;
    }else{
        this->last->setNext(data);
        data->setPrev(this->last);
    }
    this->last = data;
}
void FirstLastList::setLast(Data* data){
    if( data == NULL ){
        this->first = NULL;
        this->last = NULL;
    }else{
        data->setNext(NULL);
        this->last = data;
    }
}
void FirstLastList::deleteFirst(){
    if(this->first == NULL){ // this is a empty link
        return;// this->first->getMapSingleRecord();
    }

    if( this->first->getNext() == NULL) { // this link has only one element
        delete(this->first);
        this->last = NULL;
        this->first = NULL;
        return;
    }
    Data * temp = this->first;
    this->first->getNext()->setPrev( NULL );
    this->first = (this->first->getNext());
    delete(temp);
}

void FirstLastList::deleteLast(){
    if(this->first == NULL){
        return;
    }
    if( this->first == this->last ){
        delete(this->last);
        this->first = NULL;
        this->last = NULL;
        return;
    }
    Data * temp = this->last;
    this->last = (last->getPrev());
    this->last->setNext(NULL);
    delete(temp);
}

Data* FirstLastList::getFirst() const{
    return this->first;
}
Data* FirstLastList::getLast() const {
    return this->last;
}
