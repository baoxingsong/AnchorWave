/*
 * =====================================================================================
 *
 *       Filename:  FirstLastList.h
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

This is a link list data structure designed for storage and operation of genotype variant records.
Since the original variant records are from base-by-base comperation.
There are some variants should be merged together.
The merging process has a lot of variant records deletion and insertions.
The link list data structure could deal with records deletion and insertion very efficiently.

**************************************************************************/


#ifndef ANNOTATIONLIFTOVER_FIRSTLASTLIST_H
#define ANNOTATIONLIFTOVER_FIRSTLASTLIST_H

#include "variant.h"

class Data{
private:
    Variant mapSingleRecord;
    Data *next;
    Data *prev;
public:
    Data(Variant& mapSingleRecord);
    //        Data();
    const Variant & getMapSingleRecord() const;
    Data* getNext() const;
    Data* getPrev() const;
    void setNext( Data* data );
    void setPrev( Data* data);
    void setMapSingleRecord( Variant& mapSingleRecord);
};


class FirstLastList {
    private:
        Data *first = NULL;
        Data *last = NULL;
    public:
    //        void insertFirst(Data* data);
        void insertLast(Data* data);
        void setLast(Data* data);
        void deleteFirst();
        void deleteLast();
        Data* getFirst() const;
        Data* getLast() const;
};


#endif //ANNOTATIONLIFTOVER_FIRSTLASTLIST_H
