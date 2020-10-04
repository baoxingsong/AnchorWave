// =====================================================================================
//
//       Filename:  MsaWindow.h
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

 The MSA based variants recalling pipeline, cut the genome sequence into windows firstly,
    and them perform MSA, and recall variants

The class is the window for genome sequence cutting.
 If we perform gene structure based MSA.
If the window includes a whole CDS record, then there is no need for overlap with neighboring window.
 That is what preExtend and postEntend used for.

*************************************************************************/


#ifndef ANNOTATIONLIFTOVER_MSAWINDOW_H
#define ANNOTATIONLIFTOVER_MSAWINDOW_H

#include "REGION.h"
#include <string>

// for MSA pre
class Window {
    private:
        int start;
        int end;
        std::string transcript_name;
        REGION region;
        bool preExtand;
        bool postExtand;
    public:
        Window(const int & _start, const int & _end,  const bool & _preExtand, const bool & _postExtand);
        Window(const int & _start, const int & _end, const std::string & _transcript, const REGION & _region, const bool & _preExtand, const bool & _postExtand);
        const int & getStart() const;
        const int & getEnd() const;
        const std::string & get_transcript_name() const;
        const REGION & get_region() const;
        const bool & getPreExtand() const;
        const bool & getPostExtand() const;
        bool operator<( const Window& window ) const {
            if( start < window.start ){
                return true;
            }else if (start==window.start && end<window.end ){
                return true;
            }else {
                return false;
            }
        }
        bool operator>( const Window& window ) const {
            if( start > window.start ){
                return true;
            }else if (start==window.start && end>window.end ){
                return true;
            }else {
                return false;
            }
        }
        bool operator==( const Window& window ) const {
            if (start==window.start && end==window.end ){
                return true;
            }else {
                return false;
            }
        }
};


#endif //ANNOTATIONLIFTOVER_MSAWINDOW_H
