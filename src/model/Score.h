/*
 * =====================================================================================
 *
 *       Filename:  Score.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  04/04/2019 19:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************


Reading the score from configure and provide score query functions

 ************************************************************************/

#pragma once

#include <algorithm>

class Scorei {
private:
    signed char **m;
public:
    Scorei(const signed char &matchingScore, const signed char &mismatchingPenalty);

    ~Scorei();

    const signed char &getScore(const signed char &a, const signed char &b) const;

    signed char **getScore() const;
};
