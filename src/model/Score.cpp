//
// Created by Baoxing Song on 2019-04-04.
//

#include "Score.h"

Scorei::Scorei(const signed char &matchingScore, const signed char &mismatchingPenalty) {
    int i, j;
    m = new signed char *[5];
    for (i = 0; i < 5; ++i) {
        m[i] = new signed char[5];
        std::fill_n(m[i], 5, 0);
    }

    for (i = 0; i < 5; ++i) {
        if (i != 2) {
            for (j = 0; j < 5; ++j) {
                if (j != 2) {
                    if (i == j) {
                        m[i][j] = matchingScore;
                    } else {
                        m[i][j] = mismatchingPenalty;
                    }
                }
            }
        }
    }
}

Scorei::~Scorei() {
    for (int i = 0; i < 5; ++i) {
        delete[] m[i];
    }
    delete[] m;
}

const signed char &Scorei::getScore(const signed char &a, const signed char &b) const {
    return m[a][b];
}

signed char **Scorei::getScore() const {
    return m;
}
