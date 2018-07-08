//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_KERNEL_RANKING_H_
#define STATISTICS_KERNEL_RANKING_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"

namespace cl {

/**
 * Compute the competition ranking ("1224" ranking).
 *
 * In competition ranking, items that compare equal receive the same ranking
 * number, and then a gap is left in the ranking numbers.
 */
template <typename Iterator>
void CompetitionRanking(Iterator first, Iterator last, Array<int>* rank) {
    assert(rank);

    int n = CountElements(first, last);
    rank->resize(n);

    if (n == 0) return;

    Array<int> seq;
    IndexSort(first, last, &seq);

    (*rank)[seq[0]] = 0;
    int cur_rank = 0;
    for (int i = 1; i < n; ++i) {
        if (first[seq[i - 1]] < first[seq[i]]) {
            cur_rank = i;
        }
        (*rank)[seq[i]] = cur_rank;
    }
}

/**
 * Compute the dense ranking ("1223" ranking).
 *
 * In dense ranking, items that compare equal receive the same ranking number,
 * and the next item(s) receive the immediately following ranking number.
 */
template <typename Iterator>
void DenseRanking(Iterator first, Iterator last, Array<int>* rank) {
    assert(rank);

    int n = CountElements(first, last);
    rank->resize(n);

    if (n == 0) return;

    Array<int> seq;
    IndexSort(first, last, &seq);

    (*rank)[seq[0]] = 0;
    int cur_rank = 0;
    for (int i = 1; i < n; ++i) {
        if (first[seq[i - 1]] < first[seq[i]]) {
            ++cur_rank;
        }
        (*rank)[seq[i]] = cur_rank;
    }
}

/**
 * Compute the ordinal ranking ("1234" ranking).
 *
 * In ordinal ranking, all items receive distinct ordinal numbers, including
 * items that compare equal.
 */
template <typename Iterator>
void OrdinalRanking(Iterator first, Iterator last, Array<int>* rank) {
    assert(rank);

    int n = CountElements(first, last);
    rank->resize(n);

    if (n == 0) return;

    Array<int> seq;
    IndexSort(first, last, &seq);

    for (int i = 0; i < n; ++i) {
        (*rank)[seq[i]] = i;
    }
}

/**
 * Compute the ordinal ranking ("1 2.5 2.5 4" ranking).
 *
 * Items that compare equal receive the same ranking number, which is the mean
 * of what they would have under ordinal rankings.
 */
template <typename Iterator>
void FractionalRanking(Iterator first, Iterator last, Array<double>* rank) {
    assert(rank);

    int n = CountElements(first, last);
    rank->resize(n);

    if (n == 0) return;

    Array<int> seq;
    IndexSort(first, last, &seq);

    double sum_rank = 0.0;
    int begin = 0;
    for (int i = 1; i < n; ++i) {
        if (first[seq[i - 1]] < first[seq[i]]) {
            double mean = sum_rank / (i - begin);
            for (int j = begin; j < i; ++j) {
                (*rank)[seq[j]] = mean;
            }
            begin = i;
            sum_rank = i;
        } else {
            sum_rank += i;
        }
    }
    double mean = sum_rank / (n - begin);
    for (int i = begin; i < n; ++i) {
        (*rank)[seq[i]] = mean;
    }
}

} // namespace cl

#endif // STATISTICS_KERNEL_RANKING_H_
