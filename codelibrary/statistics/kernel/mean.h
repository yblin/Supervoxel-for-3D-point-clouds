//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_KERNEL_MEAN_H_
#define STATISTICS_KERNEL_MEAN_H_

#include <algorithm>
#include <cassert>
#include <cmath>

namespace cl {

/**
 * Get the arithmetic mean in range [first, last).
 *
 * Arithmetric mean is the sum of the list divided by its length.
 */
template <typename Iterator>
double Mean(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n > 0);

    return std::accumulate(first, last, 0.0) / n;
}

/**
 * Get the geometric mean in range [first, last).
 *
 * The geometric mean is the n-th root of the product of the list.
 */
template <typename Iterator>
double GeometricMean(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n > 0);

    double res = 1.0;
    for (Iterator i = first; i != last; ++i) {
        res *= *i;
    }

    return std::pow(res, 1.0 / n);
}

/**
 * Get the harmonic mean in range [first, last).
 *
 * Given a list {x1, x2, ..., xn}, the harmonic mean is n divided by the sum of
 * the reciprocal of each item in the list:
 *
 *    n / Sum(1 / xi).
 */
template <typename Iterator>
double HarmonicMean(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n > 0);

    double res = 0.0;
    for (Iterator i = first; i != last; ++i) {
        res += 1.0 / (*i);
    }

    return static_cast<double>(n) / res;
}

/**
 * Get the Contraharmonic mean in range [first, last).
 *
 * Given a list {x1, x2, ..., xn}, its Contraharmonic mean is:
 *
 *   Sum(xi^2) / Sum(xi).
 */
template <typename Iterator>
double ContraharmonicMean(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n > 0);

    double sum1 = 0.0, sum2 = 0.0;
    for (Iterator i = first; i != last; ++i) {
        sum1 += (*i) * (*i);
        sum2 += (*i);
    }

    return sum1 / sum2;
}

/**
 * Get the mean square in range [first, last).
 *
 * Given a list {x1, x2, ..., xn}, its mean square is:
 *
 *   Sum(xi^2) / n.
 */
template <typename Iterator>
double MeanSquare(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n > 0);

    double value = 0.0;
    for (Iterator p = first; p != last; ++p) {
        value += (*p) * (*p);
    }
    return value / n;
}

/**
 * Get the root mean square in range [first, last).
 *
 * Given a list {x1, x2, ..., xn}, its root mean square is:
 *
 *   Sqrt(Sum(xi^2) / n).
 */
template <typename Iterator>
double RootMeanSquare(Iterator first, Iterator last) {
    return std::sqrt(MeanSquare(first, last));
}

} // namespace cl

#endif // STATISTICS_KERNEL_MEAN_H_
