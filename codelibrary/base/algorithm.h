//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_ALGORITHM_H_
#define BASE_ALGORITHM_H_

#include <algorithm>
#include <cassert>
#include <climits>
#include <functional>
#include <numeric>

#include "codelibrary/base/array.h"

namespace cl {

/**
 * Clamp the value to make sure it is in the range [low, high].
 *
 * Equivalent to C++17 clamp().
 */
template <typename T>
inline const T& Clamp(const T& value, const T& low, const T& high) {
    return value < low ? low : (value > high ? high : value);
}

/**
 * Count the elements of [first, last).
 */
template <typename Iterator>
int CountElements(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n >= 0);
    assert(n <= INT_MAX && "We only support INT_MAX elements at most.");

    return static_cast<int>(n);
}

/**
 * Get the sorted indices of [first, last)
 */
template <typename Iterator, typename Compare =
          std::less<typename std::iterator_traits<Iterator>::value_type> >
void IndexSort(Iterator first, Iterator last, Array<int>* indices) {
    assert(indices);

    auto n = std::distance(first, last);
    indices->resize(n);
    std::iota(indices->begin(), indices->end(), 0);

    Compare compare;
    std::sort(indices->begin(), indices->end(), [&](int a, int b) {
              return compare(first[a], first[b]);
    });
}

} // namespace cl

#endif // BASE_ALGORITHM_H_
