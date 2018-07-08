//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_KERNEL_MEDIAN_H_
#define STATISTICS_KERNEL_MEDIAN_H_

#include <algorithm>
#include <cassert>

namespace cl {

/**
 * Get the median in range [first, last).
 */
template <typename Iterator>
const typename std::iterator_traits<Iterator>::value_type
Median(Iterator first, Iterator last) {
    assert(first != last);

    typedef typename std::iterator_traits<Iterator>::value_type T;

    Array<T> values(first, last);
    std::nth_element(values.begin(), values.begin() + values.size() / 2,
                     values.end());
    return values[values.size() / 2];
}

} // namespace cl

#endif // STATISTICS_KERNEL_MEDIAN_H_
