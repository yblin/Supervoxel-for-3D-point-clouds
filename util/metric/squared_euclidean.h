//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_METRIC_SQUARED_EUCLIDEAN_H_
#define UTIL_METRIC_SQUARED_EUCLIDEAN_H_

#include <cassert>
#include <cmath>

#include "codelibrary/base/macros.h"

namespace cl {
namespace metric {

/**
 * Squared Euclidean distance for two n-dimensional points.
 */
class SquaredEuclidean {
public:
    SquaredEuclidean() {}

    template <typename T>
    double operator() (const T& a, const T& b) const {
        int size1 = a.size();
        int size2 = b.size();
        assert(size1 == size2);

        double t = 0.0;
        for (int i = 0; i < size1; ++i) {
            t += static_cast<double>(a[i] - b[i]) * (a[i] - b[i]);
        }

        return t;
    }

private:
    DISALLOW_COPY_AND_ASSIGN(SquaredEuclidean);
};

} // namespace metric
} // namespace cl

#endif // UTIL_METRIC_SQUARED_EUCLIDEAN_H_
