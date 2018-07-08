//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_METRIC_COSINE_H_
#define UTIL_METRIC_COSINE_H_

#include <cassert>
#include <cmath>

#include "codelibrary/base/macros.h"

namespace cl {
namespace metric {

/**
 * Cosine is a measure of similarity between two vectors of an inner product
 * space that measures the cosine of the angle between them.
 *
 * The range of cosine metric is [-1, 1].
 */
class Cosine {
public:
    Cosine() {}

    template <typename T>
    double operator() (const T& a, const T& b) const {
        int size1 = a.size();
        int size2 = b.size();
        assert(size1 == size2);

        double t0 = 0.0, t1 = 0.0, t2 = 0.0;
        for (int i = 0; i < size1; ++i) {
            t0 += static_cast<double>(a[i]) * b[i];
            t1 += static_cast<double>(a[i]) * a[i];
            t2 += static_cast<double>(b[i]) * b[i];
        }

        return t0 / std::sqrt(t1 * t2);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(Cosine);
};

} // namespace metric
} // namespace cl

#endif // UTIL_METRIC_COSINE_H_
