//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_METRIC_JACCARD_H_
#define UTIL_METRIC_JACCARD_H_

#include <algorithm>
#include <cassert>
#include <cmath>

#include "codelibrary/base/macros.h"
#include "codelibrary/util/set/dynamic_bitset.h"

namespace cl {
namespace metric {

/**
 * If x = {x_1, x_2, ..., x_n} and y = {y_1, y_2, ..., y_n} are two vectors with
 * all real x_i, y_i >= 0, then their Jacard similarity coefficient is defined
 * as:
 *             Sum_i min(x_i, y_i)
 *   J(x, y) = -------------------,
 *             Sum_j max(x_i, y_i)
 *
 * and Jaccard distance:
 *
 *   D_J(x, y) = 1 - J(x, y).
 */
class Jaccard {
public:
    Jaccard() {}

    template <typename T>
    double operator() (const T& a, const T& b) const {
        int size1 = a.size();
        int size2 = b.size();
        assert(size1 == size2);

        double t1 = 0.0, t2 = 0.0;
        for (int i = 0; i < size1; ++i) {
            double m1 = std::min(a[i], b[i]);
            double m2 = std::max(a[i], b[i]);
            assert(m1 >= 0.0);

            t1 += m1;
            t2 += m2;
        }

        // If A and B are both empty, then define D_J(A, B) = 0.
        if (t2 == 0.0) return 0.0;

        return 1.0 - t1 / t2;
    }

    /**
     * Specified for Dynamic bitset.
     *
     *                                    |A n B|
     * Jaccard Distance: D_J(A, B) = 1 - ---------.
     *                                    |A U B|
     */
    double operator() (const DynamicBitset& a, const DynamicBitset& b) const {
        assert(a.size() == b.size());

        DynamicBitset t = a & b;

        int n_intersections = t.Count();
        int n_unions = a.Count() + b.Count() - n_intersections;

        // If A and B are both empty, then define D_J(A, B) = 0.
        if (n_unions == 0) return 0.0;

        return 1.0 - static_cast<double>(n_intersections) / n_unions;
    }

private:
    DISALLOW_COPY_AND_ASSIGN(Jaccard);
};

} // namespace metric
} // namespace cl

#endif // UTIL_METRIC_JACCARD_H_
