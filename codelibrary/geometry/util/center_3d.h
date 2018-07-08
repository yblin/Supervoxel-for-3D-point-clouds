//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: get the center of geometric objects.
//

#ifndef GEOMETRY_UTIL_CENTER_3D_H_
#define GEOMETRY_UTIL_CENTER_3D_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/kernel/segment_3d.h"

namespace cl {
namespace geometry {

/**
 * Compute the centroid of two 3D points.
 */
template <typename T>
RPoint3D Centroid(const Point3D<T>& p1, const Point3D<T>& p2) {
    return RPoint3D(0.5 * p1.x + 0.5 * p2.x, 0.5 * p1.y + 0.5 * p2.y,
                    0.5 * p1.z + 0.5 * p2.z);
}

/**
 * Compute the centroid of given 3D bounding box.
 */
template <typename T>
RPoint3D Centroid(const Box3D<T>& box) {
    return Centroid(RPoint3D(box.x_min(), box.y_min(), box.z_min()),
                    RPoint3D(box.y_max(), box.y_max(), box.z_max()));
}

/**
 * Compute the centroid of given 3D line segment.
 */
template <typename T>
RPoint3D Centroid(const Segment3D<T>& line) {
    return Centroid(line.lower_point(), line.upper_point());
}

/**
 * Compute the centroid of given 3D points [first, last).
 * The second template parameter is used to distinguish this function to
 * Centroid3D(Point, Point).
 */
template <typename Iterator,
          typename = typename std::enable_if<std::is_convertible<
                     typename std::iterator_traits<Iterator>::iterator_category,
                              std::input_iterator_tag>::value>::type>
RPoint3D Centroid3D(Iterator first, Iterator last) {
    int n = CountElements(first, last);
    assert(n > 0);

    RPoint3D centroid(0.0, 0.0, 0.0);
    for (Iterator p = first; p != last; ++p) {
        centroid.x += p->x;
        centroid.y += p->y;
        centroid.z += p->z;
    }
    double t = 1.0 / n;
    centroid.x *= t;
    centroid.y *= t;
    centroid.z *= t;

    return centroid;
}

/**
 * Compute weighted centroid of given 3D points [first, last) and corresponding
 * weights.
 */
template <typename Iterator>
RPoint3D Centroid3D(Iterator first, Iterator last,
                    const Array<double>& weights) {
    int n = CountElements(first, last);
    assert(n > 0);
    assert(n == weights.size());

    double sum = 0.0;
    RPoint3D centroid(0.0, 0.0, 0.0);
    int i = 0;
    for (Iterator p = first; p != last; ++p, ++i) {
        double w = weights[i];
        assert(w >= 0.0);

        centroid.x += w * p->x;
        centroid.y += w * p->y;
        centroid.z += w * p->z;
        sum += w;
    }
    assert(sum != 0.0);

    sum = 1.0 / sum;
    centroid.x *= sum;
    centroid.y *= sum;
    centroid.z *= sum;

    return centroid;
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_UTIL_CENTER_3D_H_
