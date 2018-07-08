//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: get the center of geometric objects.
//

#ifndef GEOMETRY_UTIL_CENTER_2D_H_
#define GEOMETRY_UTIL_CENTER_2D_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {
namespace geometry {

/**
 * Compute the center of the smallest circle passing through the given points.
 */
template <typename T>
RPoint2D Circumcenter(const Point2D<T>& p1, const Point2D<T>& p2,
                      const Point2D<T>& p3) {
    double x1 = p1.x, y1 = p1.y;
    double x2 = p2.x, y2 = p2.y;
    double x3 = p3.x, y3 = p3.y;

    double x0 = 0.5 * ((y3 - y1) * (y2 * y2 - y1 * y1) +
                       (y3 - y1) * (x2 * x2 - x1 * x1) -
                       (y1 - y2) * (y1 * y1 - y3 * y3) -
                       (y1 - y2) * (x1 * x1 - x3 * x3)) /
                       ((y1 - y2) * (x3 - x1) - (y3 - y1) * (x1 - x2));
    return RPoint2D(x0, (y3 * y3 - y1 * y1 - 2.0 * x0 * (x3 - x1) -
                         x1 * x1 + x3 * x3) / (2.0 * (y3 - y1)));
}

/**
 * Compute the centroid of two 2D points.
 */
template <typename T>
RPoint2D Centroid(const Point2D<T>& p1, const Point2D<T>& p2) {
    return RPoint2D(0.5 * p1.x + 0.5 * p2.x, 0.5 * p1.y + 0.5 * p2.y);
}

/**
 * Compute the centroid of given 2D line segment.
 */
template <typename T>
RPoint2D Centroid(const Segment2D<T>& line) {
    return Centroid(line.lower_point(), line.upper_point());
}

/**
 * Compute the centroid of given 2D box.
 */
template <typename T>
RPoint2D Centroid(const Box2D<T>& box) {
    return Centroid(RPoint2D(box.x_min(), box.y_min()),
                    RPoint2D(box.y_min(), box.y_max()));
}

/**
 * Compute the centroid of given 2D points [first, last).
 * The second template parameter is used to distinguish this function to
 * Centroid2D(Point, Point).
 */
template <typename Iterator,
          typename = typename std::enable_if<std::is_convertible<
                     typename std::iterator_traits<Iterator>::iterator_category,
                              std::input_iterator_tag>::value>::type>
RPoint2D Centroid2D(Iterator first, Iterator last) {
    auto n = std::distance(first, last);
    assert(n > 0);

    RPoint2D centroid(0.0, 0.0);
    for (Iterator p = first; p != last; ++p) {
        centroid.x += p->x;
        centroid.y += p->y;
    }
    centroid.x /= n;
    centroid.y /= n;

    return centroid;
}

/**
 * Compute weighted centroid of given 2D points [first, last) and corresponding
 * weights [weight_first, weight_last).
 */
template <typename Iterator>
RPoint2D Centroid2D(Iterator first, Iterator last,
                    const Array<double>& weights) {
    int n = CountElements(first, last);
    assert(n > 0);
    assert(n == weights.size());

    RPoint2D centroid(0.0, 0.0);
    double sum = 0.0;
    int i = 0;
    for (Iterator p = first; p != last; ++p, ++i) {
        double w = weights[i];
        assert(w >= 0.0);

        centroid.x += w * p->x;
        centroid.y += w * p->y;
        sum += w;
    }
    assert(sum != 0.0);

    centroid.x /= sum;
    centroid.y /= sum;

    return centroid;
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_UTIL_CENTER_2D_H_
