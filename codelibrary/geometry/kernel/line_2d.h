//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_LINE_2D_H_
#define GEOMETRY_KERNEL_LINE_2D_H_

#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {

/// 2D Line.
/**
 * A 2D line can be uniquely defined by a point and a direction vector.
 */
template <typename T>
class Line2D {
public:
    typedef T value_type;

    Line2D() {}

    /**
     * Construct from a point on the line and direction of line.
     */
    Line2D(const Point2D<T>& p, const Vector2D<T>& direction)
        : point1_(p), point2_(p + direction), direction_(direction) {}

    /**
     * Construct from two distinct points.
     */
    Line2D(const Point2D<T>& p1, const Point2D<T>& p2)
        : point1_(p1), point2_(p2), direction_(p2 - p1) {}

    /**
     * Construct from segment.
     */
    explicit Line2D(const Segment2D<T>& segment)
        : Line2D(segment.lower_point(), segment.upper_point()) {}

    const Vector2D<T>& direction() const {
        return direction_;
    }

    const Point2D<T>& point1() const {
        return point1_;
    }

    const Point2D<T>& point2() const {
        return point2_;
    }

protected:
    Point2D<T> point1_, point2_; // Two different points on the line.
    Vector2D<T> direction_;      // The direction vector of the line.
};

typedef Line2D<int>    ILine2D;
typedef Line2D<double> RLine2D;

} // namespace cl

#endif // GEOMETRY_KERNEL_LINE_2D_H_
