//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_LINE_3D_H_
#define GEOMETRY_KERNEL_LINE_3D_H_

#include "codelibrary/geometry/kernel/segment_3d.h"
#include "codelibrary/math/vector.h"

namespace cl {

/// 3D Line.
/**
 * A 3D line can be uniquely defined by a point and a direction vector.
 */
template <typename T>
class Line3D {
public:
    typedef T value_type;

    /**
     * Create a default line that passes through the origin and towrds to the
     * positive Y-axis.
     */
    Line3D() {}

    /**
     * Construct from a point on the line and direction of line.
     */
    Line3D(const Point3D<T>& p, const Vector3D<T>& direction)
        : point1_(p), point2_(p + direction), direction_(direction) {}

    /**
     * Construct from two distinct points.
     */
    Line3D(const Point3D<T>& p1, const Point3D<T>& p2)
        : point1_(p1), point2_(p2), direction_(p2 - p1) {}

    /**
     * Construct from segment.
     */
    explicit Line3D(const Segment3D<T>& segment)
        : Line3D(segment.lower_point(), segment.direction()) {}

    const Vector3D<T>& direction() const {
        return direction_;
    }

    const Point3D<T>& point1() const {
        return point1_;
    }

    const Point3D<T>& point2() const {
        return point2_;
    }

protected:
    Point3D<T> point1_, point2_; // Two different points on the line.
    Vector3D<T> direction_; // The direction vector of line.
};

typedef Line3D<int>    ILine3D;
typedef Line3D<double> RLine3D;

} // namespace cl

#endif // GEOMETRY_KERNEL_LINE_3D_H_
