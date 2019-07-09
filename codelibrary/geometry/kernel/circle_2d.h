//
// Copyright 2012 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_CIRCLE_2D_H_
#define GEOMETRY_KERNEL_CIRCLE_2D_H_

#include <cassert>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif // _USE_MATH_DEFINES

#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/kernel/box_2d.h"

namespace cl {

/// 2D Circle.
template <typename T>
class Circle2D {
public:
    typedef T value_type;

    Circle2D()
        : radius_(0) {}

    Circle2D(const Point2D<T>& center, const T& radius)
        : center_(center), radius_(radius) {
        assert(radius_ >= 0);

        bounding_box_ = Box2D<T>(center_.x - radius_, center_.x + radius_,
                                 center_.y - radius_, center_.y + radius_);
    }

    Circle2D(const T& x, const T& y, const T& radius)
        : Circle2D(Point2D<T>(x, y), radius) {}

    double area() const {
        return 0.5 * M_PI * radius_ * radius_;
    }

    /**
     * The radius of circle.
     */
    const T& radius() const {
        return radius_;
    }

    /**
     * Set a new radius.
     */
    void set_radius(double radius) {
        radius_ = radius;
    }

    /**
     * The bounding box of circle.
     */
    const Box2D<T>& bounding_box() const {
        return bounding_box_;
    }

    /**
     * The center of circle.
     */
    const Point2D<T>& center() const {
        return center_;
    }

protected:
    Box2D<T> bounding_box_; // The bounding box of circle.
    Point2D<T> center_;     // The center of circle.
    T radius_;              // The radius of circle.
};

typedef Circle2D<int>    ICircle2D;
typedef Circle2D<double> RCircle2D;

} // namespace cl

#endif // GEOMETRY_KERNEL_CIRCLE_2D_H_
