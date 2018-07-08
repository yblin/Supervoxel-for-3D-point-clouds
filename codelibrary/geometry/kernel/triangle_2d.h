//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_TRIANGLE_2D_H_
#define GEOMETRY_KERNEL_TRIANGLE_2D_H_

#include <cmath>

#include "codelibrary/geometry/kernel/point_2d.h"

namespace cl {

template <typename T>
class Triangle2D {
public:
    Triangle2D(const Point2D<T>& p1, const Point2D<T>& p2,
               const Point2D<T>& p3)
        : p1_(p1), p2_(p2), p3_(p3) {}

    const Point2D<T>& p1() const { return p1_; }
    const Point2D<T>& p2() const { return p2_; }
    const Point2D<T>& p3() const { return p3_; }

    /**
     * @return the area of the triangle.
     */
    double area() const {
        double t1 = static_cast<double>(p2_.x) * p3_.y -
                    static_cast<double>(p3_.x) * p2_.y;
        double t2 = static_cast<double>(p3_.x) * p1_.y -
                    static_cast<double>(p1_.x) * p3_.y;
        double t3 = static_cast<double>(p1_.x) * p2_.y -
                    static_cast<double>(p2_.x) * p1_.y;
        return 0.5 * std::sqrt(t1 * t1 + t2 * t2 + t3 * t3);
    }

private:
    Point2D<T> p1_, p2_, p3_;
};

typedef Triangle2D<int> ITriangle2D;
typedef Triangle2D<double> RTriangle2D;

} // namespace cl

#endif // GEOMETRY_KERNEL_TRIANGLE_2D_H_
