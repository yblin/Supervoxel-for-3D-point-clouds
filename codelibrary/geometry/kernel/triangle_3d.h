//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_TRIANGLE_3D_H_
#define GEOMETRY_KERNEL_TRIANGLE_3D_H_

#include <cmath>

#include "codelibrary/geometry/kernel/point_3d.h"

namespace cl {

template <typename T>
class Triangle3D {
public:
    Triangle3D(const Point3D<T>& p1, const Point3D<T>& p2,
               const Point3D<T>& p3)
        : p1_(p1), p2_(p2), p3_(p3) {}

    const Point3D<T>& p1() const { return p1_; }
    const Point3D<T>& p2() const { return p2_; }
    const Point3D<T>& p3() const { return p3_; }

    /**
     * Return the normal vector of the triangle. Use anti-clockwise rule.
     */
    RVector3D normal() const {
        RVector3D v1(p2_.x - p1_.x, p2_.y - p1_.y, p2_.z - p1_.z);
        RVector3D v2(p3_.x - p1_.x, p3_.y - p1_.y, p3_.z - p1_.z);
        return CrossProduct(v1, v2);
    }

    /**
     * @return the area of the triangle.
     */
    double area() const {
        RVector3D v1 = p2_ - p1_, v2 = p3_ - p1_;
        RVector3D v = CrossProduct(v1, v2);

        return 0.5 * v.norm();
    }

private:
    Point3D<T> p1_, p2_, p3_;
};

typedef Triangle3D<int> ITriangle3D;
typedef Triangle3D<double> RTriangle3D;

} // namespace cl

#endif // GEOMETRY_KERNEL_TRIANGLE_3D_H_
