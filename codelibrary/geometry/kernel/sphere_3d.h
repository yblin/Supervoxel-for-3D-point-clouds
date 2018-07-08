//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_SPHERE_3D_H_
#define GEOMETRY_KERNEL_SPHERE_3D_H_

#include <cassert>
#include <cmath>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"

namespace cl {

/// 3D Sphere.
template <typename T>
class Sphere3D {
public:
    typedef T value_type;

    Sphere3D()
        : radius_(0) {}

    Sphere3D(const Point3D<T>& center, const T& radius)
        : center_(center), radius_(radius) {
        assert(radius >= 0);

        bounding_box_ = Box3D<T>(center_.x - radius_, center_.x + radius_,
                                 center_.y - radius_, center_.y + radius_,
                                 center_.z - radius_, center_.z + radius_);
    }

    Sphere3D(const T& x, const T& y, const T& z, const T& radius)
        : Sphere3D(Point3D<T>(x, y, z), radius) {}

    /**
     * The radius of sphere.
     */
    const T& radius() const {
        return radius_;
    }

    /**
     * The bounding box of sphere.
     */
    const Box3D<T>& bounding_box() const {
        return bounding_box_;
    }

    /**
     * The center of sphere.
     */
    const Point3D<T>& center() const {
        return center_;
    }

protected:
    Box3D<T> bounding_box_; // The bounding box of sphere.
    Point3D<T> center_;     // The center of sphere.
    T radius_;              // The radius of sphere.
};

typedef Sphere3D<int>    ISphere3D;
typedef Sphere3D<double> RSphere3D;

} // namespace cl

#endif // GEOMETRY_KERNEL_SPHERE_3D_H_
