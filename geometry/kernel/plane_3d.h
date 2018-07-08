//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_PLANE_3D_H_
#define GEOMETRY_KERNEL_PLANE_3D_H_

#include "codelibrary/geometry/kernel/line_3d.h"

namespace cl {

/// 3D Plane class.
/**
 * A 3D plane can be uniquely defined by a point on plane and a normal vector.
 */
template <typename T>
class Plane3D {
public:
    typedef T value_type;

    Plane3D() {}

    Plane3D(const Point3D<T>& point, const Vector3D<T>& normal)
        : point_(point), normal_(normal) {}

    /**
     * Construct plane by three non collinear points.
     */
    Plane3D(const Point3D<T>& a, const Point3D<T>& b, const Point3D<T>& c) {
        assert(a != b && a != c && b != c);

        Vector3D<T> v1(b - a), v2(c - a);
        normal_ = CrossProduct(v1, v2);
        point_ = a;
    }

    const Vector3D<T>& normal() const { return normal_; }
    const Point3D<T>& point()   const { return point_;  }

private:
    Point3D<T> point_;   // A point on the plane.
    Vector3D<T> normal_; // The normal vector of plane.
};

typedef Plane3D<int>    IPlane3D;
typedef Plane3D<double> RPlane3D;

} // namespace cl

#endif // GEOMETRY_KERNEL_PLANE_3D_H_
