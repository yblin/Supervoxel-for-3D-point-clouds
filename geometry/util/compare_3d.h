//
// Copyright 2012 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: 3d compare for geometry objects.
//

#ifndef GEOMETRY_UTIL_COMPARE_3D_H_
#define GEOMETRY_UTIL_COMPARE_3D_H_

#include "codelibrary/geometry/kernel/point_3d.h"

namespace cl {
namespace geometry {

/**
 * Compare two 3D points along the given direction.
 */
class PointDirectionCompare3D {
public:
    template <typename T>
    explicit PointDirectionCompare3D(const Vector3D<T>& direction)
        : direction_(direction.x, direction.y, direction.z) {}

    template <typename T>
    bool operator() (const Point3D<T>& p1, const Point3D<T>& p2) const {
        RVector3D v1(p1.x, p1.y, p1.z);
        RVector3D v2(p2.x, p2.y, p2.z);
        return v1 * direction_ < v2 * direction_;
    }

private:
    RVector3D direction_;
};

/**
 * Compare two 3D points according to the their Z coordinate values.
 */
class ZCompare3D {
public:
    template <typename T>
    bool operator() (const Point3D<T>& p1, const Point3D<T>& p2) const {
        return p1.z < p2.z;
    }
};

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_UTIL_COMPARE_3D_H_
