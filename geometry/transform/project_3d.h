//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_TRANSFORM_PROJECT_3D_H_
#define GEOMETRY_TRANSFORM_PROJECT_3D_H_

#include <cmath>

#include "codelibrary/geometry/kernel/line_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/kernel/sphere_3d.h"

namespace cl {
namespace geometry {

/**
 * @return the orthogonal projection of a point 'p' on line.
 */
template <typename T1, typename T2>
RPoint3D Project(const Point3D<T1>& p, const Line3D<T2>& line) {
    const Point3D<T2>& q = line.point1();
    const Vector3D<T2>& v = line.direction();

    RVector3D v0(v.x, v.y, v.z);
    double t = v0 * v0;
    assert(t != 0.0 && "The input line is invalid.");

    RVector3D v1(static_cast<double>(p.x) - q.x,
                 static_cast<double>(p.y) - q.y,
                 static_cast<double>(p.z) - q.z);

    double b = (v0 * v1) / t;
    return RPoint3D(b * v0.x + q.x, b * v0.y + q.y, b * v0.z + q.z);
}

/**
 * @return the orthogonal projection of a point 'p' on plane.
 */
template <typename T1, typename T2>
RPoint3D Project(const Point3D<T1>& p, const Plane3D<T2>& plane) {
    const Vector3D<T2>& normal = plane.normal();
    RVector3D direction(normal.x, normal.y, normal.z);
    RLine3D line(RPoint3D(p.x, p.y, p.z), direction);
    return Project(plane.point(), line);
}

/**
 * @return the projection of a point 'p' on sphere.
 */
template <typename T1, typename T2>
RPoint3D Project(const Point3D<T1>& p, const Sphere3D<T2>& sphere) {
    double radius = sphere.radius();
    RPoint3D c(sphere.center().x, sphere.center().y, sphere.center().z);
    RPoint3D q(p.x, p.y, p.z);

    if (q == c) {
        // If p is the center of sphere, return the point that has with maximal
        // Z value.
        return c + RVector3D(0.0, 0.0, radius);
    }

    RVector3D v = q - c;
    v *= radius / v.norm();

    return c + v;
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_TRANSFORM_PROJECT_3D_H_
