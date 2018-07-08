//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: Euclidean distance between 3D geometric kernel objects.
//

#ifndef GEOMETRY_UTIL_DISTANCE_3D_H_
#define GEOMETRY_UTIL_DISTANCE_3D_H_

#include <cmath>

#include "codelibrary/geometry/kernel/line_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/segment_3d.h"

namespace cl {
namespace geometry {

/**
 * Return the squared Euclidean distance between two 3D points.
 */
template <typename T>
double SquaredDistance(const Point3D<T>& p1, const Point3D<T>& p2) {
    double t1 = static_cast<double>(p1.x) - p2.x;
    double t2 = static_cast<double>(p1.y) - p2.y;
    double t3 = static_cast<double>(p1.z) - p2.z;
    return t1 * t1 + t2 * t2 + t3 * t3;
}

/**
 * Return the Euclidean distance between two 3D points.
 */
template <typename T>
double Distance(const Point3D<T>& p1, const Point3D<T>& p2) {
    return std::sqrt(SquaredDistance(p1, p2));
}

/**
 * Return the Euclidean distance from point to plane.
 */
template <typename T>
double Distance(const Point3D<T>& p, const Plane3D<T>& plane) {
    double a = plane.normal().x;
    double b = plane.normal().y;
    double c = plane.normal().z;
    double d = -a * plane.point().x - b * plane.point().y - c * plane.point().z;

    return std::fabs(a * p.x + b * p.y + c * p.z + d) /
           std::sqrt(a * a + b * b + c * c);
}
template <typename T>
double Distance(const Plane3D<T>& plane, const Point3D<T>& p) {
    return Distance(p, plane);
}

/**
 * Return the Euclidean distance from point to segment.
 */
template <typename T>
double Distance(const Point3D<T>& p, const Segment3D<T>& segment) {
    if (segment.lower_point() == segment.upper_point()) {
        return Distance(p, segment.lower_point());
    }

    Vector3D<T> v0 = segment.upper_point() - segment.lower_point();
    Vector3D<T> v1 = p - segment.lower_point();

    double t1 = v0 * v1;
    if (t1 < 0.0) {
        return Distance(p, segment.lower_point());
    }

    double t2 = v0 * v0;
    if (t2 <= t1) {
        return Distance(p, segment.upper_point());
    }

    double b = t1 / t2;
    RPoint3D pb(b * v0.x + segment.lower_point().x,
                b * v0.y + segment.lower_point().y,
                b * v0.z + segment.lower_point().z);

    return Distance(p, pb);
}
template <typename T>
double Distance(const Segment3D<T>& seg, const Point3D<T>& p) {
    return Distance(p, seg);
}

/**
 * Return the Euclidean distance from point to line.
 */
template <typename T>
double Distance(const Point3D<T>& p, const Line3D<T>& line) {
    const Vector3D<T>& v0 = line.direction();
    Vector3D<T> v1 = p - line.point1();

    double t1 = v0 * v1;
    double t2 = v0 * v0;
    double b = t1 / t2;
    RPoint3D pb(b * v0.x + line.point1().x, b * v0.y + line.point1().y,
                b * v0.z + line.point1().z);

    return Distance(p, pb);
}
template <typename T>
double Distance(const Line3D<T>& line, const Point3D<T>& p) {
    return Distance(p, line);
}

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_UTIL_DISTANCE_3D_H_
