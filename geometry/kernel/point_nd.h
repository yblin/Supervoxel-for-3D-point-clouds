//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_POINT_3D_H_
#define GEOMETRY_KERNEL_POINT_3D_H_

#include <cassert>
#include <ostream>

#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/math/vector.h"

namespace cl {

/// 3D Point Class.
template<typename T>
class Point3D {
public:
    typedef T value_type;

    Point3D()
        : x(0), y(0), z(0) {}

    Point3D(const T& x1, const T& y1, const T& z1)
        : x(x1), y(y1), z(z1) {}

    bool operator ==(const Point3D& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator !=(const Point3D& rhs) const {
        return !(*this == rhs);
    }

    bool operator < (const Point3D& rhs) const {
        return x == rhs.x ? (y == rhs.y ? z < rhs.z : y < rhs.y) : x < rhs.x;
    }

    bool operator <=(const Point3D& rhs) const {
        return !(rhs < *this);
    }

    bool operator > (const Point3D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Point3D& rhs) const {
        return !(*this < rhs);
    }

    const Point3D& operator +=(const Vector3D<T>& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    const Point3D& operator -=(const Vector3D<T>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    /**
     * Return the i-th component value of point.
     */
    const T& operator[] (int i) const {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the reference value of the i-th component of point.
     */
    T& operator[] (int i) {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the i-th component value of point.
     */
    const T& at(int i) const {
        assert(0 <= i && i < 3);

        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the reference value of the i-th component of point.
     */
    T& at(int i) {
        assert(0 <= i && i < 3);

        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * @return the bounding box of this point.
     */
    const Box3D<T> bounding_box() const {
        return Box3D<T>(x, x, y, y, z, z);
    }

    /**
     * @return the dimension.
     */
    int size() const {
        return 3;
    }

    /**
     * Convert the point to vector.
     */
    Vector3D<T> ToVector() const {
        return Vector3D<T>(x, y, z);
    }

    friend Point3D operator +(const Point3D& lhs, const Vector3D<T>& rhs) {
        return Point3D(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    friend Point3D operator -(const Point3D<T>& lhs, const Vector3D<T>& rhs) {
        return Point3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    friend Vector3D<T> operator -(const Point3D& lhs, const Point3D& rhs) {
        return Vector3D<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    friend std::ostream& operator <<(std::ostream& os, const PointND& rhs) {
        os << rhs.ToVector();
        return os;
    }

    T x; /// X coordinate.
    T y; /// Y coordinate.
    T z; /// Z coordinate.
};

typedef Point3D<int>    IPoint3D;
typedef Point3D<double> RPoint3D;

} // namespace cl

#endif // GEOMETRY_KERNEL_POINT_3D_H_
