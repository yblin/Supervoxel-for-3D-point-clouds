//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_POINT_2D_H_
#define GEOMETRY_KERNEL_POINT_2D_H_

#include <cassert>
#include <ostream>

#include "codelibrary/geometry/kernel/box_2d.h"
#include "codelibrary/math/vector.h"

namespace cl {

/// 2D Point Class.
template<typename T>
class Point2D {
public:
    typedef T value_type;

    Point2D()
        : x(0), y(0) {}

    Point2D(const T& x1, const T& y1)
        : x(x1), y(y1) {}

    bool operator ==(const Point2D& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator !=(const Point2D& rhs) const {
        return !(*this == rhs);
    }

    bool operator <(const Point2D& rhs) const {
        return x < rhs.x || (x == rhs.x && y < rhs.y);
    }

    bool operator <=(const Point2D& rhs) const {
        return !(rhs < *this);
    }

    bool operator >(const Point2D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Point2D& rhs) const {
        return !(*this < rhs);
    }

    const Point2D& operator +=(const Vector2D<T>& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    const Point2D& operator -=(const Vector2D<T>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    /**
     * Return the i-th component value of point.
     */
    const T& operator[] (int i) const {
        return i == 0 ? x : y;
    }

    /**
     * Return the reference value of the i-th component of point.
     */
    T& operator[] (int i) {
        return i == 0 ? x : y;
    }

    /**
     * Return the i-th component value of point.
     */
    const T& at(int i) const {
        assert(0 <= i && i < 2);

        return i == 0 ? x : y;
    }

    /**
     * Return the reference value of the i-th component of point.
     */
    T& at(int i) {
        assert(0 <= i && i < 2);

        return i == 0 ? x : y;
    }

    /**
     * @return the bounding box of this point.
     */
    Box2D<T> bounding_box() const {
        return Box2D<T>(x, x, y, y);
    }

    /**
     * @return the dimension.
     */
    int size() const {
        return 2;
    }

    /**
     * Convert the point to vector.
     */
    Vector2D<T> ToVector() const {
        return Vector2D<T>(x, y);
    }

    friend Point2D operator +(const Point2D& lhs, const Vector2D<T>& rhs) {
        return Point2D<T>(lhs.x + rhs.x, lhs.y + rhs.y);
    }

    friend Point2D operator -(const Point2D& lhs, const Vector2D<T>& rhs) {
        return Point2D<T>(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    friend Vector2D<T> operator -(const Point2D& lhs, const Point2D& rhs) {
        return Vector2D<T>(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    friend std::ostream& operator <<(std::ostream& os, const Point2D& rhs) {
        os << rhs.ToVector();
        return os;
    }

    T x; /// X coordinate.
    T y; /// Y coordinate.
};

typedef Point2D<int>    IPoint2D;
typedef Point2D<double> RPoint2D;

} // namespace cl

#endif // GEOMETRY_KERNEL_POINT_2D_H_
