//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_VECTOR_VECTOR_2D_H_
#define MATH_VECTOR_VECTOR_2D_H_

#include <cassert>
#include <cmath>

#include "codelibrary/base/message.h"

namespace cl {

/// 2D Vector.
template<typename T>
class Vector2D {
public:
    typedef T value_type;

    Vector2D()
        : x(0), y(0) {}

    Vector2D(const T& x1, const T& y1)
        : x(x1), y(y1) {}

    bool operator ==(const Vector2D& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator !=(const Vector2D& rhs) const {
        return !(*this == rhs);
    }

    bool operator < (const Vector2D& rhs) const {
        return x < rhs.x || (x == rhs.x && y < rhs.y);
    }

    bool operator <=(const Vector2D& rhs) const {
        return !(rhs < *this);
    }

    bool operator > (const Vector2D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Vector2D& rhs) const {
        return !(*this < rhs);
    }

    const Vector2D& operator +=(const Vector2D& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    const Vector2D& operator -=(const Vector2D& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    const Vector2D& operator *=(const T& rhs) {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    const Vector2D& operator /=(const T& rhs) {
        assert(rhs != 0);

        x /= rhs;
        y /= rhs;
        return *this;
    }

    /**
     * @return the squared euclidean norm of the vector.
     */
    double squared_norm() const {
        return static_cast<double>(x) * x + static_cast<double>(y) * y;
    }

    /**
     * @return the euclidean norm of the vector.
     */
    double norm() const {
        return std::sqrt(squared_norm());
    }

    /**
     * Return the i-th component value of vector.
     */
    const T& operator[] (int i) const {
        return i == 0 ? x : y;
    }

    /**
     * Return the reference value of the i-th component of vector.
     */
    T& operator[] (int i) {
        return i == 0 ? x : y;
    }

    /**
     * Return the i-th component value of vector.
     */
    const T& at(int i) const {
        assert(0 <= i && i < 2);

        return i == 0 ? x : y;
    }

    /**
     * Return the reference value of the i-th component of vector.
     */
    T& at(int i) {
        assert(0 <= i && i < 2);

        return i == 0 ? x : y;
    }

    int size() const { return 2; }

    friend Vector2D operator +(const Vector2D& lhs, const Vector2D& rhs) {
        return Vector2D(lhs.x + rhs.x, lhs.y + rhs.y);
    }

    friend Vector2D operator -(const Vector2D& lhs, const Vector2D& rhs) {
        return Vector2D(lhs.x - rhs.x, lhs.y - rhs.y);
    }

    friend Vector2D operator -(const Vector2D& rhs) {
        return Vector2D(-rhs.x, -rhs.y);
    }

    friend Vector2D operator *(const T& lhs, const Vector2D& rhs) {
        return Vector2D(lhs * rhs.x, lhs * rhs.y);
    }

    friend Vector2D operator *(const Vector2D& lhs, const T& rhs) {
        return Vector2D(lhs.x * rhs, lhs.y * rhs);
    }

    friend double operator *(const Vector2D& lhs, const Vector2D& rhs) {
        return static_cast<double>(lhs.x) * rhs.x +
               static_cast<double>(lhs.y) * rhs.y;
    }

    /**
     * For Message class.
     */
    friend std::ostream& operator <<(std::ostream& os, const Vector2D& rhs) {
        Message message;
        message << "{ " << rhs.x << " " << rhs.y << " }";
        os << message;
        return os;
    }

    T x; /// X component.
    T y; /// Y component.
};

typedef Vector2D<int> IVector2D;
typedef Vector2D<double> RVector2D;

} // namespace cl

#endif // MATH_VECTOR_VECTOR_2D_H_
