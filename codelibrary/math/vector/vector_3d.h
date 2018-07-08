//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_VECTOR_VECTOR_3D_H_
#define MATH_VECTOR_VECTOR_3D_H_

#include <cassert>
#include <cmath>

#include "codelibrary/base/message.h"

namespace cl {

/// 3D Vector.
template<typename T>
class Vector3D {
public:
    typedef T value_type;

    Vector3D()
        : x(0), y(0), z(0) {}

    Vector3D(const T& x1, const T& y1, const T& z1)
        : x(x1), y(y1), z(z1) {}

    bool operator ==(const Vector3D& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator !=(const Vector3D& rhs) const {
        return !(*this == rhs);
    }

    bool operator < (const Vector3D& rhs) const {
        return x == rhs.x ? (y == rhs.y ? z < rhs.z : y < rhs.y) : x < rhs.x;
    }

    bool operator <=(const Vector3D& rhs) const {
        return !(rhs < *this);
    }

    bool operator > (const Vector3D& rhs) const {
        return rhs < *this;
    }

    bool operator >=(const Vector3D& rhs) const {
        return !(*this < rhs);
    }

    const Vector3D& operator +=(const Vector3D& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    const Vector3D& operator -=(const Vector3D& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    const Vector3D& operator *=(const T& rhs) {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    const Vector3D& operator /=(const T& rhs) {
        assert(rhs != 0);

        x /= rhs;
        y /= rhs;
        z /= rhs;
        return *this;
    }

    /**
     * @return the squared euclidean norm of the vector.
     */
    double squared_norm() const {
        return static_cast<double>(x) * x + static_cast<double>(y) * y +
               static_cast<double>(z) * z;
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
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the reference value of the i-th component of vector.
     */
    T& operator[] (int i) {
        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the i-th component value of vector.
     */
    const T& at(int i) const {
        assert(0 <= i && i < 3);

        return (i == 0) ? x : (i == 1 ? y : z);
    }

    /**
     * Return the reference value of the i-th component of vector.
     */
    T& at(int i) {
        assert(0 <= i && i < 3);

        return (i == 0) ? x : (i == 1 ? y : z);
    }

    int size() const { return 3; }

    friend Vector3D operator +(const Vector3D& lhs, const Vector3D& rhs) {
        return Vector3D(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
    }

    friend Vector3D operator -(const Vector3D& lhs, const Vector3D& rhs) {
        return Vector3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    friend Vector3D operator -(const Vector3D& rhs) {
        return Vector3D(-rhs.x, -rhs.y, -rhs.z);
    }

    friend Vector3D operator *(const T& lhs, const Vector3D& rhs) {
        return Vector3D(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
    }

    friend Vector3D operator *(const Vector3D& lhs, const T& rhs) {
        return Vector3D(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
    }

    friend double operator *(const Vector3D& lhs, const Vector3D& rhs) {
        return static_cast<double>(lhs.x) * rhs.x +
               static_cast<double>(lhs.y) * rhs.y +
               static_cast<double>(lhs.z) * rhs.z;
    }

    /**
     * For Message class.
     */
    friend std::ostream& operator <<(std::ostream& os, const Vector3D& rhs) {
        Message message;
        message << "{ " << rhs.x << " " << rhs.y << " " << rhs.z << " }";
        os << message;
        return os;
    }

    /**
     * @return the cross product of two vectors.
     *
     * The cross product of two vectors a and b is defined only in
     * three-dimensional space.
     *
     * Here we define the function outside this class to make sure this function
     * under the CL namespace. Otherwise, this function will under the global
     * namespace.
     */
    template <typename T1>
    friend Vector3D<T1> CrossProduct(const Vector3D<T1>& v1,
                                     const Vector3D<T1>& v2);

    T x; /// X component.
    T y; /// Y component.
    T z; /// Z component.
};

typedef Vector3D<int> IVector3D;
typedef Vector3D<double> RVector3D;

/**
 * @return the cross product of two vectors.
 */
template <typename T>
Vector3D<T> CrossProduct(const Vector3D<T>& v1, const Vector3D<T>& v2) {
    Vector3D<T> v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}

} // namespace cl

#endif // MATH_VECTOR_VECTOR_3D_H_
