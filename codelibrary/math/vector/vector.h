//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_VECTOR_VECTOR_H_
#define MATH_VECTOR_VECTOR_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ostream>

#include "codelibrary/base/array.h"
#include "codelibrary/base/message.h"
#include "codelibrary/math/vector/vector_2d.h"
#include "codelibrary/math/vector/vector_3d.h"

namespace cl {

/// N-dimensional Vector.
template<typename T>
class Vector : public Array<T> {
    typedef typename Array<T>::iterator Iterator;
    typedef typename Array<T>::const_iterator ConstIterator;

public:
    typedef T value_type;

    Vector()
        : Array<T>() {}

    explicit Vector(const Vector2D<T>& v)
        : Array<T>(2) {
        this->begin()[0] = v.x;
        this->begin()[1] = v.y;
    }

    explicit Vector(const Vector3D<T>& v)
        : Array<T>(3) {
        this->begin()[0] = v.x;
        this->begin()[1] = v.y;
        this->begin()[2] = v.z;
    }

    explicit Vector(int size, const T& value = T())
        : Array<T>(size, value) {}

    template <typename InputIter>
    Vector(InputIter first, InputIter last)
        : Array<T>(first, last) {}

    Vector(std::initializer_list<T> list)
        : Array<T>(list) {}

    Vector& operator = (const Vector2D<T>& v) {
        this->resize(2);
        this->begin()[0] = v.x;
        this->begin()[1] = v.y;

        return *this;
    }

    Vector& operator = (const Vector3D<T>& v) {
        this->resize(3);
        this->begin()[0] = v.x;
        this->begin()[1] = v.y;
        this->begin()[2] = v.z;

        return *this;
    }

    Vector& operator +=(const Vector& rhs) {
        assert(this->size() == rhs.size());

        // Faster.
        Iterator p1 = this->begin();
        ConstIterator p2 = rhs.begin();
        for (; this->end() - p1 >= 4; p1 += 4, p2 += 4) {
            *(p1)     += *(p2);
            *(p1 + 1) += *(p2 + 1);
            *(p1 + 2) += *(p2 + 2);
            *(p1 + 3) += *(p2 + 3);
        }
        for (; p1 != this->end(); ++p1, ++p2)
            *p1 += *p2;

        return *this;
    }

    Vector& operator -=(const Vector& rhs) {
        assert(this->size() == rhs.size());

        Iterator p1 = this->begin();
        ConstIterator p2 = rhs.begin();
        for (; this->end() - p1 >= 4; p1 += 4, p2 += 4) {
            *(p1)     -= *(p2);
            *(p1 + 1) -= *(p2 + 1);
            *(p1 + 2) -= *(p2 + 2);
            *(p1 + 3) -= *(p2 + 3);
        }
        for (; p1 != this->end(); ++p1, ++p2)
            *p1 -= *p2;

        return *this;
    }

    Vector& operator *=(const T& rhs) {
        for (Iterator p = this->begin(); p != this->end(); ++p) {
            *p *= rhs;
        }
        return *this;
    }

    /**
     * @return the squared euclidean norm of the vector.
     */
    double squared_norm() const {
        return (*this) * (*this);
    }

    /**
     * @return the euclidean norm of the vector.
     */
    double norm() const {
        return std::sqrt(squared_norm());
    }

    /**
     * @return infinity norm of this vector.
     */
    T inifity_norm() const {
        T res = 0;
        for (ConstIterator p = this->begin(); p != this->end(); ++p) {
            res = std::max(res, std::abs(*p));
        }
        return res;
    }

    /**
     * Convert vector to vector2D, require 'size_ == 2'.
     */
    Vector2D<T> ToVector2D() const {
        assert(this->size() == 2);

        return Vector2D<T>(this->begin()[0], this->begin()[1]);
    }

    /**
     * Convert vector to vector3D, Require 'size_ == 3'.
     */
    Vector3D<T> ToVector3D() const {
        assert(this->size() == 3);

        return Vector3D<T>(this->begin()[0], this->begin()[1],
                           this->begin()[2]);
    }

    friend Vector operator +(const Vector& lhs, const Vector& rhs) {
        assert(lhs.size() == rhs.size());

        ConstIterator p1 = lhs.begin(), p2 = rhs.begin();
        Vector v(lhs.size());

        Iterator p = v.begin();
        for (; v.end() - p >= 4; p1 += 4, p2 += 4, p += 4) {
            *(p)     = *(p1)     + *(p2);
            *(p + 1) = *(p1 + 1) + *(p2 + 1);
            *(p + 2) = *(p1 + 2) + *(p2 + 2);
            *(p + 3) = *(p1 + 3) + *(p2 + 3);
        }
        for (; p != v.end(); ++p1, ++p2, ++p)
            *(p) = *(p1) + *(p2);

        return v;
    }

    friend Vector operator -(const Vector& lhs, const Vector& rhs) {
        assert(lhs.size() == rhs.size());

        ConstIterator p1 = lhs.begin(), p2 = rhs.begin();
        Vector v(lhs.size());

        Iterator p = v.begin();
        for (; v.end() - p >= 4; p1 += 4, p2 += 4, p += 4) {
            *(p)     = *(p1)     - *(p2);
            *(p + 1) = *(p1 + 1) - *(p2 + 1);
            *(p + 2) = *(p1 + 2) - *(p2 + 2);
            *(p + 3) = *(p1 + 3) - *(p2 + 3);
        }
        for (; p != v.end(); ++p1, ++p2, ++p)
            *(p) = *(p1) - *(p2);

        return v;
    }

    friend Vector operator -(const Vector& rhs) {
        ConstIterator p1 = rhs.begin();
        Vector v(rhs.size());

        Iterator p = v.begin();
        for (; v.end() - p >= 4; p1 += 4, p += 4) {
            *(p)     = -*(p1);
            *(p + 1) = -*(p1 + 1);
            *(p + 2) = -*(p1 + 2);
            *(p + 3) = -*(p1 + 3);
        }
        for (; p != v.end(); ++p1, ++p)
            *(p) = -*(p1);

        return v;
    }

    friend Vector operator *(const T& lhs, const Vector& rhs) {
        ConstIterator p1 = rhs.begin();
        Vector v(rhs.size());

        Iterator p = v.begin();
        for (; v.end() - p >= 4; p1 += 4, p += 4) {
            *(p)     = *(p1)     * lhs;
            *(p + 1) = *(p1 + 1) * lhs;
            *(p + 2) = *(p1 + 2) * lhs;
            *(p + 3) = *(p1 + 3) * lhs;
        }
        for (; p != v.end(); ++p1, ++p)
            *(p) = *(p1) * lhs;

        return v;
    }

    friend Vector operator *(const Vector& lhs, const T& rhs) {
        ConstIterator p1 = lhs.begin();
        Vector v(lhs.size());

        Iterator p = v.begin();
        for (; v.end() - p >= 4; p1 += 4, p += 4) {
            *(p)     = *(p1)     * rhs;
            *(p + 1) = *(p1 + 1) * rhs;
            *(p + 2) = *(p1 + 2) * rhs;
            *(p + 3) = *(p1 + 3) * rhs;
        }
        for (; p != v.end(); ++p1, ++p)
            *(p) = *(p1) * rhs;

        return v;
    }

    friend double operator *(const Vector& lhs, const Vector& rhs) {
        assert(lhs.size() == rhs.size());

        double sum = 0.0;
        ConstIterator p1 = lhs.begin(), p2 = rhs.begin();
        for (; lhs.end() - p1 >= 4; p1 += 4, p2 += 4) {
            sum += static_cast<double>(*(p1))     * *(p2);
            sum += static_cast<double>(*(p1 + 1)) * *(p2 + 1);
            sum += static_cast<double>(*(p1 + 2)) * *(p2 + 2);
            sum += static_cast<double>(*(p1 + 3)) * *(p2 + 3);
        }
        for (; p1 != lhs.end(); ++p1, ++p2)
            sum += static_cast<double>(*p1) * *p2;

        return sum;
    }

    /**
     * For Message class.
     */
    friend std::ostream& operator <<(std::ostream& os, const Vector& rhs) {
        Message message;
        message.PrintSequence(rhs.begin(), rhs.end());
        return os << message;
    }
};

typedef Vector<int> IVector;
typedef Vector<double> RVector;

} // namespace cl

#endif // MATH_VECTOR_VECTOR_H_
