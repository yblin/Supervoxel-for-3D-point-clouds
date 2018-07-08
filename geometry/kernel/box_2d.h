//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_KERNEL_BOX_2D_H_
#define GEOMETRY_KERNEL_BOX_2D_H_

#include <algorithm>
#include <cassert>
#include <limits>

namespace cl {

/// 2D Axis-Aligned Minimum Bounding Box.
/**
 * The axis-aligned minimum bounding box for a given point set is its minimum
 * bounding box subject to the constraint that the edges of the box are parallel
 * to the (Cartesian) coordinate axis.
 */
template <typename T>
class Box2D {
public:
    typedef T value_type;

    /**
     * The default box is an invalid box.
     */
    Box2D()
        : x_min_(std::numeric_limits<T>::max()),
          x_max_(std::numeric_limits<T>::lowest()),
          y_min_(std::numeric_limits<T>::max()),
          y_max_(std::numeric_limits<T>::lowest()) {}

    /**
     * Require lower bounds no greater than upper bounds.
     */
    Box2D(const T& x_min, const T& x_max, const T& y_min, const T& y_max)
        : x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max) {
        assert(x_min_ <= x_max_ && y_min_ <= y_max_);
    }

    template <typename Iterator>
    Box2D(Iterator first, Iterator last)
        : Box2D() {
        for (Iterator i = first; i != last; ++i) {
            x_min_ = std::min(x_min_, (*i)[0]);
            x_max_ = std::max(x_max_, (*i)[0]);
            y_min_ = std::min(y_min_, (*i)[1]);
            y_max_ = std::max(y_max_, (*i)[1]);
        }
    }

    bool empty() const {
        return x_min_ > x_max_;
    }

    /**
     * The lower x value of box.
     */
    const T& x_min() const {
        return x_min_;
    }

    /**
     * The higher x value of box.
     */
    const T& x_max() const {
        return x_max_;
    }

    /**
     * The lower y value of box.
     */
    const T& y_min() const {
        return y_min_;
    }

    /**
     * The higher y value of box.
     */
    const T& y_max() const {
        return y_max_;
    }

    bool operator ==(const Box2D& rhs) const {
        return x_min_ == rhs.x_min_ && x_max_ == rhs.x_max_ &&
               y_min_ == rhs.y_min_ && y_max_ == rhs.y_max_;
    }

    bool operator !=(const Box2D& rhs) const {
        return x_min_ != rhs.x_min_ || x_max_ != rhs.x_max_ ||
               y_min_ != rhs.y_min_ || y_max_ != rhs.y_max_;
    }

    /**
     * @return the X length of box.
     * @note that x_length may be negative if x_max_ and x_min_ are extreme.
     */
    T x_length() const {
        return empty() ? 0.0 : x_max_ - x_min_;
    }

    /**
     * @return the Y length of box.
     * @note that y_length may be negative if y_max_ and y_min_ are extreme.
     */
    T y_length() const {
        return empty() ? 0.0 : y_max_ - y_min_;
    }

    /**
     * @return the bounding box of this box (itself).
     */
    const Box2D<T>& bounding_box() const {
        return *this;
    }

    /**
     * Join this box with anthor box, the result is the hull of two boxes.
     */
    void Join(const Box2D& box) {
        x_min_ = std::min(x_min_, box.x_min_);
        y_min_ = std::min(y_min_, box.y_min_);
        x_max_ = std::max(x_max_, box.x_max_);
        y_max_ = std::max(y_max_, box.y_max_);
    }

    /**
     * Intersect this box with anthor box, the result is the intersection of two
     * boxes.
     */
    void Intersect(const Box2D& box) {
        x_min_ = std::max(x_min_, box.x_min_);
        y_min_ = std::max(y_min_, box.y_min_);
        x_max_ = std::min(x_max_, box.x_max_);
        y_max_ = std::min(y_max_, box.y_max_);
    }

    /**
     * @return the minimum value of the i-th dimension.
     */
    T min(int i) const {
        assert(i >= 0 && i < 2);
        return i == 0 ? x_min_ : y_min_;
    }

    /**
     * @return the maximum value of the i-th dimension.
     */
    T max(int i) const {
        assert(i >= 0 && i < 2);
        return i == 0 ? x_max_ : y_max_;
    }

protected:
    T x_min_; // The lower x value of box.
    T x_max_; // The higher x value of box.
    T y_min_; // The lower y value of box.
    T y_max_; // The higher y value of box.
};

typedef Box2D<int>    IBox2D;
typedef Box2D<double> RBox2D;

} // namespace cl

#endif // GEOMETRY_KERNEL_BOX_2D_H_
