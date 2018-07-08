//
// Copyright 2012 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: 2d compare for geometry objects.
//

#ifndef GEOMETRY_UTIL_COMPARE_2D_H_
#define GEOMETRY_UTIL_COMPARE_2D_H_

#include "codelibrary/geometry/kernel/segment_2d.h"

namespace cl {
namespace geometry {

/**
 * Compare two 2D points along the given direction.
 */
class PointDirectionCompare2D {
public:
    template <typename T>
    explicit PointDirectionCompare2D(const Vector2D<T>& direction)
        : direction_(direction.x, direction.y) {}

    template <typename T>
    bool operator() (const Point2D<T>& p1, const Point2D<T>& p2) const {
        RVector2D v1(p1.x, p1.y);
        RVector2D v2(p2.x, p2.y);
        return v1 * direction_ < v2 * direction_;
    }

private:
    RVector2D direction_;
};

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_UTIL_COMPARE_2D_H_
