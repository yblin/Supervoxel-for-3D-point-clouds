//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_REGRESSION_LINEAR_LEAST_SQUARES_FITTING_H_
#define STATISTICS_REGRESSION_LINEAR_LEAST_SQUARES_FITTING_H_

#include <algorithm>
#include <cassert>

#include "codelibrary/geometry/kernel/line_2d.h"
#include "codelibrary/geometry/kernel/line_3d.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/util/compare_2d.h"
#include "codelibrary/geometry/util/compare_3d.h"
#include "codelibrary/geometry/transform/project_3d.h"
#include "codelibrary/statistics/principal_component_analysis_2d.h"
#include "codelibrary/statistics/principal_component_analysis_3d.h"

namespace cl {
namespace statistics {

/**
 * Linear Least Square Fitting for 2D line.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last, RLine2D* line) {
    assert(first != last);
    assert(line);

    PrincipalComponentAnalysis2D pca(first, last);

    *line = RLine2D(pca.centroid(), pca.eigenvectors()[0]);
    return pca.eigenvalues()[1] / pca.eigenvalues()[0];
}

/**
 * Weighted Linear Least Square Fitting for 2D line.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last,
                                 const Array<double>& weights,
                                 RLine2D* line) {
    assert(first != last);

    PrincipalComponentAnalysis2D pca(first, last, weights);

    if (pca.eigenvalues()[0] == pca.eigenvalues()[1]) {
        // A default line along X axis which goes through the centroid.
        *line = RLine2D(pca.centroid(), RVector2D(1.0, 0.0));
        return 1.0;
    } else {
        *line = RLine2D(pca.centroid(), pca.eigenvectors()[0]);
        return pca.eigenvalues()[1] / pca.eigenvalues()[0];
    }
}

/**
 * Linear Least Square Fitting for 2D line segment.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last,
                                 RSegment2D* seg) {
    assert(first != last);

    RLine2D line;
    double metric = LinearLeastSquaresFitting(first, last, &line);

    geometry::PointDirectionCompare2D compare(line.direction());
    RPoint2D p1 = geometry::Project(*std::min_element(first, last, compare),
                                    line);
    RPoint2D p2 = geometry::Project(*std::max_element(first, last, compare),
                                    line);
    *seg = RSegment2D(p1, p2);

    return metric;
}

/**
 * Weighted Linear Least Square Fitting for 2D line segment.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last,
                                 const Array<double>& weights,
                                 RSegment2D* seg) {
    assert(first != last);

    RLine2D line;
    double metric = LinearLeastSquaresFitting(first, last, weights, &line);

    geometry::PointDirectionCompare2D compare(line.direction());
    RPoint2D p1 = geometry::Project(*std::min_element(first, last, compare),
                                    line);
    RPoint2D p2 = geometry::Project(*std::max_element(first, last, compare),
                                    line);
    *seg = RSegment2D(p1, p2);

    return metric;
}

/**
 * Linear Least Square Fitting for 3D line.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last, RLine3D* line) {
    assert(first != last);

    PrincipalComponentAnalysis3D pca(first, last);

    if (pca.eigenvalues()[0] == pca.eigenvalues()[1] &&
        pca.eigenvalues()[0] == pca.eigenvalues()[2]) {
        // A default line along X axis which goes through the centroid.
        *line = RLine3D(pca.centroid(), RVector3D(1.0, 0.0, 0.0));
        return 1.0;
    } else {
        *line = RLine3D(pca.centroid(), pca.eigenvectors()[0]);
        return pca.eigenvalues()[1] / pca.eigenvalues()[0];
    }
}

/**
 * Linear Least Square Fitting for 3D segment.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last,
                                 RSegment3D* seg) {
    assert(first != last);

    RLine3D line;
    double metric = LinearLeastSquaresFitting(first, last, &line);

    geometry::PointDirectionCompare3D compare(line.direction());
    RPoint3D p1 = geometry::Project(*std::min_element(first, last, compare),
                                    line);
    RPoint3D p2 = geometry::Project(*std::max_element(first, last, compare),
                                    line);
    *seg = RSegment3D(p1, p2);

    return metric;
}

/**
 * Linear Least Square Fitting for 3D plane.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last,
                                 RPlane3D* plane) {
    assert(first != last);

    PrincipalComponentAnalysis3D pca(first, last);

    if (pca.eigenvalues()[0] == pca.eigenvalues()[1] &&
        pca.eigenvalues()[1] == pca.eigenvalues()[2]) {
        // A default horizontal plane that goes through the centroid.
        *plane = RPlane3D(pca.centroid(), RVector3D(0.0, 0.0, 1.0));
        return 1.0;
    } else {
        if (pca.eigenvectors()[2].squared_norm() == 0.0) {
            RVector3D v = CrossProduct(pca.eigenvectors()[0],
                                       pca.eigenvectors()[1]);
            if (v.squared_norm() == 0.0) {
                *plane = RPlane3D(pca.centroid(), RVector3D(0.0, 0.0, 1.0));
            } else {
                *plane = RPlane3D(pca.centroid(), v);
            }
        } else {
            *plane = RPlane3D(pca.centroid(), pca.eigenvectors()[2]);
        }
        return pca.eigenvalues()[2] / pca.eigenvalues()[1];
    }
}

/**
 * Weighted Linear Least Square Fitting for 3D plane.
 *
 * @return fitting quality: 0 is best, 1 is worst.
 */
template <typename Iterator>
double LinearLeastSquaresFitting(Iterator first, Iterator last,
                                 const Array<double>& weights,
                                 RPlane3D* plane) {
    assert(first != last);

    PrincipalComponentAnalysis3D pca(first, last, weights);

    if (pca.eigenvalues()[0] == pca.eigenvalues()[1] &&
        pca.eigenvalues()[1] == pca.eigenvalues()[2]) {
        // A default horizontal plane that goes through the centroid.
        *plane = RPlane3D(pca.centroid(), RVector3D(0.0, 0.0, 1.0));
        return 1.0;
    } else {
        if (pca.eigenvectors()[2].squared_norm() == 0.0) {
            RVector3D v = CrossProduct(pca.eigenvectors()[0],
                                       pca.eigenvectors()[1]);
            if (v.squared_norm() == 0.0) {
                *plane = RPlane3D(pca.centroid(), RVector3D(0.0, 0.0, 1.0));
            } else {
                *plane = RPlane3D(pca.centroid(), v);
            }
        } else {
            *plane = RPlane3D(pca.centroid(), pca.eigenvectors()[2]);
        }
        return pca.eigenvalues()[2] / pca.eigenvalues()[1];
    }
}

} // namespace statistics
} // namespace cl

#endif // STATISTICS_REGRESSION_LINEAR_LEAST_SQUARES_FITTING_H_
