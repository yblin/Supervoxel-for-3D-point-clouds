//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_POINT_CLOUD_PCA_ESTIMATE_NORMALS_H_
#define GEOMETRY_POINT_CLOUD_PCA_ESTIMATE_NORMALS_H_

#include <algorithm>
#include <cassert>
#include <cfloat>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif // _USE_MATH_DEFINES

#include "codelibrary/base/array.h"
#include "codelibrary/statistics/regression/linear_least_squares_fitting.h"
#include "codelibrary/util/metric/squared_euclidean.h"
#include "codelibrary/util/tree/kd_tree.h"

namespace cl {
namespace geometry {
namespace point_cloud {

/**
 * Estimates normal direction by weighted linear least squares firtting of a
 * plane over the input points.
 *
 * The output normal is normalize to the unit length, and its orientation is
 * randomly assigned.
 *
 * Because we only need to compute the least eigen value and eigen vector of the
 * covariance matrix, we don't choose PCA in the implementation.
 */
template <typename Iterator>
void PCAEstimateNormal(Iterator first, Iterator last,
                       const Array<double>& weights, RVector3D* normal) {
    assert(first != last);
    assert(normal);

    RPoint3D centroid = geometry::Centroid3D(first, last, weights);
    double a00 = 0.0, a01 = 0.0, a02 = 0.0, a11 = 0.0, a12 = 0.0, a22 = 0.0;
    int i = 0;
    double sum = 0.0;
    for (Iterator p = first; p != last; ++p, ++i) {
        double x = p->x - centroid.x;
        double y = p->y - centroid.y;
        double z = p->z - centroid.z;
        double w = weights[i];

        a00 += w * x * x;
        a01 += w * x * y;
        a02 += w * x * z;
        a11 += w * y * y;
        a12 += w * y * z;
        a22 += w * z * z;

        sum += w;
    }

    double t = 1.0 / sum;
    a00 = a00 * t;
    a01 = a01 * t;
    a02 = a02 * t;
    a11 = a11 * t;
    a12 = a12 * t;
    a22 = a22 * t;

    // Computing the least eigenvalue of the covariance matrix.
    double q = (a00 + a11 + a22) / 3.0;
    double pq = (a00 - q) * (a00 - q) + (a11 - q) * (a11 - q) +
                (a22 - q) * (a22 - q) +
                2.0 * (a01 * a01 + a02 * a02 + a12 * a12);
    pq = std::sqrt(pq / 6.0);
    double mpq = std::pow(1.0 / pq, 3.0);
    double det_b = mpq * ((a00 - q) * ((a11 - q) * (a22 - q) - a12 * a12) -
                          a01 * ( a01 * (a22 - q) - a12 * a02) +
                          a02 * ( a01 * a12 - (a11 - q) * a02));
    double r = 0.5 * det_b;
    double phi = 0.0;
    if (r <= -1.0)
        phi = M_PI / 3.0;
    else if (r >= 1.0)
        phi = 0.0;
    else
        phi = std::acos(r) / 3.0;
    double eig = q + 2.0 * pq * std::cos(phi + M_PI * (2.0 / 3.0));

    // Computing the corresponding eigenvector.
    normal->x = a01 * a12 - a02 * (a11 - eig);
    normal->y = a01 * a02 - a12 * (a00 - eig);
    normal->z = (a00 - eig) * (a11 - eig) - a01 * a01;

    // Normalize.
    double norm = normal->norm();
    if (norm == 0.0) {
        *normal = RVector3D(0.0, 0.0, 1.0);
    } else {
        *normal *= 1.0 / norm;
    }
}

/**
 * Estimates normal direction by linear least squares firtting of a plane over
 * the input points.
 *
 * The output normal is normalize to the unit length, and its orientation is
 * randomly assigned.
 */
template <typename Iterator>
void PCAEstimateNormal(Iterator first, Iterator last, RVector3D* normal) {
    Array<double> weights(std::distance(first, last), 1.0);
    PCAEstimateNormal(first, last, weights, normal);
}

/**
 * Estimates normal directions of the [first, last) range of 3D points by
 * linear least squares fitting of a plane over the k nearest neighbors.
 *
 * The output normals are normalize to the unit length, and their orientation
 * are randomly assigned.
 *
 * @param[in]  kd_tree  - the input points are stored in the KD tree.
 * @param[in]  k        - used to define the k-nearest neighbors.
 * @param[out] normals  - the output normals.
 */
template <typename Point>
void PCAEstimateNormals(const KDTree<Point>& kd_tree, int k,
                        Array<RVector3D>* normals) {
    assert(!kd_tree.empty());
    assert(k > 0);
    assert(normals);

    int n = kd_tree.size();
    k = std::min(k, n);

    normals->resize(n);

    const Array<Point>& points = kd_tree.points();
    Array<Point> neighbors;

    for (int i = 0; i < n; ++i) {
        kd_tree.FindKNearestNeighbors(points[i], k, &neighbors);
        PCAEstimateNormal(neighbors.begin(), neighbors.end(), &(*normals)[i]);
    }
}

/**
 * @see PCAEstimateNormals(const KDTree<Point>&, int, Array<RVector3D>*).
 */
template <typename Iterator>
void PCAEstimateNormals(Iterator first, Iterator last, int k,
                        Array<RVector3D>* normals) {
    typedef typename std::iterator_traits<Iterator>::value_type Point;
    KDTree<Point> kd_tree(first, last);
    PCAEstimateNormals(kd_tree, k, normals);
}

/**
 * Orientation aware PCA normal estimation.
 *
 * This function will re-orient the normal vector for each point according to
 * the neighbors whose normal vector has the same orientation.
 *
 * @param[in]  kd_tree          - the input points are stored in the KD tree.
 * @param[in]  k                - used to define the k-nearest neighbors.
 * @param[out] normals          - the output normals.
 */
template <typename Point>
void OrientationAwarePCAEstimateNormals(const KDTree<Point>& kd_tree, int k,
                                        Array<RVector3D>* normals) {
    assert(!kd_tree.empty());
    assert(k > 0);
    assert(normals);
    assert(normals->size() == kd_tree.size());

    int n = kd_tree.size();
    k = std::min(k, n);

    const Array<Point>& points = kd_tree.points();
    Array<int> neighbors;

    for (int i = 0; i < n; ++i) {
        kd_tree.FindKNearestNeighbors(points[i], k, &neighbors);

        Array<Point> neighbor_points;
        neighbor_points.reserve(k);
        for (int j = 0; j < k; ++j) {
            if ((*normals)[i] * (*normals)[neighbors[j]] >= 0.0) {
                neighbor_points.push_back(points[neighbors[j]]);
            }
        }

        RVector3D normal;
        PCAEstimateNormal(neighbor_points.begin(), neighbor_points.end(),
                          &normal);
        if (normal * (*normals)[i] < 0.0) {
            (*normals)[i] = -normal;
        } else {
            (*normals)[i] =  normal;
        }
    }
}

/**
 * Similar to the previous one.
 */
template <typename Iterator>
void OrientationAwarePCAEstimateNormals(Iterator first, Iterator last, int k,
                                        Array<RVector3D>* normals) {
    typedef typename std::iterator_traits<Iterator>::value_type Point;

    KDTree<Point> kd_tree(first, last);
    OrientationAwarePCAEstimateNormals(kd_tree, k, normals);
}

} // namespace point_cloud
} // namespace geometry
} // namespace cl

#endif // GEOMETRY_POINT_CLOUD_PCA_ESTIMATE_NORMALS_H_
