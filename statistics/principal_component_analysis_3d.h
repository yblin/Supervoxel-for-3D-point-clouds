//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_PRINCIPAL_COMPONENT_ANALYSIS_3D_H_
#define STATISTICS_PRINCIPAL_COMPONENT_ANALYSIS_3D_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/macros.h"
#include "codelibrary/geometry/kernel/plane_3d.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/util/center_3d.h"
#include "codelibrary/math/matrix.h"

namespace cl {
namespace statistics {

/// Principal Component Analysis for 3D data.
/**
 * Principal Component Analysis (PCA) is an orthogonal linear transformation
 * that transforms the data to a new coordinate system such that the greatest
 * variance by some projection of the data comes to lie on the first coordinate
 * (called the first principal component), the second greatest variance on the
 * second coordinate, and so on.
 */
class PrincipalComponentAnalysis3D {
public:
    /**
     * Compute PCA by covariance method.
     */
    template <typename Iterator>
    PrincipalComponentAnalysis3D(Iterator first, Iterator last) {
        assert(first != last);

        // Get the covariance matrix.
        centroid_ = geometry::Centroid3D(first, last);
        double a00 = 0.0, a01 = 0.0, a02 = 0.0, a11 = 0.0, a12 = 0.0, a22 = 0.0;
        for (Iterator p = first; p != last; ++p) {
            double x = p->x - centroid_.x;
            double y = p->y - centroid_.y;
            double z = p->z - centroid_.z;

            a00 += x * x;
            a01 += x * y;
            a02 += x * z;
            a11 += y * y;
            a12 += y * z;
            a22 += z * z;
        }

        RMatrix mat(3, 3);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(0, 2) = mat(2, 0) = a02;
        mat(1, 1) = a11;
        mat(1, 2) = mat(2, 1) = a12;
        mat(2, 2) = a22;

        // Get the eigenvectors of covariance matrix.
        matrix::SymmetricEigen3(mat, &eigenvalues_, &eigenvectors_);
    }

    /**
     * Compute weighted PCA by covariance method.
     */
    template <typename Iterator>
    PrincipalComponentAnalysis3D(Iterator first, Iterator last,
                                 const Array<double>& weights) {
        assert(first != last);

        // Get the covariance matrix.
        centroid_ = geometry::Centroid3D(first, last, weights);
        double a00 = 0.0, a01 = 0.0, a02 = 0.0, a11 = 0.0, a12 = 0.0, a22 = 0.0;
        int i = 0;
        for (Iterator p = first; p != last; ++p, ++i) {
            double x = p->x - centroid_.x;
            double y = p->y - centroid_.y;
            double z = p->z - centroid_.z;
            double w = weights[i];

            a00 += w * x * x;
            a01 += w * x * y;
            a02 += w * x * z;
            a11 += w * y * y;
            a12 += w * y * z;
            a22 += w * z * z;
        }

        RMatrix mat(3, 3);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(0, 2) = mat(2, 0) = a02;
        mat(1, 1) = a11;
        mat(1, 2) = mat(2, 1) = a12;
        mat(2, 2) = a22;

        // Get the eigenvectors of covariance matrix.
        matrix::SymmetricEigen3(mat, &eigenvalues_, &eigenvectors_);
    }

    const RMatrix& covariance_matrix() const {
        return covariance_matrix_;
    }

    const Array<double>& eigenvalues() const {
        return eigenvalues_;
    }

    const Array<RVector3D>& eigenvectors() const {
        return eigenvectors_;
    }

    const RPoint3D& centroid() const {
        return centroid_;
    }

private:
    // The centroid point of the input data.
    RPoint3D centroid_;

    // Covariance matrix of the input data.
    RMatrix covariance_matrix_;

    // Eigenvalues of the covariance matrix, they are sorted in decreasing
    // order.
    Array<double> eigenvalues_;

    // Corresponding eigenvectors.
    Array<RVector3D> eigenvectors_;

    DISALLOW_COPY_AND_ASSIGN(PrincipalComponentAnalysis3D);
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_PRINCIPAL_COMPONENT_ANALYSIS_3D_H_
