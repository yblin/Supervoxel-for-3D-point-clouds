//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef STATISTICS_PRINCIPAL_COMPONENT_ANALYSIS_2D_H_
#define STATISTICS_PRINCIPAL_COMPONENT_ANALYSIS_2D_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/macros.h"
#include "codelibrary/geometry/kernel/point_2d.h"
#include "codelibrary/geometry/util/center_2d.h"
#include "codelibrary/math/matrix.h"

namespace cl {
namespace statistics {

/// Principal Component Analysis for 2D data.
/**
 * Principal Component Analysis (PCA) is an orthogonal linear transformation
 * that transforms the data to a new coordinate system such that the greatest
 * variance by some projection of the data comes to lie on the first coordinate
 * (called the first principal component), the second greatest variance on the
 * second coordinate, and so on.
 */
class PrincipalComponentAnalysis2D {
public:
    /**
     * Compute PCA by covariance method.
     */
    template <typename Iterator>
    PrincipalComponentAnalysis2D(Iterator first, Iterator last) {
        assert(first != last);

        // Get the covariance matrix.
        centroid_ = geometry::Centroid2D(first, last);
        double a00 = 0.0, a01 = 0.0, a11 = 0.0;
        for (Iterator p = first; p != last; ++p) {
            double x = p->x - centroid_.x;
            double y = p->y - centroid_.y;

            a00 += x * x;
            a01 += x * y;
            a11 += y * y;
        }

        RMatrix mat(2, 2);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(1, 1) = a11;

        // Get the eigenstd::vectors of covariance matrix.
        matrix::SymmetricEigen2(mat, &eigenvalues_, &eigenvectors_);
    }

    /**
     * Compute weighted PCA by covariance method.
     */
    template <typename Iterator>
    PrincipalComponentAnalysis2D(Iterator first, Iterator last,
                                 const Array<double>& weights) {
        assert(first != last);

        // Get the covariance matrix.
        centroid_ = geometry::Centroid2D(first, last, weights);
        double a00 = 0.0, a01 = 0.0, a11 = 0.0;
        int i = 0;
        for (Iterator p = first; p != last; ++p, ++i) {
            double x = p->x - centroid_.x;
            double y = p->y - centroid_.y;
            double w = weights[i];

            a00 += w * x * x;
            a01 += w * x * y;
            a11 += w * y * y;
        }

        RMatrix mat(2, 2);
        mat(0, 0) = a00;
        mat(0, 1) = mat(1, 0) = a01;
        mat(1, 1) = a11;

        // Get the eigenstd::vectors of covariance matrix.
        matrix::SymmetricEigen2(mat, &eigenvalues_, &eigenvectors_);
    }

    const RMatrix& covariance_matrix() const {
        return covariance_matrix_;
    }

    const Array<double>& eigenvalues() const {
        return eigenvalues_;
    }

    const Array<RVector2D>& eigenvectors() const {
        return eigenvectors_;
    }

    const RPoint2D& centroid() const {
        return centroid_;
    }

private:
    // The centroid point of the input data.
    RPoint2D centroid_;

    // Covariance matrix of the input data.
    RMatrix covariance_matrix_;

    // Eigenvalues of the covariance matrix, they are sorted in decreasing
    // order.
    Array<double> eigenvalues_;

    // Corresponding eigenstd::vectors.
    Array<RVector2D> eigenvectors_;

    DISALLOW_COPY_AND_ASSIGN(PrincipalComponentAnalysis2D);
};

} // namespace statistics
} // namespace cl

#endif // STATISTICS_PRINCIPAL_COMPONENT_ANALYSIS_2D_H_
