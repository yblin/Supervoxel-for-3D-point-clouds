//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_MATRIX_MATRIX_UTIL_H_
#define MATH_MATRIX_MATRIX_UTIL_H_

#include "codelibrary/math/matrix/matrix.h"

namespace cl {
namespace matrix {

/**
 * Return true if the given matrix is symmetric (self-adjoint).
 */
template <typename T>
bool IsSymmetric(const Matrix<T>& mat) {
    if (mat.n_rows() != mat.n_columns()) return false;

    int n = mat.n_rows();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (mat(i, j) != mat(j, i)) return false;
        }
    }
    return true;
}

/**
 * Return true if the given matrix is a upper Hessenberg matrix.
 */
template <typename T>
bool IsUpperHessenberg(const Matrix<T>& mat) {
    if (mat.n_rows() != mat.n_columns()) return false;

    int n = mat.n_rows();
    for (int i = 2; i < n; ++i) {
        for (int j = 0; j + 1 < i; ++j) {
            if (mat(i, j) != 0) return false;
        }
    }
    return true;
}

/**
 * Balancing a matrix for calculation of eigenvalues and eigenvectors.
 *
 * Parlett B N, Reinsch C. Balancing a matrix for calculation of eigenvalues and
 * eigenvectors [J]. Numerische Mathematik, 1969, 13(4):293-304.
 */
template <typename T>
void Balance(const Matrix<T>& a, RMatrix* balanced) {
    assert(balanced);
    assert(a.n_rows() == a.n_columns());

    const int n = a.n_rows();
    RMatrix off_diagonal(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            off_diagonal(i, j) = a(i, j);
        }
        off_diagonal(i, i) = 0.0;
    }

    // gamma <= 1 controls how much a change in the scaling has to lower the
    // 1-norm of the matrix to be accepted.
    //
    // gamma = 1 seems to lead to cycles (numerical issues), so we set it
    // slightly lower.
    const double gamma = 0.9;

    // Greedily scale row/column pairs until there is no change.
    bool scaling_has_changed;
    do {
        scaling_has_changed = false;

        for (int i = 0; i < n; ++i) {
            double row_norm = 0.0, col_norm = 0.0;
            for (int j = 0; j < n; ++j) {
                row_norm += std::fabs(off_diagonal(i, j));
                col_norm += std::fabs(off_diagonal(j, i));
            }

            // Decompose row_norm/col_norm into mantissa * 2^exponent,
            // where 0.5 <= mantissa < 1. Discard mantissa (return value of
            // frexp), as only the exponent is needed.
            int exponent = 0;
            std::frexp(row_norm / col_norm, &exponent);
            exponent /= 2;

            if (exponent != 0) {
                const double scaled_col_norm = std::ldexp(col_norm, exponent);
                const double scaled_row_norm = std::ldexp(row_norm, -exponent);
                if (scaled_col_norm + scaled_row_norm <
                    gamma * (col_norm + row_norm)) {
                    // Accept the new scaling. (Multiplication by powers of 2
                    // should not introduce rounding errors (ignoring
                    // non-normalized numbers and over- or underflow))
                    scaling_has_changed = true;
                    double scale1 = std::ldexp(1.0, -exponent);
                    double scale2 = std::ldexp(1.0,  exponent);
                    for (int j = 0; j < n; ++j) {
                        off_diagonal(i, j) *= scale1;
                        off_diagonal(j, i) *= scale2;
                    }
                }
            }
        }
    } while (scaling_has_changed);

    for (int i = 0; i < n; ++i) {
        off_diagonal(i, i) = a(i, i);
    }
    *balanced = off_diagonal;
}

} // namespace matrix
} // namespace cl

#endif // MATH_MATRIX_MATRIX_UTIL_H_
