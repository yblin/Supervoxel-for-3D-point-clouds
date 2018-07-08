//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_MATRIX_DECOMPOSE_TRIDIAGONAL_DECOMPOSE_H_
#define MATH_MATRIX_DECOMPOSE_TRIDIAGONAL_DECOMPOSE_H_

#include <cassert>
#include <limits>

#include "codelibrary/math/matrix/matrix.h"
#include "codelibrary/math/matrix/matrix_util.h"

namespace cl {
namespace matrix {

/**
 * A symmetric matrix A can be factorized by similarity transformations into the
 * form:
 *      A = U T U',
 * where U is an orthogonal matrix and T is a symmetric tridiagonal matrix.
 *
 * See HessenbergDecompose for more details.
 *
 * @param[in]  mat_a        - the given symmetric matrix.
 * @param[out] diagonal     - the diagonal of the resulting tridiagonal matrix.
 * @param[out] off_diagonal - the off-diagonal starting at off_diagonal[1]
 *                            (off_diagonal[0] is set to 0).
 * @param[out] mat_u        - the transformation matrix, U, where A = UTU'.
 */
template <typename T>
void TridiagonalDecompose(const Matrix<T>& mat_a,
                          Array<double>* diagonal,
                          Array<double>* off_diagonal,
                          RMatrix* mat_u = NULL) {
    assert(!mat_a.empty());
    assert(IsSymmetric(mat_a));
    assert(diagonal);

    int n = mat_a.n_rows();
    diagonal->resize(n);
    off_diagonal->resize(n);

    Matrix<T> u = mat_a;

    for (int i = n - 1; i >= 1; --i) {
        int l = i - 1;
        double h = 0.0, scale = 0.0;
        if (l > 0) {
            for (int k = 0; k <= l; ++k) {
                scale += std::abs(u(i, k));
            }
            if (scale == 0) {
                (*off_diagonal)[i] = u(i, l);
            } else {
                for (int k = 0; k <= l; ++k) {
                    u(i, k) /= scale;
                    h += u(i, k) * u(i, k);
                }
                double f = u(i, l);
                double g = f >= 0.0 ? -std::sqrt(h) : std::sqrt(h);
                (*off_diagonal)[i] = scale * g;
                h -= f * g;
                u(i, l) = f - g;

                f = 0.0;
                for (int j = 0; j <= l; ++j) {
                    u(j, i) = u(i, j) / h;

                    g = 0.0;
                    for (int k = 0; k <= j; ++k) {
                        g += u(j, k) * u(i, k);
                    }
                    for (int k = j + 1; k <= l; ++k) {
                        g += u(k, j) * u(i, k);
                    }
                    (*off_diagonal)[j] = g / h;
                    f += (*off_diagonal)[j] * u(i, j);
                }

                double hh = f / (h + h);
                for (int j = 0; j <= l; ++j) {
                    f = u(i, j);
                    (*off_diagonal)[j] = g = (*off_diagonal)[j] - hh * f;
                    for (int k = 0; k <= j; ++k) {
                        u(j, k) -= (f * (*off_diagonal)[k] + g * u(i, k));
                    }
                }
            }
        } else {
            (*off_diagonal)[i] = u(i, l);
        }
        (*diagonal)[i] = h;
    }

    (*diagonal)[0] = 0.0;
    (*off_diagonal)[0] = 0.0;

    for (int i = 0; i < n; ++i) {
        if ((*diagonal)[i] != 0.0) {
            for (int j = 0; j < i; ++j) {
                double g = 0.0;
                for (int k = 0; k < i; ++k) {
                    g += u(i, k) * u(k, j);
                }
                for (int k = 0; k < i; ++k) {
                    u(k, j) -= g * u(k, i);
                }
            }
        }

        (*diagonal)[i] = u(i, i);
        u(i, i) = 1.0;

        for (int j = 0; j < i; ++j) {
            u(j, i) = u(i, j) = 0.0;
        }
    }

    if (mat_u) *mat_u = u;
}

/**
 * Similar to the previous one, but store the diagonal and off_diagonal in the
 * matrix.
 *
 * @param[in]  mat_a - the given real symmetric matrix.
 * @param[out] mat_t - the resulting tridiagonal matrix.
 * @param[out] mat_u - the (optional) transformation matrix, U, where A = UTU'.
 */
template <typename T>
void TridiagonalDecompose(const Matrix<T>& mat_a, RMatrix* mat_t,
                          RMatrix* mat_u = NULL) {
    assert(!mat_a.empty());
    assert(matrix::IsSymmetric(mat_a));
    assert(mat_t);

    int n = mat_a.n_rows();

    Array<double> diagonal, off_diagonal;
    TridiagonalDecompose(mat_a, &diagonal, &off_diagonal, mat_u);

    mat_t->Assign(n, n, 0.0);
    for (int i = 0; i < n; ++i) {
        (*mat_t)(i, i) = diagonal[i];
        if (i > 0) {
            (*mat_t)(i - 1, i) = off_diagonal[i];
        }
        if (i < n - 1) {
            (*mat_t)(i + 1, i) = off_diagonal[i + 1];
        }
    }
}

} // namespace matrix
} // namespace cl

#endif // MATH_MATRIX_DECOMPOSE_TRIDIAGONAL_DECOMPOSE_H_
