//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: Compute the eigen-system for symmetric matrix.
//
// Reference:
//   Kopp J. Efficient numerical diagonalization of hermitian 3x3 matrices[J].
//   International Journal of Modern Physics C, 2008, 19(03): 523-548.
//

#ifndef MATH_MATRIX_SYMMETRIC_EIGEN_H_
#define MATH_MATRIX_SYMMETRIC_EIGEN_H_

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <type_traits>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"
#include "codelibrary/base/equal.h"
#include "codelibrary/math/matrix/matrix.h"
#include "codelibrary/math/matrix/decompose/tridiagonal_decompose.h"

namespace cl {
namespace matrix {

/**
 * Compute the eigenvalues and eigenvectors for 2x2 real symmetric matrix.
 *
 * @param[in]  mat          - the given matrix.
 * @param[out] eigenvalues  - the eigenvalues of mat, they are sorted in
 *                            decreasing order.
 * @param[out] eigenvectors - the normalized eigenvectors of mat.
 */
template <typename T>
void SymmetricEigen2(const Matrix<T>& mat,
                     Array<double>* eigenvalues,
                     Array<RVector2D>* eigenvectors = NULL) {
    static_assert(std::is_floating_point<T>::value,
                  "template argument is not a floating point type");

    assert(mat.n_rows() == 2 && mat.n_columns() == 2);
    assert(IsSymmetric(mat));
    assert(eigenvalues);

    eigenvalues->resize(2);
    Array<T>& w = *eigenvalues;

    // [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
    // [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
    double a = mat(0, 0), b = mat(0, 1), c = mat(1, 1);
    double sm = a + c;
    double df = a - c;
    double rt = std::sqrt(df * df + b * b * 4.0);

    if (sm > 0.0) {
        w[0] = (sm + rt) / 2.0;
        double t = 1.0 / w[0];
        w[1] = (a * t) * c - (b * t) * b;
    } else if (sm < 0.0) {
        w[1] = (sm - rt) / 2.0;
        double t = 1.0 / w[1];
        w[0] = (a * t) * c - (b * t) * b;
    } else {
        // This case needs to be treated separately to avoid div by 0.
        w[0] =  rt / 2.0;
        w[1] = -rt / 2.0;
    }

    if (eigenvectors) {
        // Compute eigenvectors.
        double cs = df > 0.0 ? df + rt : df - rt;
        double t = 0.0, sn = 0.0;
        if (std::fabs(cs) > std::fabs(b) * 2.0) {
            t = -b * 2.0 / cs;
            sn = 1.0 / std::sqrt(t * t + 1);
            cs = t * sn;
        } else if (std::fabs(b) == 0) {
            cs = 1.0;
            sn = 0.0;
        } else {
            t = - cs / b / 2.0;
            cs = 1.0 / std::sqrt(t * t + 1);
            sn = t * cs;
        }

        if (df > 0.0) {
            t = cs;
            cs = -sn;
            sn = t;
        }

        eigenvectors->resize(2);
        (*eigenvectors)[0] = RVector2D(cs, sn);
        (*eigenvectors)[1] = RVector2D(-sn, cs);
    }
}

/**
 * Get the eigenvalues and eigenvectors for 3x3 real symmetric matrix.
 *
 * @param[in]  mat          - the given matrix.
 * @param[out] eigenvalues  - the eigenvalues of mat, they are sorted in
 *                            decreasing order.
 * @param[out] eigenvectors - the normalized eigenvectors of mat.
 */
template <typename T>
void SymmetricEigen3(const Matrix<T>& mat,
                     Array<double>* eigenvalues,
                     Array<RVector3D>* eigenvectors = NULL) {
    static_assert(std::is_floating_point<T>::value,
                  "template argument is not a floating point type");

    assert(mat.n_rows() == 3 && mat.n_columns() == 3);
    assert(IsSymmetric(mat));
    assert(eigenvalues);

    // Compute the eigenvalues.
    eigenvalues->resize(3);
    Array<double>& w = *eigenvalues;

    double a00 = mat(0, 0), a01 = mat(0, 1), a02 = mat(0, 2);
    double a11 = mat(1, 1), a12 = mat(1, 2), a22 = mat(2, 2);

    double de = a01 * a12;
    double dd = a01 * a01;
    double ee = a12 * a12;
    double ff = a02 * a02;

    double m  = a00 + a11 + a22;
    double c1 = a00 * a11 + a00 * a22 + a11 * a22 - (dd + ee + ff);
    double c0 = dd * a22 + ee * a00 + ff * a11 - a00 * a11 * a22 -
                a02 * de * 2.0;

    double p = m * m - c1 * 3.0;
    double q = m * (p - c1 * 1.5) - c0 * 13.5;
    double sqrt_p = std::sqrt(std::fabs(p));

    double phi = (c1 * c1 * (p - c1) * 0.25 + c0 * (q + c0 * 6.75)) * 27.0;
    phi = std::atan2(std::sqrt(std::fabs(phi)), q) / 3.0;

    double c = sqrt_p * std::cos(phi);
    double s = sqrt_p * std::sin(phi) / std::sqrt(3.0);

    w[1]  = (m - c) / 3.0;
    w[2]  = w[1] + s;
    w[0]  = w[1] + c;
    w[1] -= s;
    if (eigenvectors) {
        // Compute eigen vectors.
        eigenvectors->resize(3);

        Array<RVector3D>& v = *eigenvectors;

        double max_eigenvalue = std::max(std::fabs(w[0]), std::fabs(w[1]));
        max_eigenvalue = std::max(max_eigenvalue, std::fabs(w[2]));
        double epsilon = std::numeric_limits<double>::epsilon();
        double thresh = epsilon * max_eigenvalue * 8.0;
        thresh *= thresh;

        // Prepare calculation of eigenvectors.
        double n0tmp = a01 * a01 + a02 * a02;
        double n1tmp = a01 * a01 + a12 * a12;
        v[1][0] = a01 * a12 - a02 * a11;
        v[1][1] = a02 * a01 - a12 * a00;
        v[1][2] = a01 * a01;

        // Calculate first eigenvector by the formula:
        //   v[0] = (A - w[0]).e1 x (A - w[0]).e2.
        a00 -= w[0];
        a11 -= w[0];
        v[0][0] = v[1][0] + a02 * w[0];
        v[0][1] = v[1][1] + a12 * w[0];
        v[0][2] = a00 * a11 - v[1][2];

        double norm = v[0][0] * v[0][0] + v[0][1] * v[0][1] + v[0][2] * v[0][2];
        double n0    = n0tmp + a00 * a00;
        double n1    = n1tmp + a11 * a11;
        double error = n0 * n1;

        if (n0 <= thresh) {
            v[0][0] = 1.0;
            v[0][1] = 0.0;
            v[0][2] = 0.0;
        } else if (n1 <= thresh) {
            v[0][0] = 0.0;
            v[0][1] = 1.0;
            v[0][2] = 0.0;
        } else if (norm < (64.0 * epsilon * 64.0 * epsilon) * error) {
            // If angle between mat_t[0] and mat_t[1] is too small, don't use
            // cross product, but calculate v ~ (1, -A0/A1, 0).

            T t = a01 * a01;
            T f = -a00 / a01;
            if (a11 * a11 > t) {
                t =  a11 * a11;
                f = -a01 / a11;
            }
            if (a12 * a12 > t)
                f = -a02 / a12;

            double norm = 1.0 / std::sqrt(f * f + 1.0);
            v[0][0] = norm;
            v[0][1] = f * norm;
            v[0][2] = 0.0;
        } else {
            // This is the standard branch.
            v[0] *= std::sqrt(1.0 / norm);
        }

        // Prepare calculation of second eigenvector.
        double t = w[0] - w[1];
        if (std::fabs(t) > epsilon * max_eigenvalue * 8.0) {
            // For non-degenerate eigenvalue, calculate second eigenvector by
            // the formula
            //   v[1] = (A - w[1]).e1 x (A - w[1]).e2.
            a00 += t;
            a11 += t;
            v[1][0] = v[1][0] + a02 * w[1];
            v[1][1] = v[1][1] + a12 * w[1];
            v[1][2] = a00 * a11 - v[1][2];

            double norm  = v[1][0] * v[1][0] + v[1][1] * v[1][1] +
                           v[1][2] * v[1][2];
            double n0    = n0tmp + a00 * a00;
            double n1    = n1tmp + a11 * a11;
            double error = n0 * n1;

            if (n0 <= thresh) {
                v[1][0] = 1.0;
                v[1][1] = 0.0;
                v[1][2] = 0.0;
            } else if (n1 <= thresh) {
                v[1][0] = 0.0;
                v[1][1] = 1.0;
                v[1][2] = 0.0;
            } else if (norm < (epsilon * epsilon * 4096.0) * error) {
                // If angle between mat_t[0] and mat_t[1] is too small,
                // don't use cross product, but calculate v ~ (1, -A0/A1, 0).
                double t = a01 * a01;
                double f = -a00 / a01;
                if (a11 * a11 > t) {
                    t =  a11 * a11;
                    f = -a01 / a11;
                }
                if (a12 * a12 > t) f = -a02 / a12;

                double norm = 1.0 / std::sqrt(f * f + 1.0);
                v[1][0] = norm;
                v[1][1] = f * norm;
                v[1][2] = 0.0;
            } else {
                // This is the standard branch.
                v[1] *= std::sqrt(1.0 / norm);
            }
        } else {
            // For degenerate eigenvalue, calculate second eigenvector according
            // to:
            //   v[1] = v[0] x (A - w[1]).e[i].

            // Reset the mat_t to mat.
            a00 += w[0];
            a11 += w[0];
            RMatrix mat_t(3, 3);
            mat_t(0, 0) = a00, mat_t(0, 1) = a01, mat_t(0, 2) = a02;
            mat_t(1, 0) = a01, mat_t(1, 1) = a11, mat_t(1, 2) = a12;
            mat_t(2, 0) = a02, mat_t(2, 1) = a12, mat_t(2, 2) = a22;

            int i = 0;
            for (i = 0; i < 3; ++i) {
                mat_t(i, i) -= w[1];

                T n0 = mat_t(0, i) * mat_t(0, i) + mat_t(1, i) * mat_t(1, i) +
                       mat_t(2, i) * mat_t(2, i);
                if (n0 > thresh) {
                    v[1][0] = v[0][1] * mat_t(2, i) - v[0][2] * mat_t(1, i);
                    v[1][1] = v[0][2] * mat_t(0, i) - v[0][0] * mat_t(2, i);
                    v[1][2] = v[0][0] * mat_t(1, i) - v[0][1] * mat_t(0, i);

                    T norm = v[1][0] * v[1][0] + v[1][1] * v[1][1] +
                             v[1][2] * v[1][2];
                    if (norm > epsilon * epsilon * 65536.0 * n0) {
                        // Accept cross product only if the angle between the
                        // two vectors was not too small.
                        v[1] *= std::sqrt(1.0 / norm);
                        break;
                    }
                }
            }

            if (i == 3) {
                // This means that any vector orthogonal to v[0] is an EV.
                for (int j = 0; j < 3; ++j) {
                    if (v[0][j] != 0) {
                        // Find nonzero element of v[0] and swap it with the
                        // next one.
                        int k = (j + 1) % 3;
                        T inv_norm = T(1) / std::sqrt(v[0][j] * v[0][j] +
                                                      v[0][k] * v[0][k]);
                        v[1][j] =  v[0][k] * inv_norm;
                        v[1][k] = -v[0][j] * inv_norm;
                        v[1][(j + 2) % 3] = 0.0;
                        break;
                    }
                }
            }
        }

        // Calculate third eigenvector according to v[2] = v[0] x v[1].
        v[2][0] = v[0][1] * v[1][2] - v[0][2] * v[1][1];
        v[2][1] = v[0][2] * v[1][0] - v[0][0] * v[1][2];
        v[2][2] = v[0][0] * v[1][1] - v[0][1] * v[1][0];

        // Sort in decreasing order.
        if (w[0] < w[1]) {
            std::swap(w[0], w[1]);
            std::swap(v[0], v[1]);
        }
        if (w[1] < w[2]) {
            std::swap(w[1], w[2]);
            std::swap(v[1], v[2]);
            if (w[0] < w[1]) {
                std::swap(w[0], w[1]);
                std::swap(v[0], v[1]);
            }
        }
    } else {
        // Sort in decreasing order.
        if (w[0] < w[1])
            std::swap(w[0], w[1]);
        if (w[1] < w[2]) {
            std::swap(w[1], w[2]);
            if (w[0] < w[1]) {
                std::swap(w[0], w[1]);
            }
        }
    }
}

/**
 * Get the eigenvalues and eigenvectors for real symmetric matrix.
 *
 * @param[in]  mat          - the given matrix.
 * @param[out] eigenvalues  - the eigenvalues of mat, they are sorted by
 *                            decreasing order.
 * @param[out] eigenvectors - the normalized eigenvectors of mat.
 */
template <typename T>
bool SymmetricEigen(const Matrix<T>& mat,
                    Array<double>* eigenvalues,
                    Array<RVector>* eigenvectors = NULL) {
    static_assert(std::is_floating_point<T>::value,
                  "template argument is not a floating point type");

    assert(!mat.empty() && "You are using an empty matrix.");
    assert(IsSymmetric(mat));
    assert(eigenvalues);

    int n = mat.n_rows();
    if (n == 1) {
        eigenvalues->resize(1);
        (*eigenvalues)[0] = mat(0, 0);
        if (eigenvectors) {
            eigenvectors->resize(1);
            (*eigenvectors)[0] = Vector<T>(1, 1);
        }
        return true;
    } else if (n == 2) {
        Array<Vector2D<T> > vectors;
        SymmetricEigen2(mat, eigenvalues, &vectors);

        if (eigenvectors) {
            eigenvectors->resize(2);
            (*eigenvectors)[0] = vectors[0];
            (*eigenvectors)[1] = vectors[1];
        }

        return true;
    } else if (n == 3) {
        Array<Vector3D<T> > vectors;
        SymmetricEigen3(mat, eigenvalues, &vectors);

        if (eigenvectors) {
            eigenvectors->resize(3);
            (*eigenvectors)[0] = vectors[0];
            (*eigenvectors)[1] = vectors[1];
            (*eigenvectors)[2] = vectors[2];
        }

        return true;
    } else {
        // Using QL algorithm with implicit shifts for general case.
        // Transform the symmetric matrix into tridiagonal form.
        Array<double> diagonal, off_diagonal;
        // Used to store the eigenvectors.
        RMatrix z;
        TridiagonalDecompose(mat, &diagonal, &off_diagonal, &z);

        int n = mat.n_rows();

        for (int i = 0; i < n - 1; ++i) {
            off_diagonal[i] = off_diagonal[i + 1];
        }
        off_diagonal[n - 1] = 0.0;

        const int max_iter = 30;

        for (int l = 0; l < n - 1; ++l) {
            int iter = 0;
            while (true) {
                // Looking for a single small subdiagonal element to split the
                // matrix.
                int m;
                for (m = l; m < n - 1; ++m) {
                    double dd = std::fabs(diagonal[m]) +
                                std::fabs(diagonal[m + 1]);
                    if (Equal(std::fabs(off_diagonal[m]) + dd, dd)) break;
                }
                if (m == l) break;

                if (iter++ == max_iter) return false;

                double g = (diagonal[l + 1] - diagonal[l]) /
                           (off_diagonal[l] * 2);
                double r = std::sqrt(g * g + 1);
                double t = g > 0.0 ? r : -r;
                g = diagonal[m] - diagonal[l] + off_diagonal[l] / (g + t);

                double s = 1.0, c = 1.0;
                double p = 0.0;

                for (int i = m - 1; i >= l; --i) {
                    double f = s * off_diagonal[i];
                    double b = c * off_diagonal[i];
                    if (std::fabs(f) > std::fabs(g)) {
                        c = g / f;
                        r = std::sqrt(c * c + 1.0);
                        off_diagonal[i + 1] = f * r;
                        s = 1.0 / r;
                        c *= s;
                    } else {
                        s = f / g;
                        r = std::sqrt(s * s + 1.0);
                        off_diagonal[i + 1] = g * r;
                        c = 1.0 / r;
                        s *= c;
                    }

                    g = diagonal[i + 1] - p;
                    r = (diagonal[i] - g) * s + c * b * 2.0;
                    p = s * r;
                    diagonal[i + 1] = g + p;
                    g = c * r - b;

                    if (eigenvectors) {
                        for (int k = 0; k < n; ++k) {
                            double t = z(k, i + 1);
                            z(k, i + 1) = s * z(k, i) + c * t;
                            z(k, i)     = c * z(k, i) - s * t;
                        }
                    }
                }

                diagonal[l]    -= p;
                off_diagonal[l] = g;
                off_diagonal[m] = 0.0;
            }
        }

        // Sort the eigenvalues.
        if (eigenvectors) {
            Array<int> seq;
            IndexSort(diagonal.begin(), diagonal.end(), &seq);
            std::reverse(seq.begin(), seq.end());

            eigenvalues->resize(n);
            eigenvectors->resize(n);
            for (int i = 0; i < n; ++i) {
                (*eigenvalues)[i] = diagonal[seq[i]];
                (*eigenvectors)[i].resize(n);

                for (int j = 0; j < n; ++j) {
                    (*eigenvectors)[i][j] = z(j, seq[i]);
                }
                (*eigenvectors)[i] *= 1.0 / (*eigenvectors)[i].norm();
            }
        } else {
            std::sort(diagonal.begin(), diagonal.end(), std::greater<T>());
            *eigenvalues = diagonal;
        }

        return true;
    }
}

} // namespace matrix
} // namespace cl

#endif // MATH_MATRIX_SYMMETRIC_EIGEN_H_
