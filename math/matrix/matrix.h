//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_MATRIX_MATRIX_H_
#define MATH_MATRIX_MATRIX_H_

#include <algorithm>
#include <cassert>
#include <limits>
#include <iomanip>
#include <ostream>
#include <string>

#include "codelibrary/base/message.h"
#include "codelibrary/math/matrix/triplet_matrix.h"
#include "codelibrary/math/vector.h"
#include "codelibrary/util/array/array_2d.h"

namespace cl {

/// Matrix class.
template <typename T>
class Matrix : public Array2D<T> {
    static const int BLOCK_SIZE = 256;

public:
    typedef T ValueType;

    Matrix()
        : Array2D<T>() {}

    Matrix(int n_rows, int n_columns, const T& value = T())
        : Array2D<T>(n_rows, n_columns, value) {}

    template <typename InputIter>
    Matrix(int n_rows, int n_columns, InputIter first, InputIter last)
        : Array2D<T>(n_rows, n_columns, first, last) {}

    Matrix& operator += (const Matrix& rhs) { return Add(*this, rhs);      }
    Matrix& operator -= (const Matrix& rhs) { return Subtract(*this, rhs); }
    Matrix& operator *= (const T& rhs)      { return Multiply(*this, rhs); }
    Matrix& operator *= (const Matrix& rhs) { return Multiply(*this, rhs); }

    /**
     * Swap two rows.
     */
    void SwapRows(int row1, int row2) {
        assert(0 <= row1 && row1 < n_rows());
        assert(0 <= row2 && row2 < n_rows());

        if (row1 == row2) return;

        int offset1 = n_rows() * row1;
        int offset2 = n_rows() * row2;
        for (int i = 0; i < n_columns(); ++i) {
            std::swap(this->data_[offset1 + i], this->data_[offset2 + i]);
        }
    }

    /**
     * Swap two columns.
     */
    void SwapColumns(int column1, int column2) {
        assert(0 <= column1 && column1 < n_columns());
        assert(0 <= column2 && column2 < n_columns());

        if (column1 == column2) return;

        for (int i = 0; i < n_rows(); ++i) {
            std::swap(this->data_[column1], this->data_[column2]);
            column1 += n_columns();
            column2 += n_columns();
        }
    }

    /**
     * Set the matrix to be indentity.
     */
    void SetIdentity() {
        this->Fill(0);

        int min = std::min(n_rows(), n_columns());
        for (int i = 0; i < min; ++i) {
            this->operator()(i, i) = 1;
        }
    }

    /**
     * Return the transpose matrix.
     */
    Matrix Transpose() const {
        Matrix mat_transpose(n_columns(), n_rows());

        for (int i = 0; i < mat_transpose.n_rows(); ++i) {
            for (int j = 0; j < mat_transpose.n_columns(); ++j) {
                mat_transpose(i, j) = this->operator()(j, i);
            }
        }

        return mat_transpose;
    }

    /**
     * Return the trace of a square matrix.
     */
    T Trace() const {
        assert(n_rows() == n_columns());

        T trace = 0;
        for (int i = 0; i < n_rows(); ++i) {
            trace += this->operator()(i, i);
        }
        return trace;
    }

    /**
     * Return a row vector.
     */
    Vector<T> RowVector(int row) const {
        assert(0 <= row && row < n_rows());

        Vector<T> v(n_columns());
        int offset = n_columns() * row;
        for (int i = 0; i < n_columns(); ++i) {
            v[i] = this->data_[offset + i];
        }
        return v;
    }

    /**
     * Return a column vector.
     */
    Vector<T> ColumnVector(int column) const {
        assert(0 <= column && column < n_columns());

        Vector<T> v(n_rows());
        for (int i = 0, k = column; i < n_rows(); ++i, k += n_columns()) {
            v[i] = this->data_[k];
        }
        return v;
    }

    /**
     * Get the sub-matrix.
     *
     * @param r1, c1 - the left top coordinate of sub-matrix is (r1, c1).
     * @param r2, c2 - the right bottom coordinate of sub-matrix is (r2, c2).
     */
    Matrix SubMatrix(int r1, int c1, int r2, int c2) const {
        assert(0 <= r1 && r1 < n_rows());
        assert(0 <= r2 && r2 < n_rows());
        assert(0 <= c1 && c1 < n_columns());
        assert(0 <= c2 && c2 < n_columns());
        assert(r1 <= r2 && c1 <= c2);

        Matrix sub(r2 - r1 + 1, c2 - c1 + 1);
        int k = 0;
        for (int i = r1; i <= r2; ++i) {
            int offset = i * n_columns();
            for (int j = c1; j <= c2; ++j) {
                sub.data_[k++] = this->data_[offset + j];
            }
        }
        return sub;
    }

    /**
     * This = a + b.
     */
    Matrix& Add(const Matrix& a, const Matrix& b) {
        assert(a.n_rows() == b.n_rows() && a.n_columns() == b.n_columns());

        this->Resize(a.n_rows(), a.n_columns());
        for (int i = 0; i < this->size(); ++i) {
            this->data_[i] = a.data_[i] + b.data_[i];
        }
        return *this;
    }

    /**
     * This = a - b.
     */
    Matrix& Subtract(const Matrix& a, const Matrix& b) {
        assert(a.n_rows() == b.n_rows() && a.n_columns() == b.n_columns());

        this->Resize(a.n_rows(), a.n_columns());
        for (int i = 0; i < this->size(); ++i) {
            this->data_[i] = a.data_[i] - b.data_[i];
        }
        return *this;
    }

    /**
     * This = a * b.
     */
    Matrix& Multiply(const Matrix& a, const T& b) {
        this->Resize(a.n_rows(), a.n_columns());
        for (int i = 0; i < this->size(); ++i) {
            this->data_[i] = a.data_[i] * b;
        }
        return *this;
    }

    /**
     * This = a * b, use naive implementation.
     */
    Matrix& Multiply(const Matrix& a, const Matrix& b) {
        assert(a.n_columns() == b.n_rows());

        if (this == &a || this == &b) {
            Matrix c;
            BlockMultiply(a, b, &c);
            return *this = c;
        } else {
            BlockMultiply(a, b, this);
            return *this;
        }
    }

    /**
     * Convert matrix to triplets.
     */
    void ToTripletMatrix(TripletMatrix<T>* triplet_matrix) const {
        assert(triplet_matrix);

        triplet_matrix->clear();
        triplet_matrix->Resize(this->size1_, this->size2_);
        for (int i = 0; i < this->size1_; ++i) {
            for (int j = 0; j < this->size2_; ++j) {
                triplet_matrix->InsertTriplet(i, j, this->data_(i, j));
            }
        }
    }

    int n_rows()    const { return this->size1_; }
    int n_columns() const { return this->size2_; }

    friend Matrix operator -(const Matrix& rhs) {
        Matrix matrix(rhs.n_rows(), rhs.n_columns());
        for (int i = 0; i < rhs.size_; ++i) {
            matrix.data_[i] = -rhs.data_[i];
        }
        return matrix;
    }

    friend Matrix operator +(const Matrix& lhs, const Matrix& rhs) {
        Matrix t(lhs);
        t += rhs;
        return t;
    }

    friend Matrix operator -(const Matrix& lhs, const Matrix& rhs) {
        Matrix t(lhs);
        t -= rhs;
        return t;
    }

    friend Matrix operator *(const Matrix& lhs, const T& rhs) {
        Matrix t(lhs);
        t *= rhs;
        return t;
    }

    friend Matrix operator *(const T& lhs, const Matrix& rhs) {
        Matrix t(rhs);
        t *= lhs;
        return t;
    }

    friend Matrix operator *(const Matrix& lhs, const Matrix& rhs) {
        Matrix t;
        t.Multiply(lhs, rhs);
        return t;
    }

    friend Vector<T> operator *(const Matrix& lhs, const Vector<T>& rhs) {
        assert(lhs.n_columns() == rhs.size());

        int n = lhs.n_columns();

        Vector<T> res(lhs.n_rows(), 0);
        for (int i = 0; i < lhs.n_rows(); ++i) {
            // Using cache size 4 to accerate.
            int j = 0;
            for (; j + 3 < n; j += 4) {
                res[i] += lhs(i, j)     * rhs[j] +
                          lhs(i, j + 1) * rhs[j + 1] +
                          lhs(i, j + 2) * rhs[j + 2] +
                          lhs(i, j + 3) * rhs[j + 3];
            }
            for (; j < n; ++j) {
                res[i] += lhs(i, j) * rhs[j];
            }
        }

        return res;
    }

    friend Vector<T> operator *(const Vector<T>& lhs, const Matrix& rhs) {
        assert(lhs.size() == rhs.n_rows());

        int n = rhs.n_rows();
        int p = rhs.n_columns();

        Vector<T> res(p, 0);
        for (int i = 0; i < n; ++i) {
            int j = 0;
            for (; j + 3 < p; j += 4) {
                res[j]     += rhs(i, j)     * lhs[i];
                res[j + 1] += rhs(i, j + 1) * lhs[i];
                res[j + 2] += rhs(i, j + 2) * lhs[i];
                res[j + 3] += rhs(i, j + 3) * lhs[i];
            }
            for (; j < rhs.n_columns(); ++j) {
                res[j] += rhs(i, j) * lhs[i];
            }
        }

        return res;
    }

    /**
     * For Message class.
     */
    friend std::ostream& operator <<(std::ostream& os, const Matrix& rhs) {
        size_t max_width = 0;
        for (int i = 0; i < rhs.size(); ++i) {
            Message message;
            message << rhs.data_[i];
            max_width = std::max(max_width, message.ToString().length());
        }

        os << '[';
        for (int i = 0; i < rhs.n_rows(); ++i) {
            for (int j = 0; j < rhs.n_columns(); ++j) {
                if (i > 0 && j == 0) os << ' ';
                Message message;
                message << rhs(i, j);
                os << std::setw(max_width + 1) << message.ToString();
            }
            if (i + 1 < rhs.n_rows()) os << '\n';
        }
        os << " ]";

        return os;
    }

private:
    /**
     * Blocked matrix multiplication.
     *
     * It exploits cache to speed up. A default size of block is 256 (defined in
     * BLOCK_SIZE).
     * And loops are unrolled by 4 elements.
     */
    void BlockMultiply(const Matrix& a, const Matrix& b, Matrix* c) const {
        assert(c);

        int block_size = BLOCK_SIZE;

        int m = a.size1();
        int n = a.size2();
        int p = b.size2();
        c->Resize(m, p);
        c->Fill(0);

        for (int ii = 0; ii < m; ii += block_size) {
            int t1 = std::min(m, ii + block_size);
            for (int kk = 0; kk < n; kk += block_size) {
                int t2 = std::min(n, kk + block_size);
                for (int jj = 0; jj < p; jj += block_size) {
                    int t3 = std::min(p, jj + block_size);
                    for (int i = ii; i < t1; ++i) {
                        int k = kk;
                        for (; k + 3 < t2; k += 4) {
                            T* c0 = &(*c)(i, jj);
                            const T* b0 = &b(k, jj);
                            const T* b1 = &b(k + 1, jj);
                            const T* b2 = &b(k + 2, jj);
                            const T* b3 = &b(k + 3, jj);

                            const T& tmp0 = a(i, k);
                            const T& tmp1 = a(i, k + 1);
                            const T& tmp2 = a(i, k + 2);
                            const T& tmp3 = a(i, k + 3);
                            int j = jj;
                            for (; j + 3 < t3; j += 4) {
                                *c0++ += tmp0 * *b0++ + tmp1 * *b1++ +
                                         tmp2 * *b2++ + tmp3 * *b3++;
                                *c0++ += tmp0 * *b0++ + tmp1 * *b1++ +
                                         tmp2 * *b2++ + tmp3 * *b3++;
                                *c0++ += tmp0 * *b0++ + tmp1 * *b1++ +
                                         tmp2 * *b2++ + tmp3 * *b3++;
                                *c0++ += tmp0 * *b0++ + tmp1 * *b1++ +
                                         tmp2 * *b2++ + tmp3 * *b3++;
                            }
                            for (; j < t3; ++j) {
                                *c0++ += tmp0 * *b0++ + tmp1 * *b1++ +
                                         tmp2 * *b2++ + tmp3 * *b3++;
                            }
                        }

                        for (; k < t2; ++k) {
                            T* c0 = &(*c)(i, jj);
                            const T* b0 = &b(k, jj);
                            T tmp = a(i, k);
                            int j = jj;
                            for (; j + 3 < t3; j += 4) {
                                *c0++ += tmp * *b0++;
                                *c0++ += tmp * *b0++;
                                *c0++ += tmp * *b0++;
                                *c0++ += tmp * *b0++;
                            }
                            for (; j < t3; ++j) {
                                *c0++ += tmp * *b0++;
                            }
                        }
                    }
                }
            }
        }
    }
};

typedef Matrix<int> IMatrix;
typedef Matrix<double> RMatrix;

} // namespace cl

#endif // MATH_MATRIX_MATRIX_H_
