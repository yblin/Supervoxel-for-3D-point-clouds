//
// Copyright 2017 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef MATH_MATRIX_TRIPLET_MATRIX_H_
#define MATH_MATRIX_TRIPLET_MATRIX_H_

#include <cassert>

#include "codelibrary/base/array.h"

namespace cl {

/**
 * Matrix in triplet form. It is mainly used for data transmission of sparse
 * matrix.
 */
template <typename T>
class TripletMatrix {
public:
    struct Triplet {
        Triplet() {}

        Triplet(int r, int c, const T& v)
            : row(r), column(c), value(v) {}

        int row, column;
        T value;
    };

    TripletMatrix()
        : n_rows_(0), n_columns_(0) {}

    TripletMatrix(int n_rows, int n_columns)
        : n_rows_(n_rows), n_columns_(n_columns) {}

    template <typename Iterator>
    TripletMatrix(int n_rows, int n_columns, Iterator first, Iterator last) {
        assert(n_rows >= 0 && n_columns >= 0);

        n_rows_ = n_rows;
        n_columns_ = n_columns;
        for (Iterator p = first; p != last; ++p) {
            InsertTriplet(p->row, p->column, p->value);
        }
    }

    /**
     * Construct the matrix from triplets. The dimensions of matrix are
     * automatically set.
     */
    template <typename Iterator>
    TripletMatrix(Iterator first, Iterator last) {
        n_rows_ = 0;
        n_columns_ = 0;
        for (Iterator p = first; p != last; ++p) {
            assert(p->row >= 0 && p->column >= 0);

            n_rows_ = std::max(p->row + 1, n_rows_);
            n_columns_ = std::max(p->column + 1, n_rows_);
            triplets_.push_back(*p);
        }
    }

    /**
     * Resize the matrix.
     */
    void Resize(int n_rows, int n_columns) {
        assert(n_rows >= 0);
        assert(n_columns >= 0);

        if (n_rows_ > n_rows || n_columns_ > n_columns) {
            int k = 0;
            for (const Triplet& triplet : triplets_) {
                if (triplet.row < n_rows && triplet.column < n_columns) {
                    triplets_[k++] = triplet;
                }
            }
            triplets_.resize(k);
        }
        n_rows_ = n_rows;
        n_columns_ = n_columns;
    }

    /**
     * Insert a triplet.
     */
    void InsertTriplet(int row, int column, const T& value) {
        assert(row >= 0 && row < n_rows_);
        assert(column >= 0 && column < n_columns_);

        Triplet triplet(row, column, value);
        triplets_.push_back(triplet);
    }

    void clear() {
        n_rows_ = 0;
        n_columns_ = 0;
        triplets_.clear();
    }

    /**
     * @return all of triplets.
     */
    const Array<Triplet>& triplets() const {
        return triplets_;
    }

    int n_rows() const {
        return n_rows_;
    }

    int n_columns() const {
        return n_columns_;
    }

private:
    int n_rows_;
    int n_columns_;
    Array<Triplet> triplets_;
};

} // namespcae cl

#endif // MATH_MATRIX_TRIPLET_MATRIX_H_
