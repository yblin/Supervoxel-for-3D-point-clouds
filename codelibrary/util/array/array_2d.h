//
// Copyright 2011 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_ARRAY_ARRAY_2D_H_
#define UTIL_ARRAY_ARRAY_2D_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>

#include "codelibrary/base/array.h"

namespace cl {

/// 2D Array.
/**
 * This class is as most efficient as C-style 2d array.
 * User can access the data of Array2D by () operator.
 */
template<typename T>
class Array2D {
public:
    typedef T*       Iterator;
    typedef const T* ConstIterator;

    Array2D()
        : size1_(0), size2_(0) {}

    Array2D(int size1, int size2, const T& value = T())
        : size1_(size1), size2_(size2) {
        CheckDimension(size1, size2);

        data_.resize(size1_ * size2_, value);
    }

    template <typename InputIterator>
    Array2D(int size1, int size2, InputIterator first, InputIterator last)
        : size1_(size1), size2_(size2), data_(first, last) {
        CheckDimension(size1_, size2_);

        assert(size1 * size2 == data_.size());
    }

    /**
     * Clear the data of 2d array and set rows and columns to zero.
     */
    void clear() {
        data_.clear();
        size1_ = 0;
        size2_ = 0;
    }

    /**
     * Resize the 2d array.
     */
    void Resize(int size1, int size2) {
        CheckDimension(size1, size2);

        if (size2 == size2_) {
            size1_ = size1;
            data_.resize(size1_ * size2_);
        } else {
            Array<T> data = data_;
            data_.resize(size1 * size2);

            int min_size1 = std::min(size1_, size1);
            int min_size2 = std::min(size2_, size2);
            for (int i = 0; i < min_size1; ++i) {
                std::copy_n(data.begin() + i * size2_, min_size2,
                            data_.begin() + i * size2);
            }
            size1_ = size1;
            size2_ = size2;
        }
    }

    /**
     * Resize the 2d array with given value.
     */
    void Resize(int size1, int size2, const T& value) {
        CheckDimension(size1, size2);

        if (size2 == size2_) {
            size1_ = size1;
            data_.resize(size1_ * size2_, value);
        } else {
            Array<T> data = data_;
            data_.assign(size1 * size2, value);

            int min_size1 = std::min(size1_, size1);
            int min_size2 = std::min(size2_, size2);
            for (int i = 0; i < min_size1; ++i) {
                std::copy_n(data.begin() + i * size2_, min_size2,
                            data_.begin() + i * size2);
            }
            size1_ = size1;
            size2_ = size2;
        }
    }

    /**
     * Return the first dimension of 2d array.
     */
    int size1() const {
        return size1_;
    }

    /**
     * Return the second dimension of 2d array.
     */
    int size2() const {
        return size2_;
    }

    /**
     * Return the total elements of 2d array.
     */
    int size() const {
        return data_.size();
    }

    /**
     * Return the reference of raw data.
     */
    T* data() {
        return data_.data();
    }

    /**
     * Return the const reference of raw data.
     */
    const T* data() const {
        return data_.data();
    }

    /**
     * Return true if 2d array is empty.
     */
    bool empty() const {
        return data_.empty();
    }

    /**
     * Assign the data to v.
     */
    void Assign(int size1, int size2, const T& v) {
        CheckDimension(size1, size2);
        data_.assign(size1 * size2, v);
        size1_ = size1;
        size2_ = size2;
    }

    /**
     * Fill the data of 2d array with given value.
     */
    void Fill(const T& v) {
        std::fill(data_.begin(), data_.end(), v);
    }

    /**
     * Swap 2D array.
     */
    void Swap(Array2D* v) {
        assert(v);

        std::swap(size1_, v->size1_);
        std::swap(size2_, v->size2_);
        data_.Swap(v->data_);
    }

    /**
     * This operator allows for easy, array-style, data access.
     * Note that, we do not check the range of index here. Use at() for checked
     * access.
     *
     * @return the value at position (i, j).
     */
    const T& operator() (int i, int j) const {
        return data_[i * size2_ + j];
    }

    /**
     * This operator allows for easy, array-style, data access.
     * Note that, we do not check the range of index here. Use at() for checked
     * access.
     *
     * @return the reference at position (i, j).
     */
    T& operator() (int i, int j) {
        return data_[i * size2_ + j];
    }

    /**
     * @return the value at position (i, j).
     */
    const T& at(int i, int j) const {
        assert(0 <= i && i < size1_);
        assert(0 <= j && j < size2_);

        return data_[i * size2_ + j];
    }

    /**
     * @return the reference at position (i, j).
     */
    T& at(int i, int j) {
        assert(0 <= i && i < size1_);
        assert(0 <= j && j < size2_);

        return data_[i * size2_ + j];
    }

    Iterator begin()            { return data_.begin(); }
    Iterator end()              { return data_.end();   }
    ConstIterator begin() const { return data_.begin(); }
    ConstIterator end()   const { return data_.end();   }

protected:
    /**
     * Check if the given dimension is valid.
     */
    static void CheckDimension(int size1, int size2) {
        assert(size1 >= 0 && "size1 cannot be negative.");
        assert(size2 >= 0 && "size2 cannot be negative.");
        assert((size2 == 0 || size1 <= INT_MAX / size2) &&
               "The given dimensions are too large.");
    }

    int size1_;     // The first dimension.
    int size2_;     // The second dimension.
    Array<T> data_; // The data of 2d array.
};

} // namespace cl

#endif // UTIL_ARRAY_ARRAY_2D_H_
