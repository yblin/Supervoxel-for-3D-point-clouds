//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef BASE_EQUAL_H_
#define BASE_EQUAL_H_

#include <cassert>
#include <cstdint>
#include <type_traits>

namespace cl {

namespace internal {

/**
 * This class represents an single-precision IEEE floating-point number.
 *
 * The purpose of this class is to do more sophisticated number comparison.
 * Due to round-off error, etc, it's very unlikely that two floating-points will
 * be equal exactly. Hence a naive comparison by the == operation often doesn't
 * work.
 */
class Float {
    typedef uint32_t BitType;

public:
    static const int MAX_ULP = 4;
    static const BitType SIGN_MASK     = 0x80000000UL; // 0 bit.
    static const BitType EXPONENT_MASK = 0x7F800000UL; // 1 ~ 8 bits.
    static const BitType MANTISSA_MASK = 0x007FFFFFUL; // 9 ~ 31 bits.

    explicit Float(float num) {
        u_.value = num;
    }

    /**
     * @return the sign bit of this number.
     */
    BitType sign_bit() const { return SIGN_MASK & u_.bits; }

    /**
     * @return the exponent bits of this number.
     */
    BitType exponent_bits() const { return EXPONENT_MASK & u_.bits; }

    /**
     * @return the mantissa bits of this number (9 ~ 31 bits).
     */
    BitType mantissa_bits() const { return MANTISSA_MASK & u_.bits; }

    /**
     * @return true iff this is NAN (not a number).
     */
    bool is_nan() const {
        return exponent_bits() == EXPONENT_MASK && mantissa_bits() != 0;
    }

    /**
     * Return true iff this number is at most MAX_ULP ULP's away from rhs.
     *
     * In particular, this function:
     *   - returns false if either number is (or both are) NAN.
     *   - treats really large numbers as almost equal to infinity.
     *   - thinks +0.0 and -0.0 are 0 DLP's apart.
     */
    bool operator == (const Float& rhs) const {
        // The IEEE standard says that any comparison operation involving
        // a NAN must return false.
        if (is_nan() || rhs.is_nan()) return false;

        return ULPDistance(u_.bits, rhs.u_.bits) <= MAX_ULP;
    }

private:
    /**
     * Convert an integer from the sign-and-magnitude representation to the
     * biased representation.
     */
    static BitType SignAndMagnitudeToBiased(BitType sam) {
        if (SIGN_MASK & sam) {
            // sam represents a negative number.
            return ~sam + 1;
        } else {
            // sam represents a positive number.
            return SIGN_MASK | sam;
        }
    }

    /**
     * Given two numbers in the sign-and-magnitude representation, return the
     * distance of ULP between them as an unsigned number.
     */
    static BitType ULPDistance(BitType sam1, BitType sam2) {
        BitType biased1 = SignAndMagnitudeToBiased(sam1);
        BitType biased2 = SignAndMagnitudeToBiased(sam2);
        return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
    }

    // The data type used to store the actual floating-point number.
    union FloatUnion {
        float value;   // The raw floating-point number.
        BitType bits;  // The bits that represent the number.
    };

    FloatUnion u_;
};

/**
 * This class represents an double-precision IEEE floating-point number.
 *
 * The purpose of this class is to do more sophisticated number comparison.
 * Due to round-off error, etc, it's very unlikely that two floating-points will
 * be equal exactly. Hence a naive comparison by the == operation often doesn't
 * work.
 */
class Double {
    typedef uint64_t BitType;

public:
    // The maximum error of a single floating-point operation is 0.5
    // units in the last place.  On Intel CPU's, all floating-point
    // calculations are done with 80-bit precision, while double has 64
    // bits.  Therefore, 4 should be enough for ordinary use.
    static const int MAX_ULP = 4;

    static const BitType SIGN_MASK     = 0x8000000000000000ULL; // 0 bit.
    static const BitType EXPONENT_MASK = 0x7FF0000000000000ULL; // 1 ~ 11 bits.
    static const BitType MANTISSA_MASK = 0x000FFFFFFFFFFFFFULL; // 12 ~ 63 bits.

    explicit Double(double num) {
        u_.value = num;
    }

    /**
     * @return the sign bit of this number.
     */
    BitType sign_bit() const { return SIGN_MASK & u_.bits; }

    /**
     * @return the exponent bits of this number.
     */
    BitType exponent_bits() const { return EXPONENT_MASK & u_.bits; }

    /**
     * @return the mantissa bits of this number (9 ~ 31 bits).
     */
    BitType mantissa_bits() const { return MANTISSA_MASK & u_.bits; }

    /**
     * @return true iff this is NAN (not a number).
     */
    bool is_nan() const {
        return exponent_bits() == EXPONENT_MASK && mantissa_bits() != 0;
    }

    /**
     * Return true iff this number is at most MAX_ULP ULP's away from rhs.
     *
     * In particular, this function:
     *   - returns false if either number is (or both are) NAN.
     *   - treats really large numbers as almost equal to infinity.
     *   - thinks +0.0 and -0.0 are 0 DLP's apart.
     */
    bool operator == (const Double& rhs) const {
        // The IEEE standard says that any comparison operation involving
        // a NAN must return false.
        if (is_nan() || rhs.is_nan()) return false;

        return ULPDistance(u_.bits, rhs.u_.bits) <= MAX_ULP;
    }

private:
    /**
     * Convert an integer from the sign-and-magnitude representation to the
     * biased representation.
     */
    static BitType SignAndMagnitudeToBiased(BitType sam) {
        if (SIGN_MASK & sam) {
            // sam represents a negative number.
            return ~sam + 1;
        } else {
            // sam represents a positive number.
            return SIGN_MASK | sam;
        }
    }

    /**
     * Given two numbers in the sign-and-magnitude representation, return the
     * distance of ULP between them as an unsigned number.
     */
    static BitType ULPDistance(BitType sam1, BitType sam2) {
        BitType biased1 = SignAndMagnitudeToBiased(sam1);
        BitType biased2 = SignAndMagnitudeToBiased(sam2);
        return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
    }

    // The data type used to store the actual floating-point number.
    union DoubleUnion {
        double value;  // The raw floating-point number.
        BitType bits;  // The bits that represent the number.
    };

    DoubleUnion u_;
};

} // namespace internal

/**
 * Check if two objects are equal.
 */
template <typename T1, typename T2>
inline bool Equal(const T1& lhs, const T2& rhs) {
    return lhs == rhs;
}

/**
 * Return true iff 'lhs' is at most 4 away from 'rhs'.
 *
 * In particular, this function:
 *   - returns false if either number is (or both are) NAN;
 *   - treats really large numbers as almost equal to infinity;
 *   - thinks +0.0 and -0.0 are 0 DLP's apart.
 */
inline bool Equal(float lhs, float rhs) {
    return internal::Float(lhs) == internal::Float(rhs);
}

/**
 * Similar to the previous one.
 */
inline bool Equal(double lhs, double rhs) {
    return internal::Double(lhs) == internal::Double(rhs);
}

/**
 * Check if the difference between lhs and rhs doesn't exceed the given absolute
 * error.
 */
template <typename T>
bool Equal(const T& lhs, const T& rhs, const T& absolute_error) {
    static_assert(std::is_arithmetic<T>::value,
                  "template argument is not a numeric type");

    return (lhs - rhs <= absolute_error) && (rhs - lhs <= absolute_error);
}

/**
 * Check if two ranges are equal.
 */
template <typename T1, typename T2>
bool Equal(T1 first1, T1 last1, T2 first2, T2 last2) {
    for (; first1 != last1 || first2 != last2; ++first1, ++first2) {
        if ((first1 == last1 || first2 == last2)) {
            return false;
        }

        if (!Equal(*first1, *first2)) return false;
    }

    return true;
}

} // namespace cl

#endif // BASE_EQUAL_H_
