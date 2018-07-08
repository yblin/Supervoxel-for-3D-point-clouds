//
// Copyright 2015 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//
// Description: Printf variants that place their output in a C++ string.
//
// Usage:
//    string result = StringPrintf("%d %s\n", 10, "hello");
//    SStringPrintf(&result, "%d %s\n", 10, "hello");

#ifndef BASE_STRING_PRINTF_H_
#define BASE_STRING_PRINTF_H_

#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <string>

namespace cl {

/**
 * Lower-level routine that takes a va_list and appends to a specified string.
 * All other routines are just convenience wrappers around it.
 */
inline void StringAppendV(std::string* dst, const char* format, va_list ap) {
    // First try with a small fixed size buffer.
    char space[1024];

    // It's possible for methods that use a va_list to invalidate
    // the data in it upon use.  The fix is to make a copy
    // of the structure before using it and use that copy instead.
    va_list backup_ap;
    va_copy(backup_ap, ap);
    int result = vsnprintf(space, sizeof(space), format, backup_ap);
    va_end(backup_ap);

    if (static_cast<size_t>(result) < sizeof(space)) {
        if (result >= 0) {
            // Normal case -- everything fit.
            dst->append(space, result);
            return;
        }

        if (result < 0) {
            // Just an error.
            return;
        }
    }

    // Increase the buffer size to the size requested by vsnprintf,
    // plus one for the closing \0.
    int length = result + 1;
    char* buf = new (std::nothrow) char[length];
    assert(buf && "Memory is not enough.");

    // Restore the va_list before we use it again
    va_copy(backup_ap, ap);
    result = vsnprintf(buf, length, format, backup_ap);
    va_end(backup_ap);

    if (result >= 0 && result < length) {
        dst->append(buf, result);
    }
    delete[] buf;
}

/**
 * @return a C++ string.
 */
inline std::string StringPrintf(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    std::string result;
    StringAppendV(&result, format, ap);
    va_end(ap);
    return result;
}

/**
 * Store result into a supplied string and return it.
 */
inline const std::string& SStringPrintf(std::string* dst,
                                        const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    dst->clear();
    StringAppendV(dst, format, ap);
    va_end(ap);
    return *dst;
}

} // namespace cl

#endif // BASE_STRING_PRINTF_H_
