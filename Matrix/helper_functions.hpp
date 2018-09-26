#ifndef __HELPER_FUNCTIONS_HPP__
#define __HELPER_FUNCTIONS_HPP__

#include "math.hpp"

#if defined (__PX4_NUTTX) || defined (__PX4_QURT)
#include <px4_defines.h>  //TODO: define a defines.h of ourself
#endif

#define M_PI (3.14159265358979323846f)

namespace matrix
{

template<typename Type>
bool is_finite(Type x) {
#if defined (__PX4_NUTTX)
    return PX4_ISFINITE(x);
#elif defined (__PX4_QURT)
    return __builtin_isfinite(x);
#else
    return isfinite(x);
#endif
}

template<typename Type>
Type wrap_pi(Type x)
{
    if (!is_finite(x)) {
        return x;
    }

    int c = 0;

    while (x >= Type(M_PI)) {
        x -= Type(2 * M_PI);

        if (c++ > 100) {
            return INFINITY;
        }
    }

    c = 0;

    while (x < Type(-M_PI)) {
        x += Type(2 * M_PI);

        if (c++ > 100) {
            return INFINITY;
        }
    }

    return x;
}

template<typename Type>
Type wrap_2pi(Type x)
{
    if (!is_finite(x)) {
        return x;
    }

    int c = 0;

    while (x >= Type(2 * M_PI)) {
        x -= Type(2 * M_PI);

        if (c++ > 100) {
            return INFINITY;
        }
    }

    c = 0;

    while (x < Type(0)) {
        x += Type(2 * M_PI);

        if (c++ > 100) {
            return INFINITY;
        }
    }

    return x;
}

}

#endif
