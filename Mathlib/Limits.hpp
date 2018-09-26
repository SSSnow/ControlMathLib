/****************************************************************************
 *
 *   Copyright (c) 2013 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file Limits.hpp
 *
 * Limiting / constrain helper functions
 */

#ifndef __LIMITS_HPP__
#define __LIMITS_HPP__

#ifndef M_PI
#define M_PI (3.14159265358979323846f)
#endif

//this should be defined in stdint.h, but seems to be missing in the ARM toolchain (5.2.0)
#ifndef UINT64_C
# if __WORDSIZE == 64
#  define UINT64_C(c)	c ## UL
# else
#  define UINT64_C(c)	c ## ULL
# endif
#endif


namespace math
{

template<typename _Tp>
const _Tp &min(const _Tp &a, const _Tp &b)
{
	return (a < b) ? a : b;
}

template<typename _Tp>
const _Tp &max(const _Tp &a, const _Tp &b)
{
	return (a > b) ? a : b;
}

template<typename _Tp>
const _Tp &constrain(const _Tp &val, const _Tp &min_val, const _Tp &max_val)
{
	return (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
}

template<typename _Tp>
inline bool isInRange(const _Tp &val, const _Tp &min_val, const _Tp &max_val)
{
	return (min_val <= val) && (val <= max_val);
}

template<typename T>
const T radians(const T degrees)
{
	return degrees * (static_cast<T>(M_PI) / static_cast<T>(180));
}

template<typename T>
const T degrees(const T radians)
{
	return radians * (static_cast<T>(180) / static_cast<T>(M_PI));
}

}

#endif
