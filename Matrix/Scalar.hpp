/**
 * @file Scalar.hpp
 *
 * Defines conversion of matrix to scalar.
 *
 * @author James Goppert <james.goppert@gmail.com>
 */

#ifndef __SCALAR_HPP__
#define __SCALAR_HPP__

#include "math.hpp"

namespace matrix
{

template<typename Type>
class Scalar
{
public:
    Scalar() : _value()
    {
    }

    Scalar(const Matrix<Type, 1, 1> & other)
    {
        _value = other(0,0);
    }

    Scalar(Type other)
    {
        _value = other;
    }

    operator Type &()
    {
        return _value;
    }

    operator Matrix<Type, 1, 1>() const {
        Matrix<Type, 1, 1> m;
        m(0, 0) = _value;
        return m;
    }

    operator Vector<Type, 1>() const {
        Vector<Type, 1> m;
        m(0) = _value;
        return m;
    }

private:
    Type _value;

};

typedef Scalar<float> Scalarf;

} // namespace matrix
#endif
/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
