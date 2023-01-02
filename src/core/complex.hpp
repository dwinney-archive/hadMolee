// Extensions of the std::complex library for easier use throughout 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <complex>

namespace hadMolee
{
    // We use complex numbers a lot so we define this shortened data type
    using complex = std::complex<double>;

    // Unit imaginary and real
    const complex I  (0., 1.);
    const complex R  (1., 0.);

    // ---------------------------------------------------------------------------
    // Overload of operators for a int/bool and complex<double>

    inline complex operator * (const bool & a, const complex & b)
    {
        if (a)
            return b;
        else
            return complex(0, 0);
    };
    inline complex operator * (const complex & b, const bool & a)
    {
        if (a)
            return b;
        else
            return complex(0, 0);
    };

    // Same thing but with int
    inline complex operator * (const int & a, const complex & b)
    {
        if (a != 0) return complex(double(a) * real(b), double(a) * imag(b));
        else        return 0.;
    };
    inline complex operator * (const complex & b, const int & a)
    {
        if (a != 0) return complex(double(a) * real(b), double(a) * imag(b));
        else        return 0.;
    };

    // Also allow division
    inline complex operator/(const complex&c, const int& z)
    {
        return (1./z)*c;
    };

    inline complex operator+(const complex&c, const int& z)
    {
        return c + R*z;
    };

    inline complex operator+(const int& z, const complex & c)
    {
        return R*z + c;
    };

    inline complex operator-(const complex&c, const int& z)
    {
        return c - R*z;
    };

    inline complex operator-(const int& z, const complex & c)
    {
        return R*z - c;
    };

    // This makes it so we always default to complex regardless of whether the input is an int or double
    template<typename T>
    complex csqrt(T x){ return sqrt(x * R); };
};

#endif