// Abstract class for an amplitude of a particles scalar propagator
// as well as some implementations for the Y and Z mesons
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "propagator.hpp"

complex<double> hadronic_propagator::eval(double E)
{
    // The self-energy corrections contains two parts:
    // from the m1 m2 -> X -> m1 m2 resonant reactions 
    // and from the elastic scattering
    auto sigma = [&] (double e)
    {
        return  self_energy(E) * (_binding*_binding - 2.*(E-_M)*_elastic);
    };

    // Calculate the partial derivative at E = M
    auto Resigma = [&] (double e)
    {
        return real(sigma(e));
    };
    double derivative = boost::math::differentiation::finite_difference_derivative(Resigma, _M);
    complex<double> sigma_tilde = sigma(E) - real(sigma(_M)) - (E - _M) * derivative;

    // Propagator dressed by constant width and self energy above
    complex<double> D = E - _M - sigma_tilde + XI*_nonmol_width / 2.;

    return XI / (2. * D);
};

complex<double> hadronic_propagator::self_energy(double E)
{
    double mu  = reduced_mass();
    double eps = E - _m1 - _m2; 
    double prefactors = pow(mu, 3./2.) / (8.*PI);

    return - prefactors * sqrt(2.*abs(eps)) * ( XI*(eps >= 0) + XR*(eps < 0) );
};
