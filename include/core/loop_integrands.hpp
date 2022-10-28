// Auxillary classes used to evaluate the Feynman parameter integrands for loop integrals
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef INTEGRANDS
#define INTEGRANDS

#include "constants.hpp"

namespace hadMolee
{
    // Auxillary class to help interfacing with the cubature library
    class triangle_integrand
    {
        // -------------------------------------------------------------------

        public:

        // Empty constructor cus we dont need anything
        triangle_integrand(){};

        // Get masses squared from the parent triangle function
        inline void update_masses(std::array<double,3> ex_m2, std::array<double,3> in_m2, std::array<double,3> in_w)
        {
            _ema2 = ex_m2[0]; _emb2 = ex_m2[1]; _emc2 = ex_m2[2];
            _ima2 = in_m2[0]; _imb2 = in_m2[1]; _imc2 = in_m2[2];
            _wa   =  in_w[0]; _wb   =  in_w[1]; _wc   =  in_w[2];
        };

        // Evaluate the integrand at fixed values of the feynman parameters
        inline std::complex<double> eval(double x1, double x2, double x3)
        {
            std::complex<double> D = x1*(_ima2 - XI*sqrt(_ima2)*_wa) 
                            + x2*(_imb2 - XI*sqrt(_imb2)*_wb) 
                            + x3*(_imc2 - XI*sqrt(_imc2)*_wc) 
                            - x1*x3*_emb2
                            - x1*x2*_emc2 
                            - x2*x3*_ema2;
            return 1. / (D - XI*_eps);
        };

        // Set the numerical iepsilon used
        inline void set_ieps(double e){ _eps = e; };

        // -------------------------------------------------------------------

        private:

        // Default epsilon
        double _eps = EPS;

        // Need to be able to access the masses
        double _ema2, _emb2, _emc2;
        double _ima2, _imb2, _imc2;
        double _wa, _wb, _wc;
    };

    // Auxillary class to help interfacing with the cubature library
    class box_integrand
    {
        // -------------------------------------------------------------------

        public:

        // Empty constructor cus we dont need anything
        box_integrand(){};

        // Get masses squared from the parent triangle function
        inline void update_masses(std::array<double,2> st, std::array<double,4> ex_m2, std::array<double,4> in_m2, std::array<double,4> in_w)
        {
            _p01 = ex_m2[0]; _p12 = ex_m2[1]; _p23 = ex_m2[2]; _p03 = ex_m2[3];
            _m0  = in_m2[0]; _m1  = in_m2[1]; _m2  = in_m2[2]; _m3  = in_m2[3];
            _w0  =  in_w[0]; _w1  =  in_w[1]; _w2  =  in_w[2]; _w3  =  in_w[3];

            _p02    = st[0];
            _p13    = st[1];
        };

        // Evaluate the integrand at fixed values of the feynman parameters
        inline std::complex<double> eval(double x0, double x1, double x2, double x3)
        {
            std::complex<double> Y01, Y02, Y03, Y12, Y13, Y23;

            std::complex<double> M0 = _m0 - XI*sqrt(_m0*XR)*_w0;
            std::complex<double> M1 = _m1 - XI*sqrt(_m1*XR)*_w1;
            std::complex<double> M2 = _m2 - XI*sqrt(_m2*XR)*_w2;
            std::complex<double> M3 = _m3 - XI*sqrt(_m3*XR)*_w3;

            Y01 = M0 + M1 - _p01;
            Y02 = M0 + M2 - _p02;
            Y03 = M0 + M3 - _p03;
            Y12 = M1 + M2 - _p12;
            Y13 = M1 + M3 - _p13;
            Y23 = M2 + M3 - _p23;        

            std::complex<double> D = x0*x0*M0 + x1*x1*M1 + x2*x2*M2 + x3*x3*M3
                            + x0*x1*Y01  + x0*x2*Y02 + x0*x3*Y03
                                        + x1*x2*Y12 + x1*x3*Y13 
                                                    + x2*x3*Y23;
            return pow(D - XI*_eps, -2.);
        };

        // Set the numerical iepsilon used
        inline void set_ieps(double e){ _eps = e; };

        // -------------------------------------------------------------------

        private:

        // Default epsilon
        double _eps = EPS;

        // Need to be able to access the masses
        double _p01, _p02, _p03, _p12, _p13, _p23; 
        double _m0, _m1, _m2, _m3;
        double _w0, _w1, _w2, _w3;
    };
};

#endif