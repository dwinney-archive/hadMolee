// Abstract class for an amplitude of a particles scalar propagator
// as well as some implementations for the Y and Z mesons
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef BREITWIGNER
#define BREITWIGNER

#include "constants.hpp"

namespace hadMolee
{

    // First the abstract class that can be called from anywhere
    class breit_wigner
    {
        public: 

        // Default constructor only requires a 'bare' mass 
        breit_wigner(int mode)
        : _mode(mode)
        {};

        breit_wigner(int mode, double mass, double width = 0.)
        : _mode(mode), _mass(mass), _width(width)
        {};

        // Main function, evaluates as a function of one variable
        complex eval(double x)
        {
            switch (_mode)
            {
                case (kRelativistic):    return rel_eval(x);
                case (kNonrelativistic): return nonrel_eval(x);
                default:                 return std::nan("");
            }
        };
        inline double squared(double x){ return norm(eval(x)); };

        // Set the mass and width of the state
        inline void set_mass(double m, double w = 0){ _mass = m; _width = w; };

        // Different options for evaluation
        static const int kRelativistic    = 0;
        static const int kNonrelativistic = 1;

        protected:

        // Whether we are using the relativistic BW or not
        int _mode;

        complex nonrel_eval(double E)
        {
            complex D =  (E - _mass) + I * _width/2. + IEPS;
            return 1. / (2. * _mass * D);
        };  

        complex rel_eval(double s)
        {
            complex D =  (s - _mass*_mass) + I * _mass*_width + IEPS;
            return 1. / D;
        };

        // Pole mass is only really needed parameter 
        double _mass;

        // But may also allow a width
        double _width;
    };
};

#endif