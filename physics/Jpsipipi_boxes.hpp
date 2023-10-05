// Specific implementations of amplitudes relevant for the Jpsi pi pi final state

// Author:       Daniel Winney (2023)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef JPSIPIPI_BOXES_HPP
#define JPSIPIPI_BOXES_HPP

#include <memory>

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "box.hpp"
#include "Y(4260).hpp"

namespace hadMolee
{
    class Jpsipipi_BoxI : public amplitude_base
    {
        public:

        Jpsipipi_BoxI(amplitude_key key, kinematics xkinem, lineshape Y, std::string id = "Box I")
        : amplitude_base(key, xkinem, Y, 0, "Jpsipipi_BoxI", id), 
          _B(box::kLoopTools),
          _Y(get_molecular_component(Y))
        {
            _B.set_internal_masses(_internal);
            _B.set_external_masses(_external);
            _B.add_width(0, W_D1);
        };

        // The reduced amplitude is a D-wave amplitude
        inline complex reduced_amplitude(cartesian_index i, cartesian_index j)
        {
            return _C/* (3.*phat(i)*phat(j) - delta(i,j)) */;
        };

        private:

        inline void recalculate()
        {
            // Y mass and couplings
            double y   = _Y->molecular_coupling();
            double M_Y = sqrt(_s);

            // Update floating masses in the box and evaluate
            _B.set_external_masses({_W, M_PION, M_JPSI, M_PION});
            _B.set_internal_masses( _sab, _sbc );
            _vecB = _B->eval_vector();
        };

        // Scalar prefactors 
        complex _C = 0.;  

        // Box functions
        // Fixed masses for all the particles in the box
        box _B;
        std::array<double, 4> _internal = {M_D1, M_D, M_DSTAR, M_DSTAR}; 
        std::array<double, 4> _external = {4.2, M_PION, M_JPSI, M_PION};

        // Array to store the vector components of the box integral
        std::array<complex,4> _vecB;

        // We need access to the D1D molecular nature of the Y state
        molecule _Y;
    };

};

#endif