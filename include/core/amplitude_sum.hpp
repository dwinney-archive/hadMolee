// Class which allows any number of amplitude sharing the same kinematics to be added together
// the result is treated as an amplitude itself
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef SUM
#define SUM

#include "amplitude.hpp"
#include "charmoniumlike.hpp"
#include "reaction_kinematics.hpp"

namespace hadMolee
{
    class amplitude_sum : public amplitude
    {        
        // -----------------------------------------------------------------------
        public:

        // Basic constructor like parent class
        amplitude_sum(reaction_kinematics * xkinem, charmoniumlike * V, std::string identifer = "amplitude_sum")
        : amplitude(xkinem, V, 0, "amplitude_sum", identifer)
        {};

            // Constructor with a vector already set up
        amplitude_sum(reaction_kinematics * xkinem,  charmoniumlike * V, std::vector<amplitude*> vec, std::string identifer = "amplitude_sum")
        : amplitude(xkinem, V, 0, "amplitude_sum", identifer)
        {
            // Add the given vector to the store pointers 
            add(vec);
        };

        // Add a new amplitude to the vector
        inline void add(amplitude * new_amp)
        {
            if (check_compatibility(new_amp))
            {
                _amps.push_back(new_amp);
                set_nParams(N_parameters() + new_amp->N_parameters());
            };
        };

        // Add a whole vector of new amplitudes
        inline void add(std::vector<amplitude*> new_amps)
        {
            for (amplitude* amp : new_amps)
            {
                add(amp);
            };
        };

        // Allocate an aggragated vector of parameters to individual amplitudes
        inline void set_parameters(std::vector<double> x)
        {
            check_nParams(x);
        
            int N = 0;
            // Allocate parameters to the individual amplitudes
            for (amplitude* amp : _amps)
            {
                // Extract the subvector corresponding to the i-th amplitude
                auto start = x.begin() + N;
                auto end   = x.begin() + N + amp->N_parameters();
                std::vector<double> pars(start, end);

                amp->set_parameters(pars);
                N += amp->N_parameters();
            };

            // At the end check that the number of params allocated is the same as expected
            // However is this happpens the damage is done so we just send out a warning...
            if (N != N_parameters())
            {
                std::cout << "Warning! amplitude_sum::set_params() : Number of parameters allocated doesnt match those expected..." << std::endl;
            };
        };

        std::complex<double> reduced_amplitude(cartesian_index i, cartesian_index j);

        // -----------------------------------------------------------------------
        private:

        // Sums hold a vector of pointers to the constituent amplitudes
        std::vector<amplitude*> _amps;

        // Make sure all amplitudes being summed are compatible
        bool check_compatibility(amplitude* amp);
    };
};

#endif