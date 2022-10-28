// Useful utilities for debugging and implementing error messages
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DEBUG
#define DEBUG

#include <iostream>
#include <iomanip>

namespace hadMolee 
{
    // ---------------------------------------------------------------------------
    // DEBUGGING Messages

    // Default spacing value
    const int DEBUG_SPACING = 15;

    // Functions for printing to screen instead of having to copy this line all the time
    template<typename T>
    inline void debug(T x)
    {
        std::cout << x << std::endl;
    };

    template<typename T, typename F>
    inline void debug(T x, F y)
    {
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y << std::endl;
    };

    template<typename T, typename F, typename G>
    inline void debug(T x, F y, G z)
    {
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y;
        std::cout << std::left << std::setw(DEBUG_SPACING) << z << std::endl;
    };

    template<typename T, typename F, typename G, typename H>
    inline void debug(T x, F y, G z, H a)
    {
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y;
        std::cout << std::left << std::setw(DEBUG_SPACING) << z;
        std::cout << std::left << std::setw(DEBUG_SPACING) << a << std::endl;
    };

    template<typename T, typename F, typename G, typename H, typename I>
    inline void debug(T x, F y, G z, H a, I b)
    {
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y;
        std::cout << std::left << std::setw(DEBUG_SPACING) << z;
        std::cout << std::left << std::setw(DEBUG_SPACING) << a;
        std::cout << std::left << std::setw(DEBUG_SPACING) << b << std::endl;
    };

    // ---------------------------------------------------------------------------
    // ERROR Messages

    // Throw an error message then quits code 
    inline void fatal()
    {
        std::cout << std::left << "FATAL ERROR! Quiting..." << std::endl;
        exit( EXIT_FAILURE );
    };

    // Error message with location and reason messages too
    inline void fatal(std::string location, std::string reason = "")
    {
        std::cout << std::left << "FATAL ERROR! " + location + ": " + reason << std::endl;
        std::cout << std::left << "Quiting..." << std::endl;

        exit( EXIT_FAILURE );
    };

    // Warning message does not exit code but throws a message up
    inline void warning(std::string location, std::string message = "")
    {
        std::cout << std::left << "WARNING! " + location + ": " + message << std::endl;
    };
};

#endif
