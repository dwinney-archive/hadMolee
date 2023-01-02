// Methods for printing messages out to the command-line easily
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef PRINT_HPP
#define PRINT_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>

namespace hadMolee
{
    // Full width is the length of how wide I want the output to be
    static const int TEXTWIDTH = 62;

    // Default spacing between multiple columns
    const int PRINT_SPACING = 15;

    // Default precisions when converting non-ints to strings
    const int STRING_PRECISION = 3;

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

    // ---------------------------------------------------------------------------
    // Non-error printing

     // Output an empty line to the terminal
    inline void empty_line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    inline void divider()
    {
        std::cout << "--------------------------------------------------------------" << std::endl;
    };

    inline void dashed_divider()
    {
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    };

    template<typename T>
    inline void print(T x)
    {
        std::cout << std::boolalpha << std::left;  
        std::cout << std::setw(PRINT_SPACING) << x << std::endl;
    };

    template <typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left; 
        std::cout << std::setw(PRINT_SPACING) << first;
        print(rest...);
    } 

    // Print out a vector element-wise
    template<typename T>
    inline void print(std::vector<T> v)
    {
        std::cout << std::boolalpha; 
        for (auto vi : v)
        {
            std::cout << std::left << std::setw(PRINT_SPACING) << vi;
        };
        std::cout << std::endl;
    };

    // Print a string centered on the terminal 
    inline void centered(std::string words)
    {
        int x = words.length();
        int gap_width = (TEXTWIDTH - x)/2;
        std::cout << std::left << std::setw(gap_width) << "" << std::setw(x) << words << std::setw(gap_width) << "" << std::endl;
    };

    // ---------------------------------------------------------------------------
    // String operations    

    // Produce a string with the format "name = value units"

    template <typename T>
    inline std::string var_def(std::string name, T value, std::string units = "")
    {
        std::stringstream ss;
        ss << std::setprecision(STRING_PRECISION) << name + " = " << value << " " + units;
        return ss.str();
    };
};

#endif