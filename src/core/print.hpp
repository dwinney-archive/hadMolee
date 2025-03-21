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
    const std::string UNIT_DIV = std::string(PRINT_SPACING, '-');

    // Default precisions when converting non-ints to strings
    const int STRING_PRECISION = 3;

    // ---------------------------------------------------------------------------
    // Non-error printing

     // Output an empty line to the terminal
    inline void empty_line()
    {
        std::cout << std::endl;
    };

     // Output an empty line to the terminal
    inline void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    inline void divider()
    {
        std::cout << "--------------------------------------------------------------" << std::endl;
    };

    inline void divider(int n)
    {
        std::string div;
        for (int i = 0; i < n; i++)
        {
            div = div + UNIT_DIV;
        }
        std::cout << div << std::endl;
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