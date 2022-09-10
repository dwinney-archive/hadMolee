// Useful utilities for debugging and implementing error messages
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef DEBUG
#define DEBUG

#include <iostream>
#include <iomanip>

using namespace std;

// ---------------------------------------------------------------------------
// DEBUGGING Messages

// Default spacing value
const int DEBUG_SPACING = 15;

// Functions for printing to screen instead of having to copy this line all the time
template<typename T>
inline void debug(T x)
{
    cout << x << endl;
};

template<typename T, typename F>
inline void debug(T x, F y)
{
    cout << left << setw(DEBUG_SPACING) << x;
    cout << left << setw(DEBUG_SPACING) << y << endl;
};

template<typename T, typename F, typename G>
inline void debug(T x, F y, G z)
{
    cout << left << setw(DEBUG_SPACING) << x;
    cout << left << setw(DEBUG_SPACING) << y;
    cout << left << setw(DEBUG_SPACING) << z << endl;
};

template<typename T, typename F, typename G, typename H>
inline void debug(T x, F y, G z, H a)
{
    cout << left << setw(DEBUG_SPACING) << x;
    cout << left << setw(DEBUG_SPACING) << y;
    cout << left << setw(DEBUG_SPACING) << z;
    cout << left << setw(DEBUG_SPACING) << a << endl;
};

// ---------------------------------------------------------------------------
// ERROR Messages

// Throw an error message then quits code 
inline void fatal()
{
    cout << left << "FATAL ERROR! Quiting..." << endl;
    exit( EXIT_FAILURE );
};

// Error message with location and reason messages too
inline void fatal(string location, string reason = "")
{
    cout << left << "FATAL ERROR! " + location + ": " + reason << endl;
    cout << left << "Quiting..." << endl;

    exit( EXIT_FAILURE );
};

// Warning message does not exit code but throws a message up
inline void warning(string location, string message = "")
{
    cout << left << "WARNING! " + location + ": " + message << endl;
};


#endif
