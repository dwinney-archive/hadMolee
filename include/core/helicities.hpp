// Methods involving helicty combinations
// These are hard-coded in arrays so that they can easily accessed with a single index
//
// Author:       Daniel Winney (2022)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#ifndef HELIC_COMBO
#define HELIC_COMBO

#include "debug.hpp"

using namespace std;

#include <algorithm>
#include <iostream>
#include <vector>
#include <array>

// Assuming we have a V -> PPV final state we have either 6 or 9 helicity combinations
// Depending on if the initial state particle can have a longitudinal polarization

const vector<array<int,2>> MASSLESS_HELICITIES =
{
//  { V,  a}, // i
    { 1,  1}, // 0
    { 1,  0}, // 1
    { 1, -1}, // 2
    {-1,  1}, // 3
    {-1,  0}, // 4
    {-1, -1}, // 5
};

const vector<array<int,2>> MASSIVE_HELICITIES =
{
//  { V,  a}, // i
    { 1,  1}, // 0
    { 1,  0}, // 1
    { 1, -1}, // 2
    { 0,  1}, // 3
    { 0,  0}, // 4
    { 0, -1}, // 5
    {-1,  1}, // 6
    {-1,  0}, // 7
    {-1, -1}, // 8
};

// Access the above vectors by checking if we want a photon or not
inline vector<array<int,2>> get_helicities(bool is_photon = false)
{
    if (is_photon) return MASSLESS_HELICITIES;
    else           return MASSIVE_HELICITIES;
};
inline array<int,2> get_helicities(int n, bool is_photon = false)
{
    if (is_photon) return MASSLESS_HELICITIES[n];
    else           return MASSIVE_HELICITIES[n];
};

// Print out a string displaying a given set of helicities
inline string print_helicities(array<int,2> lam)
{
    array<string,2> lams;
    for (int i = 0; i < 2; i++)
    {
        switch (lam[i])
        {
            case  1: lams[i] = "#plus";    break;
            case  0: lams[i] = "0";        break;
            case -1: lams[i] = "#minus";   break;
            default: continue;
        };
    };
    string hels = "{" + lams[0] + "," + lams[1] + "}";
    return hels;  
};
// Same thing but with a an index
inline string print_helicities(int n, bool is_photon = true)
{
    array<int,2> hels = get_helicities(n, is_photon);
    return print_helicities(hels);
};

// Input an int get out the helicity combination corresponding to that index
inline int find_index( array<int,2> helicities, bool is_photon)
{
    vector<array<int,2>> hels;
    (is_photon) ? (hels = MASSLESS_HELICITIES) : (hels = MASSIVE_HELICITIES );

    auto iterator = find(hels.begin(), hels.end(), helicities);
    if (iterator != hels.end())
    {
        return iterator - hels.begin();
    }
    
    warning("find_index", "Cannot find index for " + print_helicities(helicities) + ". Returning -1...");
    return -1;
};

#endif