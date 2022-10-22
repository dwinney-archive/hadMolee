// Utility methods for converting digitized data files to nicer formats for use in fits
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include <fstream>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Take in three files for the central, and error values and reformat them into a single file
// The final argument bool is whether or not to round down to "nearest event"
void reformat_digitized(string central, string lower_error, string upper_error, string output, bool round = false);


// ---------------------------------------------------------------------------
// Take in path to a .dat file and import it as an array of vectors
// Output will be an array of vectors containing each column
std::array<std::vector<double>,4> import_data(string filename);