// Utility methods for converting digitized data files to nicer formats for use in fits
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@scnu.edu.cn
// ---------------------------------------------------------------------------

#include "data_formatting.hpp"

// ---------------------------------------------------------------------------
// Take in three files for the central, and error values and reformat them into a single file
// The final argument bool is whether or not to round down to "nearest event"
void reformat_digitized(string central_value, string lower_error, string upper_error, string output, bool rounded)
{
    // Open the three files and import them into vectors 
    std::vector<double> E, central, upper, lower;

    // Open the three files
    std::ifstream central_file(central_value.c_str());
    if (central_file.fail())
    {
        warning("reformat_digitized", "Couldn't open file " + central_value + " ! Continuing without reformatting...");
        return;
    }

    std::ifstream lower_file(lower_error.c_str());
    if (lower_file.fail())
    {
        warning("reformat_digitized", "Couldn't open file " + lower_error + " ! Continuing without reformatting...");
        return;
    }

     std::ifstream upper_file(upper_error.c_str());
    if (upper_file.fail())
    {
        warning("reformat_digitized", "Couldn't open file " + upper_error + " ! Continuing without reformatting...");
        return;
    }
    
    // Import all the data
    double x, y;
    while (central_file >> x >> y)
    {
        E.push_back(x); 
        central.push_back(y);
    }
    while (lower_file >> x >> y)
    {
        lower.push_back(y); 
    }
    while (upper_file >> x >> y)
    {
        upper.push_back(y); 
    }

    // Close everything up
    central_file.close();
    upper_file.close();
    lower_file.close();

    // Check theyre all the correct size
    if ((E.size() != upper.size()) || (E.size() != lower.size()))
    {
        warning("reformat_digitized", "Input files sizes dont match! Continuing without reformatting...");
        return;
    };

    // Then output to file
    std::ofstream output_file;
    output_file.open(output.c_str());

    if (!rounded)
    {
        for (int i = 0; i < E.size(); i++)
        {
            output_file << left;
            output_file << setw(15) << E[i];
            output_file << setw(15) << central[i];
            output_file << setw(15) << central[i] - lower[i];
            output_file << setw(15) << upper[i]   - central[i] << endl;
        }
    }
    else
    {
        for (int i = 0; i < E.size(); i++)
        {
            output_file << left;
            output_file << setw(15) << E[i];
            output_file << setw(15) << abs(round(central[i]));
            output_file << setw(15) << abs(round(central[i] - lower[i]));
            output_file << setw(15) << abs(round(upper[i] - central[i])) << endl;
        }
    }
    output_file.close();


    return;
};

// ---------------------------------------------------------------------------
// Take in path to a .dat file and import it as an array of vectors
// Output will be an array of vectors containing each column
std::array<std::vector<double>,4> import_data(string filename)
{
    // Output will contain column data from here
    std::vector<double> E, central, upper, lower;

    // Open the file
    std::ifstream file(filename.c_str());
    if (file.fail())
    {
        warning("import_data", "Couldn't open file " + filename + " ! Continuing without reformatting...");
        return {{}};
    }

    double e, cen, up, low;
    while (file >> e >> cen >> up >> low)
    {
        E.push_back(e); 
        central.push_back(cen);
        upper.push_back(up);
        lower.push_back(low);
    }
    file.close();

    // Check theyre all the correct size
    if ((E.size() != central.size()) || (E.size() != upper.size()) || (E.size() != lower.size()))
    {
        warning("import_data", "Input files sizes dont match! Continuing without reformatting...");
        return {{}};
    };

    return {E, central, upper, lower};
};