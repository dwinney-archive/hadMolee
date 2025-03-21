// Class and methods for handling data sets used for fitting
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef DATA_SET_HPP
#define DATA_SET_HPP

#include "constants.hpp"
#include "elementwise.hpp"

#include <fstream>
#include <sstream>

namespace hadMolee
{

    // Importing data sets we'll need to be able to find the /data/ directory from the 
    // top level one. Thus we need to be able to access the HADMOLEE environment variable
    inline std::string hadMolee_dir()
    {
       // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("HADMOLEE");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("import_data", "Cannot find environment variable HADMOLEE!", "");
        }
        return std::string(env);  
    };

    // Import a set of data with N columns with relative path 
    // and full path hadMolee_dir/ + rel_path
    template<int N> 
    inline std::array<std::vector<double>,N> import_data(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = hadMolee_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data", "Cannot open file " + file_path + "!", result);
        };

        // Import data!
        std::string line;
        while (std::getline(infile, line))
        {   
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            for (int i = 0; i < N; i++)
            {
                double x;
                is >> x;
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // Similar to above except that the data is transposed, i.e. rows are the "categories"
    // and the columns are data points. We specify the number of rows in this case
    template<int N>
    inline std::array<std::vector<double>, N> import_transposed(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = hadMolee_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data", "Cannot open file " + file_path + "!", result);
        };

        // Import data!
        for (int i = 0; i < N; i++)
        {   
            std::string line;
            std::getline(infile, line);
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            double x;
            while(is >> x)
            {
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // If data file has more columns than are actually needed,
    // import with import_data and use this to throw out all but the desired columns
    template<int Nin, int Nout> 
    inline std::array<std::vector<double>,Nout> reshape_data(std::array<std::vector<double>,Nin> data, std::array<int,Nout> to_keep)
    {
        std::array<std::vector<double>, Nout> result;

        for (int i = 0; i < Nout; i++)
        {
            result[i] = data[ to_keep[i] ];
        };
        return result;
    };

    // Make sure all the vectors are the correct size
    template<int S>
    inline int check(std::array<std::vector<double>,S> data, std::string id)
    {
        // Grab the size of the first entry
        int N = data[0].size();
        
        // And compare to the rest
        for (auto column : data)
        {
            if (column.size() != N)
            {
                warning("data_set", "Input vectors of " + id + " have mismatching sizes!");
                return 0;
            };
        };

        return N;
    };

    // Define different archetypes of data
    enum data_type { integrated_data, differential_data};

    struct data_set
    {
        // Number of data points
        int _N = 0;

        std::string _id = "data_set";

        data_type _type;
        
        // Vectors to store energy and momentum transfer variables and observable
        std::vector<double> _w, _t, _obs, _obserr;

        // Other possible vectors to store things like bin sizes, etc
        std::array<std::vector<double>,2> _werr, _terr;

        // Whether the values stored in _w correspond to invariant energy W = sqrt(s) (false)
        // or lab frame energy Egamma (true)
        bool _lab = false;

        // Whether the momentum transfer values stored in _t correspond to invariant t (false)
        // or t' = t - t_min (true)
        bool _tprime = false;

        // Whether saves values in _t are positive or negative t
        // i.e. -t (true) vs t (false)
        bool _negt   = false;

        // For a differential set it may be useful to have an average s 
        double _avg_w = 0;

        // If we want a data entry in t he legend
        bool _add_to_legend = false;
    };
};

#endif