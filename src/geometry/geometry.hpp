# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../ioFunctions/ioFunctions.hpp"

class Geometry
{
private:

    bool read_geometry;

    std::string xps_file;
    std::string rps_file;
    std::string sps_file;

    std::vector<float> linspace(float xi, float xf, int n);

public:

    int nsrc;
    int nrec;
    int nrel;

    int * sInd = nullptr;
    int * iRec = nullptr;
    int * fRec = nullptr;

    float * xsrc = nullptr;
    float * zsrc = nullptr;
    float * xrec = nullptr;
    float * zrec = nullptr;

    Geometry(std::string parameters);     
};

# endif