# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../ioFunctions/ioFunctions.hpp"

class Geometry
{
public:

    int nsrc;
    int nrec;
    int nrel;

    int * sInd = nullptr;
    int * iRec = nullptr;
    int * fRec = nullptr;

    int * spread = nullptr;

    float * xsrc = nullptr;
    float * zsrc = nullptr;

    float * xrec = nullptr;
    float * zrec = nullptr;

    std::string parameters;

    void set_parameters();     
};

# endif
