# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../admin/admin.hpp"

class Geometry
{
public:

    int nsrc;
    int nrec;

    float * xsrc = nullptr;
    float * zsrc = nullptr;

    float * xrec = nullptr;
    float * zrec = nullptr;

    std::string parameters;

    void set_parameters();     
};

# endif
