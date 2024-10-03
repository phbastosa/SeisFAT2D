# ifndef MODELING_HPP
# define MODELING_HPP

# include "../geometry/geometry.hpp"

class Modeling
{
private:

protected:

    int spread, max_spread;

    virtual void set_specifications() = 0;

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

public:

    float dx, dz, dt;
    int nxx, nzz, matsize;
    int nt, nx, nz, nb, nPoints;
    int srcId, recId, sIdx, sIdz;

    float * S = nullptr;
    float * T = nullptr;
    float * Vp = nullptr;
    float * Vs = nullptr;
    float * Rho = nullptr;

    float * synthetic_data = nullptr;

    Geometry * geometry;

    std::string parameters;

    void set_parameters();
    void show_information();    

    virtual void initialization() = 0;
    virtual void forward_solver() = 0;
};

# endif