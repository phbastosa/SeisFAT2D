# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../geometry/geometry.hpp"

class Eikonal
{
protected:

    virtual void set_specifications() = 0;

public:

    float dx, dz;
    int srcId, recId;
    int nxx, nzz, matsize;
    int nx, nz, nb, nPoints;

    float * velocity = nullptr;
    float * slowness = nullptr;
    float * eikonalT = nullptr;

    float * synthetic_data = nullptr;

    std::string data_folder;

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

    Geometry * geometry;

    std::string parameters;
    
    void set_parameters();
    void show_information();    
    void get_synthetic_data();

    virtual void initialization() = 0;
    virtual void forward_solver() = 0;
};

# endif