# ifndef MODELING_HPP
# define MODELING_HPP

# include "../geometry/geometry.hpp"

class Modeling
{
private:

    void check_geometry_bounds();

protected:

    virtual void set_eikonal_parameters() = 0;

public:

    float dx, dz;
    int shot_index;
    int sidx, sidy, sidz;
    int nxx, nzz, matsize;
    int nx, nz, nb, nPoints;

    float * velocity = nullptr;
    float * slowness = nullptr;
    float * eikonalT = nullptr;

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

    Geometry * geometry;

    std::string parameters;
    
    void set_parameters();
    
    virtual void initialization() = 0;
    virtual void forward_solver() = 0;

};

// void fast_sweeping_method_CPU();
// void inner_sweep(int i, int j, int sgnvx, int sgnvz, int sgntx, int sgntz, float dx2i, float dz2i);

# endif