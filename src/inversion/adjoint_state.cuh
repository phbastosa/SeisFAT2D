# ifndef ADJOINT_STATE_HPP
# define ADJOINT_STATE_HPP

# include "tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    int total_levels;
    int nSweeps, meshDim;
    int nThreads, nBlocks;

    float cell_area;

    float * m = nullptr;
    float * v = nullptr;

    float * d_T = nullptr;

    float * d_source = nullptr;
    float * d_adjoint = nullptr;

    float * source = nullptr;
    float * adjoint = nullptr;
    float * gradient = nullptr;

    void initialization();
    void set_specifications();

    void apply_inversion_technique();

public:

    void optimization();

};

__global__ void inner_sweep(float * T, float * adjoint, float * source, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz); 

# endif
