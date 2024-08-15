# ifndef EIKONAL_CUH
# define EIKONAL_CUH

# include "modeling.hpp"

# include <cuda_runtime.h>

class Eikonal : public Modeling
{
private:

    int nSweeps;
    int meshDim;

    float dz2i, dx2i;

    int total_levels;

    float * S = nullptr;
    float * T = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    void set_eikonal_parameters();

public:

    void initialization();
    void forward_solver();

};

__global__ void fast_sweeping_method(int z_offset, int zd, int x_offset, int xd, int nxx, int nzz);

# endif