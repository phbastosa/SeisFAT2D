# ifndef EIKONAL_CUH
# define EIKONAL_CUH

# include "modeling.hpp"

# include <cuda_runtime.h>

class Eikonal : public Modeling
{
private:

    int nSweeps;
    int meshDim;

    int total_levels;
    int blocksPerGrid;
    int threadsPerBlock;

    float * S = nullptr;
    float * T = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    void set_eikonal_parameters();

public:

    void initialization();
    void forward_solver();

};

__global__ void fast_sweeping_method(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz, float dx2i, float dz2i);

# endif