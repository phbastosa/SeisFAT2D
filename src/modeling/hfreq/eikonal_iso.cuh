# ifndef EIKONAL_ISO_CUH
# define EIKONAL_ISO_CUH

# include "eikonal.hpp"

class Eikonal_Iso : public Eikonal
{
private:

    int total_levels;
    int nSweeps, meshDim;
    int nThreads, nBlocks;

    float * d_T = nullptr;
    float * d_S = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    void set_properties();
    void set_conditions();
    
public:

    void forward_solver();
};

__global__ void inner_sweep(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz, float dx2i, float dz2i);

# endif
