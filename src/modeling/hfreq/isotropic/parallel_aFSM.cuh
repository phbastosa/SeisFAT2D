# ifndef PARALLEL_AFSM_CUH
# define PARALLEL_AFSM_CUH

# include "../eikonal.hpp"

# include <cuda_runtime.h>

class Parallel_aFSM : public Eikonal
{
private:

    // int nSweeps;
    // int meshDim;

    // int total_levels;
    // int blocksPerGrid;
    // int threadsPerBlock;

    // float * S = nullptr;
    // float * T = nullptr;

    // int * d_sgnv = nullptr;
    // int * d_sgnt = nullptr;

    // void set_specifications();

protected:

public:
    
    void set_name();
    // void initialization();
    // void forward_solver();
};

// __global__ void kernel_FSM(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz, float dx2i, float dz2i);

# endif