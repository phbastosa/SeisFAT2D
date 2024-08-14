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

    int totalLevels;
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

# endif