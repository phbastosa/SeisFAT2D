# ifndef IDLSKDM_CUH
# define IDLSKDM_CUH

# include "migration.cuh"

class IDLSKDM : public Migration
{
    void forward();
    void adjoint();

    void image_building();
    void export_outputs();
};

__global__ void forward_kernel(float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz);
__global__ void adjoint_kernel(float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz);

# endif