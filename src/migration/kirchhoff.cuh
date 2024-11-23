# ifndef KIRCHHOFF_HPP
# define KIRCHHOFF_HPP

# include "migration.hpp"
# include <cuda_runtime.h>

# define PI 3.14159265359

class Kirchhoff : public Migration
{
private:

    int nBlocks;
    int nThreads;

    float * d_Tr = nullptr;
    float * d_Ts = nullptr;
    float * d_image = nullptr;
    float * d_seismic = nullptr;

    void set_specifications();
    void run_cross_correlation();
};

__global__ void cross_correlation(float * seismic, float * Ts, float * Tr, float * image, float aperture, float cmp, int nPoints, int spread, int nz, int nt, float dt, float dx, float dz);

# endif