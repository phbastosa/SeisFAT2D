# ifndef KIRCHHOFF_HPP
# define KIRCHHOFF_HPP

# include "migration.hpp"
# include <cuda_runtime.h>

class Kirchhoff : public Migration
{
private:

    int nBlocks;
    int nThreads;

    float * d_Tr = nullptr;
    float * d_Ts = nullptr;
    float * d_image = nullptr;
    float * d_seismic = nullptr;

    void set_components();
    void run_cross_correlation();
};

__global__ void cross_correlation(float * seismic, float * Ts, float * Tr, float * image, int nPoints, int spread, int nt, float dt);

# endif