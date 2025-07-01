# ifndef MODELING_HPP
# define MODELING_HPP

# include <cuda_runtime.h>

# include "../geometry/geometry.hpp"

# define NSWEEPS 4
# define MESHDIM 2

# define COMPRESS 65535

typedef unsigned short int uintc; 

class Modeling
{
private:

    void set_eikonal();
    void set_properties();

    float cubic1d(float P[4], float dx);
    float cubic2d(float P[4][4], float dx, float dz);

protected:

    int total_levels;
    int nThreads, nBlocks;

    float dx2i, dz2i;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    virtual void set_conditions() = 0;
    
    void compression(float * input, uintc * output, int N, float &max_value, float &min_value);

public:

    float dx, dz, sx, sz;
    int nxx, nzz, matsize;
    int nx, nz, nb, nPoints;
    int srcId, recId, sIdx, sIdz;

    float * S = nullptr;
    float * T = nullptr;

    float * d_T = nullptr;
    float * d_S = nullptr;

    float * seismogram = nullptr;

    int max_spread;
    Geometry * geometry;
    
    std::string parameters;
    std::string data_folder;
    std::string modeling_type;
    std::string modeling_name;

    void set_parameters();
    void initialization();
    void eikonal_solver();
    void set_shot_point();
    void show_information();    
    void compute_seismogram();

    void copy_slowness_to_device();

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);
    
    virtual void time_propagation() = 0;

    void export_seismogram();
};

__global__ void time_set(float * T, int matsize);

__global__ void time_init(float * T, float * S, float sx, float sz, float dx, 
                          float dz, int sIdx, int sIdz, int nzz, int nb);

__global__ void inner_sweep(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, 
                            int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, 
                            float dx, float dz, float dx2i, float dz2i);

# endif
