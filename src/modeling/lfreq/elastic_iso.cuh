# ifndef ELASTIC_ISO_CUH
# define ELASTIC_ISO_CUH

# include "wavefield.cuh"

class elastic_Iso : public Wavefield
{
private:

    float * M = nullptr;
    float * L = nullptr;
    float * B = nullptr;
    float * P = nullptr;

    float * d_M = nullptr;
    float * d_L = nullptr;
    float * d_B = nullptr;
    float * d_P = nullptr;

    float * d_Vx = nullptr;
    float * d_Vz = nullptr;
    float * d_Txx = nullptr;
    float * d_Tzz = nullptr;
    float * d_Txz = nullptr;

    void set_properties();
    void set_conditions();

public:

    void initialization();
    void forward_solver();
};

__global__ void compute_pressure(float * P, float * Vx, float * Vz, float * Txx, float * Tzz, float * Txz, float * M, float * L, float * wavelet, int sIdx, int sIdz, int tId, int nt, int nxx, int nzz, float dx, float dz, float dt);
__global__ void compute_velocity(float * Vx, float * Vz, float * Txx, float * Tzz, float * Txz, float * B, int nxx, int nzz, float dx, float dz, float dt);

# endif