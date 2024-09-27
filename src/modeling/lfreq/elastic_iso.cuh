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

    void forward_solver();
};

# endif