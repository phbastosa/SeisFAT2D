# ifndef WAVEFIELD_CUH
# define WAVEFIELD_CUH

# include "../modeling.hpp"

# include <cuda_runtime.h>

class Wavefield : public Modeling
{
private:

    void set_wavelet();
    void set_boundaries();
    void set_specifications();

protected:

    float fmax;

    int tlag, nThreads;
    int sBlocks, nBlocks;

    int * rIdx = nullptr;
    int * rIdz = nullptr;

    int * current_xrec = nullptr;
    int * current_zrec = nullptr;

    float * wavelet = nullptr;

    virtual void set_conditions() = 0;
    virtual void set_properties() = 0;

public:

    virtual void initialization() = 0;
    virtual void forward_solver() = 0;
};

# endif