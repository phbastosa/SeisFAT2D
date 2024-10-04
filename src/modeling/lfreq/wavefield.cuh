# ifndef WAVEFIELD_CUH
# define WAVEFIELD_CUH

# include <complex>
# include <fftw3.h>

# include <cuda_runtime.h>

# include "../modeling.hpp"

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

    float * d1D = nullptr;
    float * d2D = nullptr;

    int * rIdx = nullptr;
    int * rIdz = nullptr;

    int * current_xrec = nullptr;
    int * current_zrec = nullptr;

    float * wavelet = nullptr;

    void define_cerjan_dampers();

    virtual void set_conditions() = 0;
    virtual void set_properties() = 0;

public:

    virtual void initialization() = 0;
    virtual void forward_solver() = 0;
};

# endif