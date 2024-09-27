# ifndef WAVEFIELD_CUH
# define WAVEFIELD_CUH

# include "../modeling.hpp"

# include <cuda_runtime.h>

class Wavefield : public Modeling
{
private:

    float fmax, tlag;

    float * h_wavelet = nullptr;
    float * d_wavelet = nullptr;

    void set_wavelet();
    void set_boundaries();
    void set_specifications();

protected:

    virtual void set_conditions() = 0;
    virtual void set_properties() = 0;

public:

    void initialization();

    virtual void forward_solver() = 0;

};

# endif