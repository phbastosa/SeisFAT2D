# ifndef ADJOINT_STATE_HPP
# define ADJOINT_STATE_HPP

# include <complex>
# include <fftw3.h>

# include "tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    int i, j;

    float cell_area;

    float * source_grad = nullptr;
    float * source_comp = nullptr;

    float * adjoint_grad = nullptr;
    float * adjoint_comp = nullptr;
    
    float * gradient = nullptr;
    float * illumination = nullptr;

    void inner_sweep();
    void initialization();
    void set_specifications();
    void gradient_preconditioning();
    void apply_inversion_technique();

public:

    void optimization();

};

# endif