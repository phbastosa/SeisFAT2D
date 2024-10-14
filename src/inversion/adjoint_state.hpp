# ifndef ADJOINT_STATE_CUH
# define ADJOINT_STATE_CUH

# include "tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    int totalLevels;
    int nSweeps, meshDim;

    int i, j;

    float cell_area;

    float * source = nullptr;       
    float * adjoint = nullptr;
    float * gradient = nullptr;

    void inner_sweep();
    void initialization();
    void set_specifications();
    void gradient_preconditioning();
    void apply_inversion_technique();

public:

    void optimization();

};

# endif