# ifndef ADLSKDM_CUH
# define ADLSKDM_CUH

# include "LSKDM.cuh"

class ADLSKDM : public LSKDM
{
    void set_migration();
    void perform_forward();
    void perform_adjoint();
};

# endif