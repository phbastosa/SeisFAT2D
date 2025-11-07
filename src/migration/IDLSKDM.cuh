# ifndef IDLSKDM_CUH
# define IDLSKDM_CUH

# include "LSKDM.cuh"

class IDLSKDM : public LSKDM
{
    void set_migration();
    void perform_forward();
    void perform_adjoint();
};


# endif