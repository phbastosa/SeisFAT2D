# ifndef ADKDM_CUH
# define ADKDM_CUH

# include "KDM.cuh"

class ADKDM : public KDM
{
    void set_migration();
    void initialization();
    void perform_migration();

public:
    
    void export_outputs();
};

__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);

# endif
