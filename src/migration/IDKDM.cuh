# ifndef IDKDM_CUH
# define IDKDM_CUH

# include "KDM.cuh"

class IDKDM : public KDM
{
    void set_migration();
    void initialization();
    void perform_migration();

public:
    
    void export_outputs();
};

__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, int nxx, int nzz, int nt, int nb);

# endif
