# ifndef ODLSKDM_CUH
# define ODLSKDM_CUH

# include "migration.cuh"

class ODLSKDM : public Migration
{
    void forward();
    void adjoint();

    void image_building();
    void export_outputs();
};

# endif