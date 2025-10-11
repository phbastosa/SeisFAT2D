# ifndef ADLSKDM_CUH
# define ADLSKDM_CUH

# include "migration.cuh"

class ADLSKDM : public Migration
{
    void forward();
    void adjoint();

    void image_building();
    void export_outputs();
};

# endif