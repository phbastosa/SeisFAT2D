# ifndef ADLSKDM_CUH
# define ADLSKDM_CUH

# include "migration.cuh"

class ADLSKDM : public Migration
{
    void forward(float * m, float * d);
    void adjoint(float * d, float * m);

    void image_building();
    void export_outputs();
};

# endif