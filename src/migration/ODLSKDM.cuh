# ifndef ODLSKDM_CUH
# define ODLSKDM_CUH

# include "migration.cuh"

class ODLSKDM : public Migration
{
    void forward(float * m, float * d);
    void adjoint(float * d, float * m);

    void image_building();
    void export_outputs();
};

# endif