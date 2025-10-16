# ifndef KDM_CUH
# define KDM_CUH

# include "migration.cuh"

class KDM : public Migration
{
    void image_building();
    void export_outputs();

    void prepare_components();
};

__global__ void cross_correlation(float * S, float * Ts, float * Tr, float * data, float * trace, float * angle, float * image, float cmp, float aperture, int nTraces, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz);

# endif