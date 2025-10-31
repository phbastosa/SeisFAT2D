# ifndef IDLSKDM_CUH
# define IDLSKDM_CUH

# include "migration.cuh"

class IDLSKDM : public Migration
{
    float alpha, beta;

    float * dcal = nullptr;
    float * dres = nullptr;
    
    float * model = nullptr;
    
    float * gradient_p = nullptr;
    float * gradient_m = nullptr;
    
    float * direction = nullptr;

    std::vector<float> residuo;

    void forward(float * m, float * d);
    void adjoint(float * d, float * m);

    void prepare_components();
    void check_convergence();
    void set_initial_model();
    void compute_gradient();
    void compute_direction();
    void compute_stepLength();
    void update_reflectivity();
    void show_iteration_info();

    void image_building();
    void export_outputs();
};

__global__ void forward_kernel(float * S, float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz);
__global__ void adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz);

# endif