# ifndef LEAST_SQUARES_HPP
# define LEAST_SQUARES_HPP

# include "tomography.hpp"

class Least_Squares : public Tomography
{
private:

    int nx_tomo, nz_tomo, n_model;
    float dx_tomo, dz_tomo;

    int tk_order;
    float tk_param;

    size_t ray_path_max_samples;

    std::vector<int> iG;
    std::vector<int> jG;
    std::vector<float> vG;

    int M, N, NNZ;

    int * iA = nullptr;
    int * jA = nullptr;
    float * vA = nullptr;
    float * B = nullptr;
    float * x = nullptr; 

    void set_specifications();

    void apply_regularization();
    void solve_linear_system_lscg();
    void slowness_variation_rescaling();

    void apply_inversion_technique();
    void set_specific_parameters();
    void gradient_preconditioning();

public:

    void optimization();

};

# endif