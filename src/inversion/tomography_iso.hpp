# ifndef LEAST_SQUARES_HPP
# define LEAST_SQUARES_HPP

# include "tomography.hpp"

class Least_Squares : public Tomography
{
private:

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
    void apply_inversion_technique();

public:

    void optimization();
};

# endif
