# ifndef TOMOGRAPHY_ANI_HPP
# define TOMOGRAPHY_ANI_HPP

# include "inversion.hpp"

class Tomography_ANI : public Inversion
{
    float * E = nullptr;
    float * D = nullptr;

    void set_modeling_type();
    void set_sensitivity_matrix();
    void set_regularization();
};

# endif
