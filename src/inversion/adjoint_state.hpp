# ifndef ADJOINT_STATE_CUH
# define ADJOINT_STATE_CUH

# include "tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    void set_specifications();
    void apply_inversion_technique();

public:

    void optimization();

};

# endif