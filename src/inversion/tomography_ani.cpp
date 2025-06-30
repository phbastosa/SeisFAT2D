# include "tomography_ani.hpp"

void Tomography_ANI::set_forward_modeling()
{
    modeling = new Eikonal_ANI();
    modeling->parameters = parameters;
    modeling->set_parameters();
}