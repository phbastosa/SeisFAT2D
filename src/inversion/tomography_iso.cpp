# include "tomography_iso.hpp"

void Tomography_ISO::set_forward_modeling()
{
    modeling = new Eikonal_ISO();
    modeling->parameters = parameters;
    modeling->set_parameters();
}