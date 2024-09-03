# include "tomography.hpp"

void Tomography::set_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));
    max_slowness_variation = std::stof(catch_parameter("max_slowness_variation", parameters));

    obs_data_folder = catch_parameter("obs_data_folder", parameters);
    obs_data_prefix = catch_parameter("obs_data_prefix", parameters);

    smooth_model_per_iteration = str2bool(catch_parameter("smooth_per_iteration", parameters));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", parameters));
    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", parameters));
    
    convergence_map_folder = catch_parameter("convergence_folder", parameters);
    estimated_model_folder = catch_parameter("estimated_model_folder", parameters);

    write_model_per_iteration = str2bool(catch_parameter("export_model_per_iteration", parameters));

    set_forward_modeling();
    set_inversion_elements();
    set_specifications();
}

void Tomography::set_forward_modeling()
{
    modeling = new Parallel_aFSM();

    modeling->parameters = parameters;

    modeling->set_parameters();
}

void Tomography::set_inversion_elements()
{    
    ndata = 0;
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
        ndata += modeling->geometry->spread[shot];

    dcal = new float[ndata]();
    dobs = new float[ndata]();

    perturbation = new float[modeling->nPoints]();
}

void Tomography::import_obsData()
{
    int ptr = 0;     
    
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        float * data = new float[modeling->geometry->spread[shot]]();

        std::string path = obs_data_folder + obs_data_prefix + std::to_string(modeling->geometry->sInd[shot]) + ".bin";

        import_binary_float(path, data, modeling->geometry->spread[shot]);
    
        for (int d = ptr; d < ptr + modeling->geometry->spread[shot]; d++) 
            dobs[d] = data[d - ptr];

        ptr += modeling->geometry->spread[shot]; 

        delete[] data;
    }
}

void Tomography::forward_modeling()
{


}

void Tomography::objective_function()
{


}

void Tomography::model_update()
{


}

void Tomography::export_results()
{

    
}