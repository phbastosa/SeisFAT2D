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
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        float * data = new float[modeling->geometry->spread[shot]]();

        std::string path = obs_data_folder + obs_data_prefix + std::to_string(modeling->geometry->sInd[shot]) + ".bin";

        import_binary_float(path, data, modeling->geometry->spread[shot]);

        int skipped = shot * modeling->geometry->spread[shot];    
        for (int i = 0; i < modeling->geometry->spread[shot]; i++) 
            dobs[i + skipped] = data[i];

        delete[] data;
    }
}

void Tomography::forward_modeling()
{
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        modeling->srcId = shot;
    
        modeling->show_information();

        show_information();

        modeling->initialization();
        modeling->forward_solver();

        modeling->get_synthetic_data();

        concatenate_synthetic_data();
        
        if (iteration != max_iteration)
            apply_inversion_technique();
    }

    // if (iteration != max_iteration)
    // {
    //     gradient_preconditioning();
    //     export_gradient();
    // }    
}

void Tomography::show_information()
{
    std::cout << "\nInversion type: ";
    std::cout << inversion_method << "\n\n";

    if (iteration == max_iteration)
    { 
        std::cout << "Checking final residuo \n\n";
    }
    else
        std::cout << "Computing iteration " << iteration + 1 << " out of " << max_iteration << "\n\n";

    if (iteration > 0) std::cout << "Previous residuo: " << residuo.back() << "\n\n";   
}

void Tomography::concatenate_synthetic_data()
{
    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    for (int i = 0; i < modeling->geometry->spread[modeling->srcId]; i++) 
        dcal[i + skipped] = modeling->synthetic_data[i];    
}

void Tomography::check_convergence()
{
    float square_difference = 0.0f;

    for (int i = 0; i < ndata; i++)
        square_difference += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(square_difference));

    if ((iteration >= max_iteration))
    {
        std::cout << "\nFinal residuo: "<< residuo.back() <<"\n\n";
        converged = true;
    }
    else
    {
        iteration += 1;
        converged = false;
    }
}

void Tomography::model_update()
{


}

void Tomography::export_results()
{

    
}