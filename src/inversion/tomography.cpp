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


}