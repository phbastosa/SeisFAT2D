# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../ioFunctions/ioFunctions.hpp"

class Tomography
{
private:

protected:

    int max_iteration;

    bool write_model_per_iteration;
    bool smooth_model_per_iteration;
    int smoother_samples;
    float smoother_stdv;

    float max_slowness_variation;

    std::string obs_data_folder;
    std::string obs_data_prefix;
    std::string convergence_map_folder;
    std::string estimated_model_folder;


public:

    std::string parameters;

    void set_parameters();

};


# endif