# ifndef TOMOGRAPHY_HPP
# define TOMOGRAPHY_HPP

# include "../ioFunctions/ioFunctions.hpp"
# include "../modeling/parallel_aFSM.cuh"

class Tomography
{
protected:

    int ndata, max_iteration;

    bool write_model_per_iteration;
    bool smooth_model_per_iteration;
    int smoother_samples;
    float smoother_stdv;

    float max_slowness_variation;

    float * dcal = nullptr;
    float * dobs = nullptr;

    float * perturbation = nullptr;

    std::string obs_data_folder;
    std::string obs_data_prefix;
    std::string convergence_map_folder;
    std::string estimated_model_folder;

    Eikonal * modeling;

    void set_forward_modeling();
    void set_inversion_elements();
    
    virtual void set_specifications() = 0;

public:

    std::string parameters;

    void set_parameters();
    void import_obsData();

    void forward_modeling();
    void objective_function();

    virtual void optimization() = 0;
    
    void model_update();

    void export_results();
};


# endif