# ifndef KIRCHHOFF_HPP
# define KIRCHHOFF_HPP

# include "../modeling/modeling.hpp"

# include "../modeling/hfreq/eikonal_iso.hpp"

class Kirchhoff
{
private:

    std::string input_data_folder;
    std::string input_data_prefix;

protected:

    int nt;
    
    float dt;

    float * Tr = nullptr;
    float * Ts = nullptr;
    float * Im = nullptr;

    float * image = nullptr;
    float * seismic = nullptr;

    Modeling * modeling = nullptr;

    void initialization();
    void run_cross_correlation();
    void get_receiver_traveltimes();
    void export_receiver_traveltimes();

public:
    
    std::string parameters;

    void set_parameters();
    
    void read_seismic_data();

    void image_building();
    void export_outputs();
};

# endif