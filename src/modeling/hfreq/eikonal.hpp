# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

protected:

public:

    virtual void set_name() = 0;

    // virtual void set_specifications() = 0;

    // float dx, dz;
    // int srcId, recId;
    // int nxx, nzz, matsize;
    // int nx, nz, nb, nPoints;

    // float * velocity = nullptr;
    // float * slowness = nullptr;
    // float * eikonalT = nullptr;

    // float * synthetic_data = nullptr;

    // std::string data_folder;

    // void expand_boundary(float * input, float * output);
    // void reduce_boundary(float * input, float * output);

    // Geometry * geometry;

    // std::string parameters;
    
    // void set_parameters();
    // void show_information();    
    // void get_synthetic_data();
    // void export_synthetic_data();

    // virtual void initialization() = 0;
    // virtual void forward_solver() = 0;
};

# endif