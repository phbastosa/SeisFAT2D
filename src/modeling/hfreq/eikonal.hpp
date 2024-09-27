# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    void set_boundaries();
    void set_specifications();

protected:

    virtual void set_conditions() = 0;
    virtual void set_properties() = 0;

public:

    // std::string data_folder;

    void initialization();

    // void get_synthetic_data();
    // void export_synthetic_data();

    virtual void forward_solver() = 0;
};

# endif