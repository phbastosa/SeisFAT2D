# ifndef EIKONAL_ISO_HPP
# define EIKONAL_ISO_HPP

# include "eikonal.hpp"

class Eikonal_Iso : public Eikonal
{
private:

    int i, i1, sgntx, sgntz;
    int j, j1, sgnvx, sgnvz;

    float Sref, t1D, t2D, t3D;
    float t1d1, t1d2, t1, t2, t3; 
    float dz2i, dx2i, tv, te, tev; 

    void inner_sweep();

    void set_properties();
    void set_conditions();
    
public:

    void forward_solver();
};

# endif