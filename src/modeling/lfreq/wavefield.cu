# include "wavefield.cuh"

void Wavefield::set_specifications()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));

    fmax = std::stof(catch_parameter("max_frequency", parameters));

    set_wavelet();
    set_boundaries();
    set_properties();    
    set_conditions();    
}

void Wavefield::set_wavelet()
{



    
}

void Wavefield::set_boundaries()
{
    nb = std::stoi(catch_parameter("boundary_samples", parameters));

    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    matsize = nxx*nzz;
}

void Wavefield::initialization()
{


}