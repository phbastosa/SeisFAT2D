# include "modeling.hpp"

void Modeling::set_parameters()
{
    nb = 1;

    nx = std::stoi(catch_parameter("x_samples", parameters));        
    nz = std::stoi(catch_parameter("z_samples", parameters));        

    dx = std::stof(catch_parameter("x_spacing", parameters));
    dz = std::stof(catch_parameter("z_spacing", parameters));

    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    nPoints = nx*nz;
    matsize = nxx*nzz;

    velocity = new float[nPoints]();
    slowness = new float[matsize]();
    eikonalT = new float[matsize]();

    std::string model_file = catch_parameter("model_file", parameters);  

    import_binary_float(model_file, velocity, nPoints);

    expand_boundary(velocity, slowness);

    for (int index = 0; index < matsize; index++) 
        slowness[index] = 1.0f / slowness[index];

    geometry = new Geometry(parameters);

    check_geometry_bounds();
    
    set_eikonal_parameters();
}

void Modeling::check_geometry_bounds()
{
    for (int i = 0; i < geometry->nsrc; i++)
    {
        if ((geometry->xsrc[i] < 0) || (geometry->xsrc[i] > (nx-1)*dx) || 
            (geometry->zsrc[i] < 0) || (geometry->zsrc[i] > (nz-1)*dz))
        throw std::invalid_argument("\033[31mError: shots geometry overflow!\033[0;0m");
    }

    for (int i = 0; i < geometry->nrec; i++)
    {
        if ((geometry->xrec[i] < 0) || (geometry->xrec[i] > (nx-1)*dx) || 
            (geometry->zrec[i] < 0) || (geometry->zrec[i] > (nz-1)*dz))
        throw std::invalid_argument("\033[31mError: nodes geometry overflow!\033[0;0m");
    }
}

void Modeling::expand_boundary(float * input, float * output)
{
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            output[(i + nb) + (j + nb)*nzz] = input[i + j*nz];
        }
    }

    for (int i = 0; i < nb; i++)
    {
        for (int j = nb; j < nxx - nb; j++)
        {
            output[i + j*nzz] = output[nb + j*nzz];
            output[(nzz - i - 1) + j*nzz] = output[(nzz - nb - 1) + j*nzz];
        }
    }

    for (int i = 0; i < nzz; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            output[i + j*nzz] = output[i + nb*nzz];
            output[i + (nxx - j - 1)*nzz] = output[i + (nxx - nb - 1)*nzz];
        }
    }
}

void Modeling::reduce_boundary(float * input, float * output)
{
    for (int index = 0; index < nPoints; index++)
    {
        int x = (int) (index / nz);    
        int z = (int) (index % nz);  

        output[z + x*nz] = input[(z + nb) + (x + nb)*nzz];
    }
}

// void Modeling::fast_sweeping_method_CPU()
// {
//     float dz2i = 1.0f / (dz*dz);    
//     float dx2i = 1.0f / (dx*dx);    

//     float sgntz, sgntx, sgnvx, sgnvz;

//     sgntz = 1; sgnvz = 1;
//     for (int i = 1; i < nzz; i++)
//     {
//         sgntx = 1; sgnvx = 1;
//         for (int j = 1; j < nxx; j++)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        

//         sgntx = -1; sgnvx = 0;
//         for (int j = nxx - 2; j > -1; j--)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        
//     }

//     sgntz = -1; sgnvz = 0;
//     for (int i = nzz-2; i > -1; i--)
//     {
//         sgntx = 1; sgnvx = 1;
//         for (int j = 1; j < nxx; j++)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        

//         sgntx = -1; sgnvx = 0;
//         for (int j = nxx-2; j > -1; j--)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        
//     }        

//     sgntx = 1; sgnvx = 1;
//     for (int j = 1; j < nxx; j++) 
//     {
//         sgntz = 1; sgnvz = 1;
//         for (int i = 1; i < nzz; i++)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        

//         sgntz = -1; sgnvz = 0;
//         for (int i = nzz-2; i > -1; i--)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        
//     }        

//     sgntx = -1; sgnvx = 0;
//     for (int j = nxx-2; j > -1; j--)
//     {
//         sgntz = 1; sgnvz = 1;
//         for (int i = 1; i < nzz; i++)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        

//         sgntz = -1; sgnvz = 0;
//         for (int i = nzz-2; i > -1; i--)
//             inner_sweep(i,j,sgnvx,sgnvz,sgntx,sgntz,dx2i,dz2i);        
//     }
// }

