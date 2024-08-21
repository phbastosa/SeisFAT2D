# include "eikonal.hpp"

void Eikonal::set_parameters()
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

    data_folder = catch_parameter("modeling_output_folder", parameters); 

    std::string model_file = catch_parameter("model_file", parameters);  
    
    import_binary_float(model_file, velocity, nPoints);

    expand_boundary(velocity, slowness);

    for (int index = 0; index < matsize; index++) 
        slowness[index] = 1.0f / slowness[index];

    geometry = new Geometry(parameters);

    synthetic_data = new float[geometry->nrec]();

    set_specifications();
}

void Eikonal::expand_boundary(float * input, float * output)
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

void Eikonal::reduce_boundary(float * input, float * output)
{
    for (int index = 0; index < nPoints; index++)
    {
        int x = (int) (index / nz);    
        int z = (int) (index % nz);  

        output[z + x*nz] = input[(z + nb) + (x + nb)*nzz];
    }
}

void Eikonal::show_information()
{
    auto clear = system("clear");

    std::cout << "\033[34mSeis\033[0;0mmic \033[34mF\033[0;0mirst-\033[34mA\033[0;0mrrival \033[34mT\033[0;0moolkit \033[34m2D\033[0;0m\n\n";

    std::cout << "Model dimensions: (z = " << (nz - 1)*dz << ", x = " << (nx - 1) * dx <<") m\n\n";

    std::cout << "Shot " << srcId + 1 << " of " << geometry->nrel;

    std::cout << " at position: (z = " << geometry->zsrc[srcId] << 
                              ", x = " << geometry->xsrc[srcId] << ") m\n";
}

void Eikonal::get_synthetic_data()
{
    int spread = 0;

    for (int recId = geometry->iRec[srcId]; recId < geometry->fRec[srcId]; recId++)
    {
        float x = geometry->xrec[recId];
        float z = geometry->zrec[recId];
 
        float x1 = floorf(x / dx)*dx;  
        float x2 = floorf(x / dx)*dx + dx;

        float z1 = floorf(z / dz)*dz;  
        float z2 = floorf(z / dz)*dz + dz;

        int i1 = (int)(z1 / dz) + nb;
        int i2 = (int)(z2 / dz) + nb;

        int j1 = (int)(x1 / dx) + nb;
        int j2 = (int)(x2 / dx) + nb;

        float q11 = eikonalT[i1 + j1*nzz];
        float q12 = eikonalT[i2 + j1*nzz];
        float q21 = eikonalT[i1 + j2*nzz];
        float q22 = eikonalT[i2 + j2*nzz];

        float p0 = 1.0 / ((x2 - x1) * (z2 - z1));

        float p1 = q11 * (x2 - x) * (z2 - z);
        float p2 = q21 * (x - x1) * (z2 - z);
        float p3 = q12 * (x2 - x) * (z - z1);
        float p4 = q22 * (x - x1) * (z - z1);

        synthetic_data[spread++] = p0*(p1 + p2 + p3 + p4);
    }

    std::string data_file = data_folder + "travel_time_" + std::to_string(spread) + "_stations_shot_" + std::to_string(srcId+1) + ".bin";
    export_binary_float(data_file, synthetic_data, spread);    
}
