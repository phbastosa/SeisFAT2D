# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));

    aperture = std::stof(catch_parameter("mig_aperture", parameters));
    max_offset = std::stof(catch_parameter("max_offset", parameters));    

    input_data_folder = catch_parameter("input_data_folder", parameters);
    input_data_prefix = catch_parameter("input_data_prefix", parameters);

    output_image_folder = catch_parameter("output_image_folder", parameters);
    output_table_folder = catch_parameter("output_table_folder", parameters);

    set_modeling_type();

    modeling->parameters = parameters;
    modeling->set_parameters();

    f_image = new float[modeling->nPoints]();
    h_image = new float[modeling->matsize]();
    h_seismic = new float[nt*modeling->max_spread]();

    cudaMalloc((void**)&(d_Tr), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_image), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_seismic), nt*modeling->max_spread*sizeof(float));

    nThreads = 256;
    nBlocks = (int)((modeling->matsize + nThreads - 1) / nThreads);
}

void Migration::read_seismic_data()
{
    std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

    import_binary_float(data_path, h_seismic, nt*modeling->geometry->spread[modeling->srcId]);
    
    cudaMemcpy(d_seismic, h_seismic, nt*modeling->geometry->spread[modeling->srcId]*sizeof(float), cudaMemcpyHostToDevice);
}

void Migration::image_building()
{
    get_receiver_eikonal();
    run_cross_correlation();
}

void Migration::get_receiver_eikonal()
{
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        set_receiver_point();

        show_information();

        modeling->time_propagation();
        
        export_receiver_eikonal();
    }
}

void Migration::set_receiver_point()
{
    modeling->sx = modeling->geometry->xrec[modeling->recId];
    modeling->sz = modeling->geometry->zrec[modeling->recId];
}

void Migration::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------------------\n";
    std::cout << "                                 \033[34mSeisFAT2D\033[0;0m\n";
    std::cout << "-------------------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (modeling->nz - 1)*modeling->dz << 
                                  ", x = " << (modeling->nx - 1)*modeling->dx << ") m\n\n";

    std::cout << "Running receiver " << modeling->recId + 1 << " of " << modeling->geometry->nrec << " in total\n\n";

    std::cout << "Current receiver position: (z = " << modeling->geometry->zrec[modeling->recId] << 
                                           ", x = " << modeling->geometry->xrec[modeling->recId] << ") m\n\n";

    std::cout << "Kirchhoff Depth Migration: computing receiver travel time matrices\n";
}

void Migration::export_receiver_eikonal()
{
    cudaMemcpy(modeling->T, modeling->d_T, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float(output_table_folder + "eikonal_receiver_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->matsize);    
}

void Migration::run_cross_correlation()
{
    cudaMemset(d_image, 0.0f, modeling->nPoints*sizeof(float));

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        modeling->set_shot_point();
        modeling->show_information();
        modeling->time_propagation();

        std::cout << "\nKirchhoff depth migration: computing image matrix\n";

        read_seismic_data();

        int spread = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            float rx = modeling->geometry->xrec[modeling->recId];

            float offset = sqrtf(powf(modeling->sx - rx, 2.0f));

            if (offset < max_offset)
            {
                float cmp = modeling->sx + 0.5f*(rx - modeling->sx);

                import_binary_float(output_table_folder + "eikonal_receiver_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->matsize);
            
                cudaMemcpy(d_Tr, modeling->T, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

                cross_correlation<<<nBlocks,nThreads>>>(modeling->d_T, d_Tr, d_image, d_seismic, aperture, cmp, spread, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);
            }

            ++spread;
        }
    }
}

void Migration::export_outputs()
{
    cudaMemcpy(h_image, d_image, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

    # pragma omp parallel for
    for (int index = 0; index < modeling->matsize; index++)
    {
        int i = (int)(index % modeling->nzz);
        int j = (int)(index / modeling->nzz);

        if ((i > modeling->nb) && (i < modeling->nzz-modeling->nb) && 
            (j > modeling->nb) && (j < modeling->nxx-modeling->nb))
        {
            int xp = i + (j+1)*modeling->nzz;
            int xm = i + (j-1)*modeling->nzz;

            int zp = (i+1) + j*modeling->nzz;
            int zm = (i-1) + j*modeling->nzz;

            float dIm_dx = (h_image[xp] - 2.0f*h_image[index] + h_image[xm])/modeling->dx/modeling->dx;
            float dIm_dz = (h_image[zp] - 2.0f*h_image[index] + h_image[zm])/modeling->dz/modeling->dz;
        
            int inner = (i-modeling->nb) + (j-modeling->nb)*modeling->nz;

            f_image[inner] = (dIm_dx + dIm_dz) / modeling->geometry->nrel;
        }
    }

    export_binary_float(output_image_folder + "kirchhoff_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin", f_image, modeling->nPoints);
}

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture, float cmp, int spread, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float sigma = tanf(aperture * M_PI / 180.0f)*(i-nb)*dz;

        float value = expf(-0.5f*powf(((j-nb)*dx - cmp)/(sigma + 1e-6f), 2.0f));

        float T = Ts[index] + Tr[index]; 
    
        int tId = (int)(T / dt);

        if (tId < nt) image[index] += value * seismic[tId + spread*nt];
    }
}
