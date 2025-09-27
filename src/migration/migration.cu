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

    set_common_gathers();

    h_ODCIG = new float[modeling->nz]();
    h_ADCIG = new float[modeling->nz]();

    h_image = new float[modeling->matsize]();
    f_image = new float[modeling->nPoints]();
    h_seismic = new float[nt*modeling->max_spread]();

    cudaMalloc((void**)&(d_ODCIG), modeling->nz*sizeof(float));
    cudaMalloc((void**)&(d_ADCIG), modeling->nz*sizeof(float));

    cudaMalloc((void**)&(d_Tr), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_image), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_seismic), nt*modeling->max_spread*sizeof(float));

    nThreads = 256;
    nBlocks = (int)((modeling->matsize + nThreads - 1) / nThreads);
}

void Migration::set_common_gathers()
{
    da = 5.0f;
    nang = 12; 

    ds = modeling->geometry->xsrc[1] - modeling->geometry->xsrc[0];
    dr = modeling->geometry->xrec[1] - modeling->geometry->xrec[0];

    std::map<int,int> counts;    

    for (int srcId = 0; srcId < modeling->geometry->nrel; srcId++)
    {
        int spreadId = 0;
        for (int recId = modeling->geometry->iRec[srcId]; recId < modeling->geometry->fRec[srcId]; recId++)
        {
            int midpoint = spreadId + 2.0*(ds/dr)*srcId; 

            counts[midpoint]++;

            ++spreadId;
        }
    }

    nTraces = 0;
    ncmp = counts.size();
    partial_cmp_sum = new int[ncmp]();

    int * traces_per_cmp = new int[ncmp]();

    int i = 0;
    for (auto& count : counts) 
    {
        nTraces += count.second;
        traces_per_cmp[i++] = count.second;
    }

    for (int i = 0; i < ncmp; i++)
        for (int j = 0; j < i; j++)    
            partial_cmp_sum[i] += traces_per_cmp[j];

    ODCIG = new float[modeling->nz * nTraces]();
    ADCIG = new float[modeling->nz * ncmp*nang]();

    delete[] traces_per_cmp;
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

    std::vector<int> localCount(ncmp, 0);    

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        modeling->set_shot_point();
        modeling->show_information();
        modeling->time_propagation();

        std::cout << "\nKirchhoff depth migration: computing image matrix\n";

        read_seismic_data();

        int spreadId = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            float rx = modeling->geometry->xrec[modeling->recId];

            float offset = modeling->sx - rx;

            if (fabsf(offset) < max_offset)
            {
                float cmp = 0.5f*(modeling->sx + rx);

                import_binary_float(output_table_folder + "eikonal_receiver_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->matsize);
            
                cudaMemcpy(d_Tr, modeling->T, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

                cross_correlation<<<nBlocks,nThreads>>>(modeling->d_S, modeling->d_T, d_Tr, d_image, d_seismic, d_ODCIG, d_ADCIG, aperture, cmp, spreadId, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);

                cudaMemcpy(h_ODCIG, d_ODCIG, modeling->nz*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(h_ADCIG, d_ADCIG, modeling->nz*sizeof(float), cudaMemcpyDeviceToHost);
                
                int cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId; 
                int traceId = localCount[cmpId]++;                         
                int gatherId = partial_cmp_sum[cmpId] + traceId;    

                for (int i = 0; i < modeling->nz; i++)
                {
                    ODCIG[i + gatherId*modeling->nz] = h_ODCIG[i];

                    float angle = 180.0f*h_ADCIG[i]/M_PI/da;
                    int angleId = static_cast<int>(angle);

                    if ((angle >= 0) && (angleId < nang))  
                        ADCIG[i + angleId*modeling->nz + cmpId*nang*modeling->nz] = h_ODCIG[i];
                }
            }

            ++spreadId;
        }
    }
}

void Migration::export_outputs()
{
    cudaMemcpy(h_image, d_image, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);
    modeling->reduce_boundary(h_image, f_image);

    export_binary_float(output_image_folder + "kirchhoff_section_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin", f_image, modeling->nPoints);
    export_binary_float(output_image_folder + "kirchhoff_ODCIG_" + std::to_string(modeling->nz) + "x" + std::to_string(nTraces) + ".bin", ODCIG, modeling->nz * nTraces);
    export_binary_float(output_image_folder + "kirchhoff_ADCIG_" + std::to_string(modeling->nz) + "x" + std::to_string(ncmp*nang) + ".bin", ADCIG, modeling->nz * ncmp*nang);
}

__global__ void cross_correlation(float * S, float * Ts, float * Tr, float * image, float * seismic, float * ODCIG, float * ADCIG, float aperture, float cmp, int spread, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
    
        int tId = (int)(T / dt);

        if (tId < nt) 
        {   
            float eps = 1e-6f;

            float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
            float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

            float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
            float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
            float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

            float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
            float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + eps;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + eps;

            float dot_product = dTs_dx*dTr_dx + dTs_dz*dTr_dz;

            float angle = 0.5f*acos(dot_product / (norm_Ts * norm_Tr));

            float nvr_x = sinf(angle);
            float nvr_z = cosf(angle);

            float norm_nvr = sqrtf(nvr_x*nvr_x + nvr_z*nvr_z) + eps;

            nvr_x = nvr_x / norm_nvr;
            nvr_z = nvr_z / norm_nvr;

            float tvr_x = nvr_z;
            float tvr_z =-nvr_x;

            float ps_x = dTs_dx / norm_Ts;
            float ps_z = dTs_dz / norm_Ts;

            float pr_x = dTr_dx / norm_Tr;
            float pr_z = dTr_dz / norm_Tr;

            float cos_s = fabsf(ps_x*nvr_x + ps_z*nvr_z);
            float cos_r = fabsf(pr_x*nvr_x + pr_z*nvr_z);

            float N1 = d2Ts_dx2*tvr_x*tvr_x + 2.0f*d2Ts_dxdz*tvr_x*tvr_z + d2Ts_dz2*tvr_z*tvr_z;
            float N2 = d2Tr_dx2*tvr_x*tvr_x + 2.0f*d2Tr_dxdz*tvr_x*tvr_z + d2Tr_dz2*tvr_z*tvr_z;

            float weights = S[index]*sqrtf(cos_s*cos_r)*fabsf(N1 + N2) / (sqrtf(fabsf(N1 * N2)) + eps);

            float sigma = tanf(aperture * M_PI / 180.0f)*(i-nb)*dz;

            float value = expf(-0.5f*powf(((j-nb)*dx - cmp)/(sigma + eps), 2.0f));

            float amplitude = value * weights * seismic[tId + spread*nt];
            
            image[index] += amplitude;
            
            if (j == (int)(cmp / dx) + nb)
            {
                ADCIG[i - nb] = angle;    
                ODCIG[i - nb] = amplitude;
            }
        }    
    }
}
