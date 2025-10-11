# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));
    da = std::stof(catch_parameter("mig_dangle", parameters));

    fmax = std::stof(catch_parameter("max_frequency", parameters));
    max_it = std::stoi(catch_parameter("mig_max_it", parameters));
    aperture = std::stof(catch_parameter("mig_aperture", parameters));
    max_angle = std::stof(catch_parameter("mig_max_angle", parameters));

    input_data_folder = catch_parameter("mig_data_folder", parameters);
    input_data_prefix = catch_parameter("mig_data_prefix", parameters);
    
    tables_folder = catch_parameter("mig_tables_folder", parameters);
    outputs_folder = catch_parameter("mig_outputs_folder", parameters);

    anisotropy = str2bool(catch_parameter("anisotropy", parameters));

    if (anisotropy) modeling = new Eikonal_ANI(); 
               else modeling = new Eikonal_ISO();

    modeling->parameters = parameters;
    modeling->set_parameters();

    nang = (int)(max_angle / da) + 1;

    nBlocks = (int)((modeling->matsize + NTHREADS - 1) / NTHREADS);

    set_wavelet();
    set_gathers();

    h_trace = new float[modeling->nz]();
    h_angle = new float[modeling->nz]();
    h_image = new float[modeling->matsize](); 
    seismic = new float[nt*modeling->max_spread]();

    cudaMalloc((void**)&(d_data), nt*sizeof(float));    
    
    cudaMalloc((void**)&(d_Tr), modeling->matsize*sizeof(float));

    cudaMalloc((void**)&(d_trace), modeling->nz*sizeof(float));
    cudaMalloc((void**)&(d_angle), modeling->nz*sizeof(float));
    cudaMalloc((void**)&(d_image), modeling->matsize*sizeof(float));
}

void Migration::set_wavelet()
{
    float t0 = 2.0f*sqrtf(M_PI) / fmax;
    float fc = fmax / (3.0f * sqrtf(M_PI));

    nw = 2*((int)((t0 - 0.5f*dt) / dt));

    wavelet = new float[nw]();

    for (int wId = 0; wId < nw; wId++)
    {
        float td = wId*dt - t0;

        float arg = M_PI*M_PI*M_PI*fc*fc*td*td;

        wavelet[wId] = (1.0f - 2.0f*arg)*expf(-arg);
    }
}

void Migration::set_gathers()
{
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
    nCMP = counts.size();
    pSUM = new int[nCMP]();

    int * traces_per_cmp = new int[nCMP]();

    int i = 0;
    for (auto& count : counts) 
    {
        nTraces += count.second;
        traces_per_cmp[i++] = count.second;
    }

    for (int i = 0; i < nCMP; i++)
        for (int j = 0; j < i; j++)    
            pSUM[i] += traces_per_cmp[j];

    IMAGE = new float[modeling->nPoints]();
    ODCIG = new float[modeling->nz * nTraces]();
    ADCIG = new float[modeling->nz * nCMP*nang]();

    delete[] traces_per_cmp;
}

void Migration::set_rec_travel_times()
{
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        modeling->sx = modeling->geometry->xrec[modeling->recId];
        modeling->sz = modeling->geometry->zrec[modeling->recId];

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(modeling->T, modeling->d_T, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->matsize);
    }
}

void Migration::prepare_convolution()
{
    nfft = nextpow2(nt);

    trace_in = new float[nt]();
    trace_out = new float[nt]();

    time_trace = (double *) fftw_malloc(nfft*sizeof(double));
    time_result = (double *) fftw_malloc(nfft*sizeof(double));
    time_wavelet = (double *) fftw_malloc(nfft*sizeof(double));

    freq_trace = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));
    freq_result = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));
    freq_wavelet = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));

    trace_forward_plan = fftw_plan_dft_r2c_1d(nt, time_trace, freq_trace, FFTW_ESTIMATE);
    result_adjoint_plan = fftw_plan_dft_c2r_1d(nt, freq_result, time_result, FFTW_ESTIMATE);
    wavelet_forward_plan = fftw_plan_dft_r2c_1d(nt, time_wavelet, freq_wavelet, FFTW_ESTIMATE);

    for (int tId = 0; tId < nfft; tId++)
    {
        time_trace[tId] = 0.0;
        time_result[tId] = 0.0;
        time_wavelet[tId] = (tId < nw) ? wavelet[tId] : 0.0;
    }

    fftw_execute(wavelet_forward_plan);
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
    
    std::cout << "Computing receiver travel time matrices\n";
}
