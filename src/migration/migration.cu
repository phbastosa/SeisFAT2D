# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));
    da = std::stof(catch_parameter("mig_dangle", parameters));

    fmax = std::stof(catch_parameter("max_frequency", parameters));
    aperture = std::stof(catch_parameter("mig_aperture", parameters));
    max_angle = std::stof(catch_parameter("mig_max_angle", parameters));

    input_data_folder = catch_parameter("mig_data_folder", parameters);
    input_data_prefix = catch_parameter("mig_data_prefix", parameters);
    
    tables_folder = catch_parameter("mig_tables_folder", parameters);
    images_folder = catch_parameter("mig_images_folder", parameters);
    gathers_folder = catch_parameter("mig_gathers_folder", parameters);

    anisotropy = str2bool(catch_parameter("anisotropy", parameters));

    if (anisotropy) modeling = new Eikonal_ANI(); 
               else modeling = new Eikonal_ISO();

    modeling->parameters = parameters;
    modeling->set_parameters();

    nang = (int)(max_angle / da) + 1;

    nBlocks = (int)((modeling->matsize + NTHREADS - 1) / NTHREADS);

    set_wavelet();
    set_gathers();
    
    set_migration();

    h_data = new float[nt]();    

    h_Ts = new float[modeling->matsize]();
    h_Tr = new float[modeling->matsize]();
    
    seismic = new float[nt*modeling->max_spread]();

    cudaMalloc((void**)&(d_data), nt*sizeof(float));    
    
    cudaMalloc((void**)&(d_Ts), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_Tr), modeling->matsize*sizeof(float));
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

    nCMP = counts.size();
}

void Migration::set_src_domain()
{
    keyword = "source";
    total = std::to_string(modeling->geometry->nsrc); 
}

void Migration::set_current_src()
{
    modeling->sx = modeling->geometry->xsrc[modeling->srcId];
    modeling->sz = modeling->geometry->zsrc[modeling->srcId];

    current = std::to_string(modeling->srcId+1);
    
    xpos = format1Decimal(modeling->sx);
    zpos = format1Decimal(modeling->sz);

    current_operation = "Computing " + keyword + " travel time matrices";
}

void Migration::set_rec_domain()
{
    keyword = "receiver";
    total = std::to_string(modeling->geometry->nrec); 
}

void Migration::set_current_rec()
{
    modeling->sx = modeling->geometry->xrec[modeling->recId];
    modeling->sz = modeling->geometry->zrec[modeling->recId];

    current = std::to_string(modeling->recId+1);
    
    xpos = format1Decimal(modeling->sx);
    zpos = format1Decimal(modeling->sz);

    current_operation = "Computing " + keyword + " travel time matrices";
}

void Migration::set_src_travel_times()
{
    set_src_domain();

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        set_current_src();

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(h_Ts, modeling->d_T, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    }
}

void Migration::set_rec_travel_times()
{
    set_rec_domain();
    
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        set_current_rec();

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(h_Tr, modeling->d_T, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    }
}

void Migration::prepare_convolution()
{
    nfft = nextpow2(nt);

    time_trace = (double *) fftw_malloc(nfft*sizeof(double));
    time_wavelet = (double *) fftw_malloc(nfft*sizeof(double));

    freq_trace = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));
    freq_wavelet = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));

    trace_forward_plan = fftw_plan_dft_r2c_1d(nt, time_trace, freq_trace, FFTW_ESTIMATE);
    trace_inverse_plan = fftw_plan_dft_c2r_1d(nt, freq_trace, time_trace, FFTW_ESTIMATE);
    wavelet_forward_plan = fftw_plan_dft_r2c_1d(nt, time_wavelet, freq_wavelet, FFTW_ESTIMATE);

    for (int tId = 0; tId < nfft; tId++)
    {
        time_trace[tId] = 0.0;
        time_wavelet[tId] = (tId < nw) ? wavelet[tId] : 0.0;
    }

    fftw_execute(wavelet_forward_plan);
}

void Migration::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------\n";
    std::cout << " \033[34mSeisFAT2D\033[0;0m --------------------------------------------------------\n";
    std::cout << "-------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (modeling->nz - 1)*modeling->dz << 
                                  ", x = " << (modeling->nx - 1)*modeling->dx << ") m\n\n";

    std::cout << "Running " << keyword << " " << current << " of " << total << " in total\n\n";

    std::cout << "Current " << keyword << " position: (z = " << zpos << 
                                                    ", x = " << xpos << ") m\n\n";
    
    std::cout << current_operation << "\n";                                                   
}

void Migration::adjoint_convolution()
{
    for (int tId = 0; tId < nt; tId++)
        time_trace[tId] = (double)seismic[tId + spreadId*nt];

    fftw_execute(trace_forward_plan);

    for (int fId = 0; fId < nfft; fId++)
    {
        double a_re = freq_trace[fId][0];
        double a_im = freq_trace[fId][1];

        double b_re = freq_wavelet[fId][0];
        double b_im = freq_wavelet[fId][1];

        freq_trace[fId][0] = a_re * b_re + a_im * b_im;  
        freq_trace[fId][1] = a_im * b_re - a_re * b_im;  
    }

    fftw_execute(trace_inverse_plan);

    for (int tId = nw/2 + nw/16; tId < nt; tId++)
        h_data[tId] = (float)time_trace[tId - nw/2 - nw/16] / nfft;
}

void Migration::forward_convolution()
{


}
