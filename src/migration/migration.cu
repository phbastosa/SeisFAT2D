# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));
    da = std::stof(catch_parameter("mig_angle_spacing", parameters));

    max_angle = std::stof(catch_parameter("mig_max_angle", parameters));
    max_it = std::stoi(catch_parameter("mig_max_iteration", parameters));
    residuo_folder = catch_parameter("mig_residuo_folder", parameters);

    fmax = std::stof(catch_parameter("max_frequency", parameters));

    aperture = std::stof(catch_parameter("mig_aperture", parameters));

    max_offset = std::stof(catch_parameter("mig_max_offset", parameters));

    input_data_folder = catch_parameter("mig_data_folder", parameters);
    input_data_prefix = catch_parameter("mig_data_prefix", parameters);
    
    tables_folder = catch_parameter("mig_tables_folder", parameters);
    seismic_folder = catch_parameter("mig_seismic_folder", parameters);

    anisotropy = str2bool(catch_parameter("anisotropy", parameters));

    if (anisotropy) modeling = new Eikonal_ANI(); 
               else modeling = new Eikonal_ISO();

    modeling->parameters = parameters;
    
    modeling->set_geometry();

    set_interpolation();
    set_anisotropy();
    set_slowness();
    set_wavelet();
    set_gathers();

    modeling->set_eikonal();
    modeling->set_conditions();
    
    set_migration();

    h_data = new float[nt]();    
    h_model = new float[m_samples]();

    h_Ts = new float[modeling->matsize]();
    h_Tr = new float[modeling->matsize]();
    
    seismic = new float[nt*modeling->geometry->nrec]();

    cudaMalloc((void**)&(d_data), nt*sizeof(float));    
    cudaMalloc((void**)&(d_model), m_samples*sizeof(float));
    
    cudaMalloc((void**)&(d_Ts), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_Tr), modeling->matsize*sizeof(float));
}

void Migration::set_interpolation()
{
    old_nx = std::stoi(catch_parameter("x_samples", parameters));
    old_nz = std::stoi(catch_parameter("z_samples", parameters));

    old_dx = std::stof(catch_parameter("x_spacing", parameters));
    old_dz = std::stof(catch_parameter("z_spacing", parameters));

    new_dx = std::stof(catch_parameter("cubic_dx", parameters)); 
    new_dz = std::stof(catch_parameter("cubic_dz", parameters)); 

    new_nx = (int)((old_nx-1)*old_dx/new_dx)+1;
    new_nz = (int)((old_nz-1)*old_dz/new_dz)+1;

    old_nPoints = old_nx*old_nz;
    new_nPoints = new_nx*new_nz;

    modeling->nb = 3;

    modeling->dx = new_dx;
    modeling->dz = new_dz;

    modeling->nx = new_nx;
    modeling->nz = new_nz;

    modeling->nxx = modeling->nx + 2*modeling->nb;
    modeling->nzz = modeling->nz + 2*modeling->nb;

    modeling->nPoints = new_nPoints;    
    modeling->matsize = modeling->nxx*modeling->nzz;

    modeling->nBlocks = (int)((modeling->matsize + NTHREADS - 1) / NTHREADS);

    nBlocks = (int)((old_nPoints + NTHREADS - 1) / NTHREADS);
}

void Migration::set_anisotropy()
{
    if (anisotropy)
    {
        modeling->Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

        float * Cin = new float[old_nPoints]();
        float * Cout = new float[new_nPoints]();

        std::vector<std::string> targets = {"C11", "C13", "C15", "C33", "C35", "C55"}; 

        for (auto Cij : targets)
        {
            import_binary_float(modeling->Cijkl_folder + Cij + ".bin", Cin, old_nPoints);
            perform_cubic(Cin, Cout);
            export_binary_float(modeling->Cijkl_folder + "cubic_" + Cij + ".bin", Cout, new_nPoints); 
        }

        modeling->Cijkl_folder += "cubic_";
   
        delete[] Cin;
        delete[] Cout;
    }
}

void Migration::set_slowness()
{
    float * old_s = new float[old_nPoints]();
    float * new_s = new float[new_nPoints]();

    std::string slowness_file = catch_parameter("slowness_file", parameters);

    import_binary_float(slowness_file, old_s, old_nPoints);

    perform_cubic(old_s, new_s);

    modeling->S = new float[modeling->matsize]();

    modeling->expand_boundary(new_s, modeling->S);

    delete[] old_s;
    delete[] new_s;
}

void Migration::perform_cubic(float * input, float * output)
{
    float P[4][4];

    if ((new_nx == old_nx) && (new_nz == old_nz))
    {
        std::swap(output, input);
    }
    else
    {
        for (int index = 0; index < new_nPoints; index++)
        {
            int new_i = (int) (index % new_nz);    
            int new_j = (int) (index / new_nz); 
            
            float z = (float)(new_i) * new_dz;
            float x = (float)(new_j) * new_dx;
            
            float z0 = floorf(z / old_dz) * old_dz;
            float x0 = floorf(x / old_dx) * old_dx;
            
            float z1 = floorf(z / old_dz) * old_dz + old_dz;
            float x1 = floorf(x / old_dx) * old_dx + old_dx;

            float zd = (z - z0) / (z1 - z0);
            float xd = (x - x0) / (x1 - x0);

            int old_i = (int)(z / old_dz); 
            int old_j = (int)(x / old_dx);   

            if ((new_i > 0) && (new_i < new_nz-1) && (new_j > 0) && (new_j < new_nx-1))
            {
                for (int pIdx = 0; pIdx < 4; pIdx++)
                    for (int pIdz = 0; pIdz < 4; pIdz++)
                        P[pIdx][pIdz] = input[(old_i + pIdz - 1) + (old_j + pIdx - 1)*old_nz];

                output[new_i + new_j*new_nz] = cubic2d(P, xd, zd);
            }        
        }

        for (int i = 0; i < new_nz; i++)
        {
            int beg = 0;
            int end = new_nx-1;

            output[i + beg*new_nz] = output[i + (beg+1)*new_nz];
            output[i + end*new_nz] = output[i + (end-1)*new_nz];
        }

        for (int j = 0; j < new_nx; j++)
        {
            int beg = 0;
            int end = new_nz-1;

            output[beg + j*new_nz] = output[(beg+1) + j*new_nz];
            output[end + j*new_nz] = output[(end-1) + j*new_nz];
        }
    }
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
    nang = (int)(max_angle / da) + 1;

    ds = modeling->geometry->xsrc[1] - modeling->geometry->xsrc[0];
    dr = modeling->geometry->xrec[1] - modeling->geometry->xrec[0];

    dCMP = 0.5f*min(ds,dr);

    minCMP = 1e6f; minCMP = 1e6f;

    for (int srcId = 0; srcId < modeling->geometry->nsrc; srcId++)
    {
        for (int recId = 0; recId < modeling->geometry->nrec; recId++)
        {
            float CMP = 0.5f*(modeling->geometry->xsrc[srcId] + 
                              modeling->geometry->xrec[recId]);

            minCMP = (CMP < minCMP) ? CMP : minCMP;
            maxCMP = (CMP > maxCMP) ? CMP : maxCMP;
        }
    }

    nCMP = (int)((maxCMP - minCMP) / dCMP) + 1;
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
}

void Migration::set_src_travel_times()
{
    set_src_domain();

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        set_current_src();
    
        current_operation = "Computing " + keyword + " travel time matrices";

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(h_Ts, modeling->d_T, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    }
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
}

void Migration::set_rec_travel_times()
{
    set_rec_domain();
    
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        set_current_rec();

        current_operation = "Computing " + keyword + " travel time matrices";

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(h_Tr, modeling->d_T, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    }
}

void Migration::prepare_convolution()
{
    nfft = nextpow2(nt + nw - 1);
    nfreq = nfft/2 + 1;

    time_trace = (double *) fftw_malloc(nfft*sizeof(double));
    time_wavelet = (double *) fftw_malloc(nfft*sizeof(double));

    freq_trace = (fftw_complex *) fftw_malloc(nfreq*sizeof(fftw_complex));
    freq_wavelet = (fftw_complex *) fftw_malloc(nfreq*sizeof(fftw_complex));

    trace_forward_plan = fftw_plan_dft_r2c_1d(nfft, time_trace, freq_trace, FFTW_ESTIMATE);
    trace_inverse_plan = fftw_plan_dft_c2r_1d(nfft, freq_trace, time_trace, FFTW_ESTIMATE);
    wavelet_forward_plan = fftw_plan_dft_r2c_1d(nfft, time_wavelet, freq_wavelet, FFTW_ESTIMATE);

    for (int tId = 0; tId < nfft; tId++)
    {
        time_trace[tId] = 0.0;
        time_wavelet[tId] = (tId < nw) ? (double)wavelet[tId] : 0.0;
    }

    fftw_execute(wavelet_forward_plan);
}

void Migration::adjoint_convolution()
{
    for (int tId = 0; tId < nfft; tId++)
    {
        if (tId < nt)
        {
            time_trace[tId] = (double)h_data[tId];
            h_data[tId] = 0.0f;
        }
        else 
            time_trace[tId] = 0.0;    
    }
    
    fftw_execute(trace_forward_plan);

    for (int fId = 0; fId < nfreq; fId++)
    {
        double a_re = freq_trace[fId][0];
        double a_im = freq_trace[fId][1];

        double b_re = freq_wavelet[fId][0];
        double b_im = freq_wavelet[fId][1];

        freq_trace[fId][0] = a_re * b_re + a_im * b_im;  
        freq_trace[fId][1] = a_im * b_re - a_re * b_im;  
    }

    fftw_execute(trace_inverse_plan);

    for (int tId = nw/2 + nw/10; tId < nt; tId++)
        h_data[tId] = (float)(time_trace[tId - nw/2 - nw/10]);
}

void Migration::forward_convolution()
{
    for (int tId = 0; tId < nfft; tId++)
    {
        if (tId < nt)
        {
            time_trace[tId] = (double)h_data[tId];
            h_data[tId] = 0.0f;
        }
        else 
            time_trace[tId] = 0.0;    
    }

    fftw_execute(trace_forward_plan);

    for (int fId = 0; fId < nfreq; fId++)
    {
        double a_re = freq_trace[fId][0];
        double a_im = freq_trace[fId][1];

        double b_re = freq_wavelet[fId][0];
        double b_im = freq_wavelet[fId][1];

        freq_trace[fId][0] = a_re * b_re - a_im * b_im;
        freq_trace[fId][1] = a_re * b_im + a_im * b_re;
    }

    fftw_execute(trace_inverse_plan);

    for (int tId = 0; tId < nt; tId++)
        h_data[tId] = (float)(time_trace[tId + nw/2 + nw/10]);
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

// void Migration::dot_product_test()
// {
//     set_src_travel_times();
//     set_rec_travel_times();
//     prepare_convolution();
    
//     auto trash = system("clear");
//     std::cout << "Computing dot product test!\n\n";

//    d_samples = nt*modeling->geometry->nrec*modeling->geometry->nsrc;

//     d1 = new float[d_samples]();
//     d2 = new float[d_samples]();

//     m1 = new float[m_samples]();
//     m2 = new float[m_samples]();

//     int minValue =-100;
//     int maxValue = 100;

//     std::mt19937 prng(std::random_device{}());
//     std::uniform_int_distribution<int> dist(minValue, maxValue);

//     for (int mId = 0; mId < m_samples; mId++)
//         m1[mId] = dist(prng);

//     for (int dId = 0; dId < d_samples; dId++)
//         d2[dId] = dist(prng);
    
//     cudaMemcpy(d_model, m1, m_samples*sizeof(float), cudaMemcpyHostToDevice);

//     for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
//     {     
//         import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
//         cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
        
//         int spreadId = 0;
        
//         for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
//         {
//             cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              
//             cmp = 0.5f*(modeling->geometry->xsrc[modeling->srcId] + 
//                         modeling->geometry->xrec[modeling->recId]);

//             import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
//             cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

//             cudaMemset(d_data, 0.0f, nt * sizeof(float));

//             perform_forward();

//             cudaMemcpy(h_data, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

//             forward_convolution();

//             for (int tId = 0; tId < nt; tId++)
//             {
//                 int index = tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt;
                
//                 d1[index] = h_data[tId];    
//             }                
        
//             ++spreadId;
//         }
//     }

//     cudaMemset(d_model, 0.0f, m_samples*sizeof(float));

//     for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
//     { 
//         modeling->sx = modeling->geometry->xsrc[modeling->srcId];

//         import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
//         cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

//         int spreadId = 0;
        
//         for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
//         {
//             cmp = 0.5f*(modeling->sx + modeling->geometry->xrec[modeling->recId]);
//             cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              

//             import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
//             cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

//             for (int tId = 0; tId < nt; tId++)
//                 h_data[tId] = d2[tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt];

//             adjoint_convolution();

//             cudaMemcpy(d_data, h_data, nt * sizeof(float), cudaMemcpyHostToDevice);
            
//             perform_adjoint();

//             ++spreadId;
//         }
//     }

//     cudaMemcpy(m2, d_model, m_samples*sizeof(float), cudaMemcpyDeviceToHost);        

//     double dot_m = 0.0;
//     double dot_d = 0.0;

//     for (int mId = 0; mId < m_samples; mId++)
//         dot_m += m1[mId]*m2[mId];
    
//     for (int dId = 0; dId < d_samples; dId++)
//         dot_d += d1[dId]*d2[dId];

//     double r = (dot_d - dot_m) / (dot_d + dot_m);
    
//     std::cout << "<m1, m2> = " << dot_m << std::endl;
//     std::cout << "<d1, d2> = " << dot_d << std::endl; 
//     std::cout << "residuo = " << r << std::endl;
// }

__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dt, int nt, 
                                            float old_dx, float old_dz, float new_dx, float new_dz, int old_nx, int old_nz, 
                                            int new_nxx, int new_nzz, int nb, float aperture, float cmp)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % old_nz);
    int j = (int)(index / old_nz);

    float z = __int2float_rd(i)*old_dz;
    float x = __int2float_rd(j)*old_dx;

    float z0 = __float2int_rd(z / new_dz) * new_dz;
    float x0 = __float2int_rd(x / new_dx) * new_dx;
    
    float z1 = __float2int_rd(z / new_dz) * new_dz + new_dz;
    float x1 = __float2int_rd(x / new_dx) * new_dx + new_dx;

    float zd = (z - z0) / (z1 - z0);
    float xd = (x - x0) / (x1 - x0);

    int new_i = __float2int_rd(z / new_dz) + nb; 
    int new_j = __float2int_rd(x / new_dx) + nb;   

    const int n = 4;

    float cS[n][n];
    float cTs[n][n]; 
    float cTr[n][n]; 

    if ((new_i >= nb) && (new_i < new_nzz-nb) && 
        (new_j >= nb) && (new_j < new_nxx-nb))
    {
        for (int pIdx = 0; pIdx < n; pIdx++)
        {
            for (int pIdz = 0; pIdz < n; pIdz++)
            {
                int targetId = (new_i+pIdz-1) + (new_j+pIdx-1)*new_nzz; 

                cS[pIdx][pIdz] = S[targetId];
                cTs[pIdx][pIdz] = Ts[targetId];
                cTr[pIdx][pIdz] = Tr[targetId];
            }
        }

        float cubic_S = d_cubic2d(cS, xd, zd);
        float cubic_Ts = d_cubic2d(cTs, xd, zd);
        float cubic_Tr = d_cubic2d(cTr, xd, zd);

        float T = cubic_Ts + cubic_Tr; 
        int tId = __float2int_rd(T / dt);

        if (tId < nt)
        {
            float dTs_dz = 0.5f*(Ts[(new_i+1) + new_j*new_nzz] - Ts[(new_i-1) + new_j*new_nzz]) / new_dz;            
            float dTs_dx = 0.5f*(Ts[new_i + (new_j+1)*new_nzz] - Ts[new_i + (new_j-1)*new_nzz]) / new_dx;

            float d2Ts_dz2 = (Ts[(new_i+1) + new_j*new_nzz] - 2.0f*Ts[new_i + new_j*new_nzz] + Ts[(new_i-1) + new_j*new_nzz]) / (new_dz*new_dz); 
            float d2Ts_dx2 = (Ts[new_i + (new_j+1)*new_nzz] - 2.0f*Ts[new_i + new_j*new_nzz] + Ts[new_i + (new_j-1)*new_nzz]) / (new_dx*new_dx); 
           
            float d2Ts_dxdz = (Ts[(new_i+1) + (new_j+1)*new_nzz] - Ts[(new_i-1) + (new_j+1)*new_nzz] - 
                               Ts[(new_i+1) + (new_j-1)*new_nzz] + Ts[(new_i-1) + (new_j-1)*new_nzz]) / (4.0f*new_dx*new_dz);

            float dTr_dz = 0.5f*(Tr[(new_i+1) + new_j*new_nzz] - Tr[(new_i-1) + new_j*new_nzz]) / new_dz;            
            float dTr_dx = 0.5f*(Tr[new_i + (new_j+1)*new_nzz] - Tr[new_i + (new_j-1)*new_nzz]) / new_dx;

            float d2Tr_dz2 = (Tr[(new_i+1) + new_j*new_nzz] - 2.0f*Tr[new_i + new_j*new_nzz] + Tr[(new_i-1) + new_j*new_nzz]) / (new_dz*new_dz); 
            float d2Tr_dx2 = (Tr[new_i + (new_j+1)*new_nzz] - 2.0f*Tr[new_i + new_j*new_nzz] + Tr[new_i + (new_j-1)*new_nzz]) / (new_dx*new_dx); 
           
            float d2Tr_dxdz = (Tr[(new_i+1) + (new_j+1)*new_nzz] - Tr[(new_i-1) + (new_j+1)*new_nzz] - 
                               Tr[(new_i+1) + (new_j-1)*new_nzz] + Tr[(new_i-1) + (new_j-1)*new_nzz]) / (4.0f*new_dx*new_dz);
 
            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;

            float a1 = d2Ts_dx2; float a2 = d2Tr_dx2;
            float b1 = d2Ts_dxdz; float b2 = d2Tr_dxdz;
            float c1 = d2Ts_dz2; float c2 = d2Tr_dz2;

            float detHs = a1*c1 - b1*b1;
            float detHr = a2*c2 - b2*b2;

            detHs = fabsf(detHs) + EPS;
            detHr = fabsf(detHr) + EPS;

            float Jterm = sqrtf(detHs * detHr) * cubic_S / (norm_Ts * norm_Tr);    

            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float obliquity = sqrt(max(0.0f, cos_s) * max(0.0f, cos_r));

            float R_s = max(cubic_Ts / cubic_S, EPS);
            float R_r = max(cubic_Tr / cubic_S, EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta_s = acosf(min(1.0f,max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));
            float theta_r = acosf(min(1.0f,max(-1.0f, ux_r*nx_norm + uz_r*nz_norm)));
            float theta = 0.5f*(theta_s + theta_r);
            float R = 1.0f + 0.2f*cosf(theta); 

            float scale_const = 1.0f / (2.0f * M_PI);

            float sigma = tanf(aperture * M_PI / 180.0f)*z;        

            float taper = expf(-0.5f*((x - cmp) / (sigma + EPS))*
                                     ((x - cmp) / (sigma + EPS)));

            float weights = scale_const * Jterm * obliquity * G * R * taper;            
    
            atomicAdd(&model[index], weights * data[tId]);
        }
    }    
}

__global__ void image_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dt, int nt, 
                                            float old_dx, float old_dz, float new_dx, float new_dz, int old_nx, int old_nz, 
                                            int new_nxx, int new_nzz, int nb, float aperture, float cmp)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % old_nz);
    int j = (int)(index / old_nz);

    float z = __int2float_rd(i)*old_dz;
    float x = __int2float_rd(j)*old_dx;

    float z0 = __float2int_rd(z / new_dz) * new_dz;
    float x0 = __float2int_rd(x / new_dx) * new_dx;
    
    float z1 = __float2int_rd(z / new_dz) * new_dz + new_dz;
    float x1 = __float2int_rd(x / new_dx) * new_dx + new_dx;

    float zd = (z - z0) / (z1 - z0);
    float xd = (x - x0) / (x1 - x0);

    int new_i = __float2int_rd(z / new_dz) + nb; 
    int new_j = __float2int_rd(x / new_dx) + nb;   

    const int n = 4;

    float cS[n][n];
    float cTs[n][n]; 
    float cTr[n][n]; 

    if ((new_i >= nb) && (new_i < new_nzz-nb) && 
        (new_j >= nb) && (new_j < new_nxx-nb))
    {
        for (int pIdx = 0; pIdx < n; pIdx++)
        {
            for (int pIdz = 0; pIdz < n; pIdz++)
            {
                int targetId = (new_i+pIdz-1) + (new_j+pIdx-1)*new_nzz; 

                cS[pIdx][pIdz] = S[targetId];
                cTs[pIdx][pIdz] = Ts[targetId];
                cTr[pIdx][pIdz] = Tr[targetId];
            }
        }

        float cubic_S = d_cubic2d(cS, xd, zd);
        float cubic_Ts = d_cubic2d(cTs, xd, zd);
        float cubic_Tr = d_cubic2d(cTr, xd, zd);

        float T = cubic_Ts + cubic_Tr; 
        int tId = __float2int_rd(T / dt);

        if (tId < nt)
        {
            float dTs_dz = 0.5f*(Ts[(new_i+1) + new_j*new_nzz] - Ts[(new_i-1) + new_j*new_nzz]) / new_dz;            
            float dTs_dx = 0.5f*(Ts[new_i + (new_j+1)*new_nzz] - Ts[new_i + (new_j-1)*new_nzz]) / new_dx;

            float d2Ts_dz2 = (Ts[(new_i+1) + new_j*new_nzz] - 2.0f*Ts[new_i + new_j*new_nzz] + Ts[(new_i-1) + new_j*new_nzz]) / (new_dz*new_dz); 
            float d2Ts_dx2 = (Ts[new_i + (new_j+1)*new_nzz] - 2.0f*Ts[new_i + new_j*new_nzz] + Ts[new_i + (new_j-1)*new_nzz]) / (new_dx*new_dx); 
           
            float d2Ts_dxdz = (Ts[(new_i+1) + (new_j+1)*new_nzz] - Ts[(new_i-1) + (new_j+1)*new_nzz] - 
                               Ts[(new_i+1) + (new_j-1)*new_nzz] + Ts[(new_i-1) + (new_j-1)*new_nzz]) / (4.0f*new_dx*new_dz);

            float dTr_dz = 0.5f*(Tr[(new_i+1) + new_j*new_nzz] - Tr[(new_i-1) + new_j*new_nzz]) / new_dz;            
            float dTr_dx = 0.5f*(Tr[new_i + (new_j+1)*new_nzz] - Tr[new_i + (new_j-1)*new_nzz]) / new_dx;

            float d2Tr_dz2 = (Tr[(new_i+1) + new_j*new_nzz] - 2.0f*Tr[new_i + new_j*new_nzz] + Tr[(new_i-1) + new_j*new_nzz]) / (new_dz*new_dz); 
            float d2Tr_dx2 = (Tr[new_i + (new_j+1)*new_nzz] - 2.0f*Tr[new_i + new_j*new_nzz] + Tr[new_i + (new_j-1)*new_nzz]) / (new_dx*new_dx); 
           
            float d2Tr_dxdz = (Tr[(new_i+1) + (new_j+1)*new_nzz] - Tr[(new_i-1) + (new_j+1)*new_nzz] - 
                               Tr[(new_i+1) + (new_j-1)*new_nzz] + Tr[(new_i-1) + (new_j-1)*new_nzz]) / (4.0f*new_dx*new_dz);
 
            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;

            float a1 = d2Ts_dx2; float a2 = d2Tr_dx2;
            float b1 = d2Ts_dxdz; float b2 = d2Tr_dxdz;
            float c1 = d2Ts_dz2; float c2 = d2Tr_dz2;

            float detHs = a1*c1 - b1*b1;
            float detHr = a2*c2 - b2*b2;

            detHs = fabsf(detHs) + EPS;
            detHr = fabsf(detHr) + EPS;

            float Jterm = sqrtf(detHs * detHr) * cubic_S / (norm_Ts * norm_Tr);    

            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float obliquity = sqrt(max(0.0f, cos_s) * max(0.0f, cos_r));

            float R_s = max(cubic_Ts / cubic_S, EPS);
            float R_r = max(cubic_Tr / cubic_S, EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta_s = acosf(min(1.0f,max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));
            float theta_r = acosf(min(1.0f,max(-1.0f, ux_r*nx_norm + uz_r*nz_norm)));
            float theta = 0.5f*(theta_s + theta_r);
            float R = 1.0f + 0.2f*cosf(theta); 

            float scale_const = 1.0f / (2.0f * M_PI);

            float sigma = tanf(aperture * M_PI / 180.0f)*z;        

            float taper = expf(-0.5f*((x - cmp) / (sigma + EPS))*
                                     ((x - cmp) / (sigma + EPS)));

            float weights = scale_const * Jterm * obliquity * G * R * taper;            

            atomicAdd(&data[tId], weights * model[index]);            
        }
    }    
}

__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        
        int tId = __float2int_rd(T / dt);

        int nz = (nzz - 2*nb);
        
        if (tId < nt) 
        {
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

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

            float R_s = max(Ts[index] / S[index], EPS);
            float R_r = max(Tr[index] / S[index], EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

            float reflection_angle = 0.5f*acos((dTs_dx*dTr_dx + dTs_dz*dTr_dz) / (norm_Ts * norm_Tr));

            float ang = 180.0f*reflection_angle / M_PI;
            int aId = __float2int_rd(ang / da);  

            int mId = (i - nb) + aId*nz + cmpId*na*nz;
            
            if ((aId >= 0) && (aId < na))             
                atomicAdd(&model[mId], weights * data[tId]);
        }
    }   
}

__global__ void angle_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        
        int tId = __float2int_rd(T / dt);

        int nz = (nzz - 2*nb);
        
        if (tId < nt) 
        {
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

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

            float R_s = max(Ts[index] / S[index], EPS);
            float R_r = max(Tr[index] / S[index], EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

            float reflection_angle = 0.5f*acos((dTs_dx*dTr_dx + dTs_dz*dTr_dz) / (norm_Ts * norm_Tr));

            float ang = 180.0f*reflection_angle / M_PI;
            int aId = __float2int_rd(ang / da);  

            int mId = (i - nb) + aId*nz + cmpId*na*nz;
            
            if ((aId >= 0) && (aId < na))             
                atomicAdd(&data[tId], weights * model[mId]);
        }
    }   
}

__device__ float d_cubic1d(float P[4], float dx)
{
    return P[1] + 0.5f*dx*(P[2] - P[0] + dx*(2.0f*P[0] - 5.0f*P[1] + 4.0f*P[2] - P[3] + dx*(3.0f*(P[1] - P[2]) + P[3] - P[0])));
}

__device__ float d_cubic2d(float P[4][4], float dx, float dy)
{    
    float p[4];
    p[0] = d_cubic1d(P[0], dy);
    p[1] = d_cubic1d(P[1], dy);
    p[2] = d_cubic1d(P[2], dy);
    p[3] = d_cubic1d(P[3], dy);    
    return d_cubic1d(p, dx);
}
