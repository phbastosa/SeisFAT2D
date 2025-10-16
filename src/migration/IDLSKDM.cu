# include "IDLSKDM.cuh"

void IDLSKDM::image_building()
{
    set_src_travel_times();
    set_rec_travel_times();
    prepare_convolution();
    prepare_components();
    set_seismic_data();

    set_initial_model();

    // while (true)
    // {
        check_convergence();

    //     if (converged) break;

        compute_gradient();
        compute_direction();
        compute_stepLength();
        update_reflectivity();
    // }
}

void IDLSKDM::prepare_components()
{
    dobs = new float[nt*nTraces]();
    dcal = new float[nt*nTraces]();
    dres = new float[nt*nTraces]();
    
    model = new float[modeling->matsize]();
    
    gradient_p = new float[modeling->matsize]();
    gradient_m = new float[modeling->matsize]();
    
    direction = new float[modeling->matsize]();

    IMAGE = new float[modeling->nPoints]();

    iteration = 0;
}

void IDLSKDM::set_seismic_data()
{
    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
        
        import_binary_float(data_path, seismic, nt*modeling->max_spread);

        for (int tId = 0; tId < nt; tId++)
        {
            for (int spreadId = 0; spreadId < modeling->max_spread; spreadId++)
            {
                int index = tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt;

                dobs[index] = seismic[tId + spreadId*nt]; 
            }
        }    
    }
}

void IDLSKDM::set_initial_model()
{
    current_operation = "IDLSKDM: adjoint operator -> Generating initial model";

    adjoint(dobs, model);
}

void IDLSKDM::check_convergence()
{
    current_operation = "IDLSKDM: forward operator -> Checking convergence";

    forward(model, dcal);        

    float r = 0.0f;
    for (int index = 0; index < nt*nTraces; index++)
    {
        r += dobs[index] - dcal[index];

        dres[index] = dobs[index] - dcal[index];
    }

    residuo.push_back(sqrtf(r));
}

void IDLSKDM::compute_gradient()
{
    current_operation = "IDLSKDM: adjoint operator -> Computing gradient";

    adjoint(dres, gradient_p);
    
    float grad_norm = 0.0f;
    for (int index = 0; index < modeling->matsize; index++)
        grad_norm += gradient_p[index]*gradient_p[index];

    grad_norm = sqrtf(grad_norm);

    for (int index = 0; index < modeling->matsize; index++)
        gradient_p[index] /= grad_norm;
}

void IDLSKDM::compute_direction()
{
    beta = 0.0f;
    for (int index = 0; index < modeling->matsize; index++)
        beta += (gradient_p[index] - gradient_m[index])*gradient_p[index] / (gradient_m[index] * gradient_m[index] + 1e-6f); 

    for (int index = 0; index < modeling->matsize; index++)
        direction[index] = beta*direction[index] - gradient_p[index];
}

void IDLSKDM::compute_stepLength()
{
    current_operation = "IDLSKDM: forward operator -> Computing step length";

    forward(direction, dcal);

    float sum1 = 0.0f; 
    float sum2 = 0.0f;
    
    for (int index = 0; index < nt*nTraces; index++)
    {
        sum1 += dres[index]*dcal[index];
        sum2 += dcal[index]*dcal[index];
    }
    
    alpha = sum1 / sum2;
}

void IDLSKDM::update_reflectivity()
{
    for (int index = 0; index < modeling->matsize; index++)
    {
        model[index] = model[index] + alpha*direction[index];
        gradient_m[index] = gradient_p[index];
    }   
}

void IDLSKDM::forward(float * m, float * d)
{
    std::vector<int> localCount(nCMP, 0);    

    cudaMemcpy(d_image, m, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    keyword = "source";
    
    total = std::to_string(modeling->geometry->nsrc); 

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {     
        import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);

        cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

        float sx = modeling->geometry->xsrc[modeling->srcId];
        float sz = modeling->geometry->zsrc[modeling->srcId];

        current = std::to_string(modeling->srcId+1);
        
        xpos = format1Decimal(sx);
        zpos = format1Decimal(sz);

        show_information();
                
        int spreadId = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            cudaMemset(d_data, 0.0f, nt * sizeof(float));

            import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);

            cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

            forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_image, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);

            cudaMemcpy(trace_in, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

            for (int tId = 0; tId < nt; tId++)
                time_trace[tId] = trace_in[tId];

            fftw_execute(trace_forward_plan);

            for (int fId = 0; fId < nfft; fId++)
            {
                double a_re = freq_trace[fId][0];
                double a_im = freq_trace[fId][1];

                double b_re = freq_wavelet[fId][0];
                double b_im = freq_wavelet[fId][1];

                freq_result[fId][0] = a_re * b_re - a_im * b_im;
                freq_result[fId][1] = a_re * b_im + a_im * b_re;
            }
    
            fftw_execute(result_adjoint_plan);

            for (int tId = 0; tId < nt; tId++)
            {
                int index = tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt;
                
                d[index] = time_result[tId + nw/2] / nfft;    
            }                
        
            ++spreadId;
        }
    }
}

void IDLSKDM::adjoint(float * d, float * m)
{
    cudaMemset(d_image, 0.0f, modeling->matsize*sizeof(float));

    keyword = "source";
    
    total = std::to_string(modeling->geometry->nsrc); 

    std::vector<int> localCount(nCMP, 0);    

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    { 
        import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);

        cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

        float sx = modeling->geometry->xsrc[modeling->srcId];
        float sz = modeling->geometry->zsrc[modeling->srcId];

        current = std::to_string(modeling->srcId+1);
        
        xpos = format1Decimal(sx);
        zpos = format1Decimal(sz);

        show_information();
        
        int spreadId = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);

            cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

            for (int tId = 0; tId < nt; tId++)
            {
                int trace_extraction = tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt;
                
                time_trace[tId] = d[trace_extraction];
            }

            fftw_execute(trace_forward_plan);
    
            for (int fId = 0; fId < nfft; fId++)
            {
                double a_re = freq_trace[fId][0];
                double a_im = freq_trace[fId][1];

                double b_re = freq_wavelet[fId][0];
                double b_im = freq_wavelet[fId][1];

                freq_result[fId][0] = a_re * b_re + a_im * b_im;  
                freq_result[fId][1] = a_im * b_re - a_re * b_im;  
            }
        
            fftw_execute(result_adjoint_plan);

            for (int tId = nw/2; tId < nt; tId++)
                trace_out[tId] = (float)time_result[tId - nw/2] / nfft;

            cudaMemcpy(d_data, trace_out, nt * sizeof(float), cudaMemcpyHostToDevice);
                
            adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_image, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);

            ++spreadId;
        }
    }

    cudaMemcpy(m, d_image, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);
}

void IDLSKDM::export_outputs()
{


}

__global__ void forward_kernel(float * S, float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        int tId = __float2int_rd(T / dt);
        
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
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float absDetH = fabsf(detH) + eps;
            float J = 1.0f / sqrtf(absDetH); 

            float v = 1.0f / S[index];
            float R_s = max(Ts[index]*v, eps);
            float R_r = max(Tr[index]*v, eps);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));
            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

            atomicAdd(&data[tId], weights * m[index]);                    
        }            
    }    
}

__global__ void adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        int tId = __float2int_rd(T / dt);
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
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float absDetH = fabsf(detH) + eps;
            float J = 1.0f / sqrtf(absDetH); 

            float v = 1.0f / S[index];
            float R_s = max(Ts[index]*v, eps);
            float R_r = max(Tr[index]*v, eps);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));
            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            
                
            atomicAdd(&m[index], weights * data[tId]);
        }
    }    
}
