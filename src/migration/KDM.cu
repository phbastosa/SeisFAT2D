# include "KDM.cuh"

void KDM::image_building()
{
    set_src_travel_times();
    set_rec_travel_times();
    prepare_convolution();
    prepare_components();

    cudaMemset(d_image, 0.0f, modeling->matsize*sizeof(float));

    keyword = "source";
    
    total = std::to_string(modeling->geometry->nsrc); 

    std::vector<int> localCount(nCMP, 0);    

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
        
        import_binary_float(data_path, seismic, nt*modeling->max_spread);
     
        import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);

        cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
        
        float sx = modeling->geometry->xsrc[modeling->srcId];
        float sz = modeling->geometry->zsrc[modeling->srcId];

        current = std::to_string(modeling->srcId+1);
        
        xpos = format1Decimal(sx);
        zpos = format1Decimal(sz);

        current_operation = "Kirchhoff depth migration: adjoint operator";

        show_information();

        int spreadId = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {            
            float rx = modeling->geometry->xrec[modeling->recId];
            
            float cmp = 0.5f*(sx + rx);

            import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);

            cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

            for (int tId = 0; tId < nt; tId++)
                time_trace[tId] = (double)seismic[tId + spreadId*nt];

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

            cudaMemcpy(d_data, trace_out, nt*sizeof(float), cudaMemcpyHostToDevice);
                
            cross_correlation<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_trace, d_angle, d_image, cmp, aperture, nTraces, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);

            cudaMemcpy(h_trace, d_trace, modeling->nz*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_angle, d_angle, modeling->nz*sizeof(float), cudaMemcpyDeviceToHost);

            int cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;             
            int traceId = localCount[cmpId]++;                         
            int gatherId = pSUM[cmpId] + traceId;    

            for (int zId = 0; zId < modeling->nz; zId++)
            {
                ODCIG[zId + gatherId*modeling->nz] = h_trace[zId];

                float angle = 180.0f*h_angle[zId]/M_PI/da;
                int angleId = static_cast<int>(angle);

                if ((angleId >= 0) && (angleId < nang))  
                    ADCIG[zId + angleId*modeling->nz + cmpId*nang*modeling->nz] = h_trace[zId];
            }

            ++spreadId;
        }
    }
}

void KDM::prepare_components()
{
    IMAGE = new float[modeling->nPoints]();
    ODCIG = new float[modeling->nz * nTraces]();
    ADCIG = new float[modeling->nz * nCMP*nang]();
}

void KDM::export_outputs()
{
    cudaMemcpy(h_image, d_image, modeling->matsize * sizeof(float), cudaMemcpyDeviceToHost);
    
    modeling->reduce_boundary(h_image, IMAGE);
    
    export_binary_float("image.bin", IMAGE, modeling->nPoints);
    export_binary_float("ODCIG.bin", ODCIG, modeling->nz*nTraces);
    export_binary_float("ADCIG.bin", ADCIG, modeling->nz*nang*nCMP);
}

__global__ void cross_correlation(float * S, float * Ts, float * Tr, float * data, float * trace, float * angle, float * image, float cmp, float aperture, int nTraces, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
    
        int tId = __float2int_rd(T / dt);
        int cmpId = __float2int_rd(cmp / dx);

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

            float sigma = tanf(aperture * M_PI / 180.0f)*(i-nb)*dz;
            float taper = expf(-0.5f*powf(((j-nb)*dx - cmp) / (sigma + eps), 2.0f));

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R * taper;            
            
            float reflection_angle = 0.5f*acos((dTs_dx*dTr_dx + dTs_dz*dTr_dz) / (norm_Ts * norm_Tr));
            
            float seismic_amplitude = weights * data[tId];    
            
            atomicAdd(&image[index], seismic_amplitude);

            if (j == cmpId + nb)
            {
                angle[i - nb] = reflection_angle;    
                trace[i - nb] = seismic_amplitude;
            }                        
        }    
    }
}
