# include "IDLSKDM.cuh"

void IDLSKDM::image_building()
{
    // set_rec_travel_times();        
    // set_src_travel_times();
    // prepare_convolution();

    // prepare to implement least-squares migration




}

void IDLSKDM::forward()
{
    // for (int cmpId = 0; cmpId < nCMP; cmpId++)
    // {
    //     auto cmp = CMP.find(xCMP[cmpId]); 
        
    //     int traceId = 0;

    //     for (auto &pair : cmp->second.srcId_recId)         
    //     {
    //         cudaMemset(d_data, 0.0f, nt * sizeof(float));

    //         cudaMemcpy(d_Ts, h_Ts + pair.srcId * modeling->matsize, modeling->matsize * sizeof(float), cudaMemcpyHostToDevice);
    //         cudaMemcpy(d_Tr, h_Tr + pair.recId * modeling->matsize, modeling->matsize * sizeof(float), cudaMemcpyHostToDevice);

    //         forward_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_ID_m, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);

    //         cudaMemcpy(trace_in, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

    //         # pragma omp parallel for
    //         for (int tId = 0; tId < nt; tId++)
    //             time_trace[tId] = (double)trace_in[tId];

    //         fftw_execute(trace_forward_plan);

    //         # pragma omp parallel for
    //         for (int fId = 0; fId < nfft; fId++)
    //         {
    //             double a_re = freq_trace[fId][0];
    //             double a_im = freq_trace[fId][1];

    //             double b_re = freq_wavelet[fId][0];
    //             double b_im = freq_wavelet[fId][1];

    //             freq_result[fId][0] = a_re * b_re - a_im * b_im;
    //             freq_result[fId][1] = a_re * b_im + a_im * b_re;
    //         }
    
    //         fftw_execute(result_adjoint_plan);

    //         # pragma omp paralle for
    //         for (int tId = 0; tId < nt; tId++)
    //             trace_out[tId] = (float)time_result[tId + nw/2] / nfft;
        
    //         # pragma omp parallel for
    //         for (int tId = 0; tId < nt; tId++)
    //             cal_data[tId + traceId * nt + cmp->second.partialSum] = trace_out[tId];        

    //         ++traceId;
    //     }
    // }       
}

void IDLSKDM::adjoint()
{
    // for (int cmpId = 0; cmpId < nCMP; cmpId++)
    // {
    //     auto cmp = CMP.find(xCMP[cmpId]); 
        
    //     int traceId = 0;

    //     for (auto &pair : cmp->second.srcId_recId)         
    //     {
    //         # pragma omp parallel for
    //         for (int tId = 0; tId < nt; tId++)
    //             time_trace[tId] = (double)obs_data[tId + traceId * nt + cmp->second.partialSum];

    //         fftw_execute(trace_forward_plan);

    //         # pragma omp parallel for
    //         for (int fId = 0; fId < nfft; fId++)
    //         {
    //             double a_re = freq_trace[fId][0];
    //             double a_im = freq_trace[fId][1];

    //             double b_re = freq_wavelet[fId][0];
    //             double b_im = freq_wavelet[fId][1];

    //             freq_result[fId][0] = a_re * b_re + a_im * b_im;  
    //             freq_result[fId][1] = a_im * b_re - a_re * b_im;  
    //         }
    
    //         fftw_execute(result_adjoint_plan);

    //         # pragma omp paralle for
    //         for (int tId = nw/2; tId < nt; tId++)
    //             trace_out[tId] = (float)time_result[tId - nw/2] / nfft;

    //         cudaMemcpy(d_data, trace_out, nt * sizeof(float), cudaMemcpyHostToDevice);

    //         cudaMemcpy(d_Ts, h_Ts + pair.srcId * modeling->matsize, modeling->matsize * sizeof(float), cudaMemcpyHostToDevice);
    //         cudaMemcpy(d_Tr, h_Tr + pair.recId * modeling->matsize, modeling->matsize * sizeof(float), cudaMemcpyHostToDevice);

    //         adjoint_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_ID_m, modeling->nxx, modeling->nzz, modeling->nb, nt, dt, modeling->dx, modeling->dz);

    //         ++traceId;
    //     }
    // }       
}

void IDLSKDM::export_outputs()
{


}

// __global__ void forward_kernel(float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
// {
//     int index = blockIdx.x * blockDim.x + threadIdx.x;

//     int i = (int)(index % nzz);
//     int j = (int)(index / nzz);

//     if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
//     {
//         float T = Ts[index] + Tr[index]; 
//         int tId = __float2int_rd(T / dt);
//         if (tId < nt) atomicAdd(&data[tId], m[index]);                    
//     }    
// }

// __global__ void adjoint_kernel(float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
// {
//     int index = blockIdx.x * blockDim.x + threadIdx.x;

//     int i = (int)(index % nzz);
//     int j = (int)(index / nzz);

//     if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
//     {
//         float T = Ts[index] + Tr[index]; 
//         int tId = __float2int_rd(T / dt);
//         if (tId < nt) atomicAdd(&m[index], data[tId]);
//     }    
// }
