# include "ADKDM.cuh"

void ADKDM::set_migration()
{
    h_model = new float[modeling->nz*nCMP*nang]();
    
    cudaMalloc((void**)&(d_model), modeling->nz*nCMP*nang*sizeof(float));
}

void ADKDM::initialization()
{
    cudaMemset(d_model, 0.0f, modeling->nz*nCMP*nang*sizeof(float));    
}

void ADKDM::perform_migration()
{
   angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADKDM::export_outputs()
{
    cudaMemcpy(h_model, d_model, modeling->nz*nCMP*nang*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float("ADKDM.bin", h_model, modeling->nz*nCMP*nang);    
}

__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        int nz = (nzz - 2*nb);

        float T = Ts[index] + Tr[index]; 
        
        int tId = __float2int_rd(T / dt);
        
        if (tId < nt) 
        {
            float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
            float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

            float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
            float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

            float norm_dTs = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_dTr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;

            float reflection_angle = 0.5f*acos((dTs_dx*dTr_dx + dTs_dz*dTr_dz) / (norm_dTs * norm_dTr));

            float ang = 180.0f*reflection_angle / M_PI;
            int aId = __float2int_rd(ang / da);  
            
            int mId = (i - nb) + aId*nz + cmpId*na*nz;
            
            atomicAdd(&model[mId], data[tId]);
        }
    }   
}
