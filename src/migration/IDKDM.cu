# include "IDKDM.cuh"

void IDKDM::set_migration()
{
    h_model = new float[modeling->nPoints]();

    cudaMalloc((void**)&(d_model), modeling->nPoints*sizeof(float));
}

void IDKDM::initialization()
{
    cudaMemset(d_model, 0.0f, modeling->nPoints*sizeof(float));    
}

void IDKDM::perform_migration()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, modeling->nxx, modeling->nzz, nt, modeling->nb);
}

void IDKDM::export_outputs()
{
    cudaMemcpy(h_model, d_model, modeling->nPoints*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float("IDKDM.bin", h_model, modeling->nPoints);
}

__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, int nxx, int nzz, int nt, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        int nz = nzz - 2*nb;

        float T = Ts[index] + Tr[index]; 
        
        int tId = __float2int_rd(T / dt);
        
        int mId = (i - nb) + (j - nb)*nz;

        if (tId < nt) atomicAdd(&model[mId], data[tId]);
    }
}
