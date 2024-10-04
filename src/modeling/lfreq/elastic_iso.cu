# include "elastic_iso.cuh"

void elastic_Iso::set_properties()
{
    std::string vp_file = catch_parameter("vp_model_file", parameters);
    std::string vs_file = catch_parameter("vs_model_file", parameters);
    std::string rho_file = catch_parameter("rho_model_file", parameters);

    float * vp = new float[nPoints]();
    float * vs = new float[nPoints]();
    float * rho = new float[nPoints]();

    Vp = new float[matsize]();
    Vs = new float[matsize]();
    Rho = new float[matsize]();

    import_binary_float(vp_file, vp, nPoints);
    import_binary_float(vs_file, vs, nPoints);
    import_binary_float(rho_file, rho, nPoints);

    expand_boundary(vp, Vp);
    expand_boundary(vs, Vs);
    expand_boundary(rho, Rho);

    delete[] vp;
    delete[] vs;
    delete[] rho;
}

void elastic_Iso::set_conditions()
{
    M = new float[matsize]();
    L = new float[matsize]();
    B = new float[matsize]();
    P = new float[matsize]();

    for (int index = 0; index < matsize; index++)
    {
        M[index] = Rho[index]*Vs[index]*Vs[index];
        L[index] = Rho[index]*Vp[index]*Vp[index] - 2.0f*M[index];
        B[index] = 1.0f / Rho[index];
    }

    synthetic_data = new float[nt*max_spread]();
    cudaMalloc((void**)&(seismogram), nt*max_spread*sizeof(float));

    cudaMalloc((void**)&(d_M), matsize*sizeof(float));
    cudaMalloc((void**)&(d_L), matsize*sizeof(float));
    cudaMalloc((void**)&(d_B), matsize*sizeof(float));
    cudaMalloc((void**)&(d_P), matsize*sizeof(float));

    cudaMalloc((void**)&(d_Vx), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Vz), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Txx), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Tzz), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Txz), matsize*sizeof(float));

    cudaMemcpy(d_M, M, matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L, L, matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, matsize*sizeof(float), cudaMemcpyHostToDevice);
}

void elastic_Iso::initialization()
{
    cudaMemset(d_P, 0.0f, matsize*sizeof(float));
    cudaMemset(d_Vx, 0.0f, matsize*sizeof(float));
    cudaMemset(d_Vz, 0.0f, matsize*sizeof(float));
    cudaMemset(d_Txx, 0.0f, matsize*sizeof(float));
    cudaMemset(d_Tzz, 0.0f, matsize*sizeof(float));
    cudaMemset(d_Txz, 0.0f, matsize*sizeof(float));

    sIdx = (int)(geometry->xsrc[geometry->sInd[srcId]] / dx) + nb;
    sIdz = (int)(geometry->zsrc[geometry->sInd[srcId]] / dz) + nb;

    spread = 0;

    for (recId = geometry->iRec[srcId]; recId < geometry->fRec[srcId]; recId++)
    {
        current_xrec[spread] = (int)(geometry->xrec[recId] / dx) + nb;
        current_zrec[spread] = (int)(geometry->zrec[recId] / dz) + nb;

        spread++;
    }

    sBlocks = (int)(spread / nThreads) + 1; 

    cudaMemcpy(rIdx, current_xrec, spread*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(rIdz, current_zrec, spread*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemset(seismogram, 0.0f, nt*spread*sizeof(float));
}

void elastic_Iso::forward_solver()
{
    for (int tId = 0; tId < tlag + nt; tId++)
    {
        compute_pressure<<<nBlocks, nThreads>>>(d_Vx, d_Vz, d_Txx, d_Tzz, d_Txz, d_P, d_M, d_L, wavelet, sIdx, sIdz, tId, nt, nxx, nzz, dx, dz, dt);
        cudaDeviceSynchronize();

        compute_velocity<<<nBlocks, nThreads>>>(d_Vx, d_Vz, d_Txx, d_Tzz, d_Txz, d_B, d1D, d2D, nb, nxx, nzz, dx, dz, dt);
        cudaDeviceSynchronize();

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdz, seismogram, spread, tId, tlag, nt, nzz);     
        cudaDeviceSynchronize();
    }

    cudaMemcpy(P, d_P, matsize*sizeof(float), cudaMemcpyDeviceToHost);
    export_binary_float("P_shot_" + std::to_string(srcId+1) + ".bin", P, matsize);

    cudaMemcpy(synthetic_data, seismogram, nt*spread*sizeof(float), cudaMemcpyDeviceToHost);
    export_binary_float("seismogram_shot_" + std::to_string(srcId+1) + ".bin", synthetic_data, nt*spread);
}

__global__ void compute_pressure(float * Vx, float * Vz, float * Txx, float * Tzz, float * Txz, float * P, float * M, float * L, float * wavelet, int sIdx, int sIdz, int tId, int nt, int nxx, int nzz, float dx, float dz, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((index == 0) && (tId < nt))
    {
        Txx[sIdz + sIdx*nzz] += wavelet[tId] / (dx * dz); 
        Tzz[sIdz + sIdx*nzz] += wavelet[tId] / (dx * dz);
    }

    if ((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4)) 
    {           
        float dVx_dx = (75.0f*(Vx[i + (j-3)*nzz] - Vx[i + (j+4)*nzz]) + 
                      1029.0f*(Vx[i + (j+3)*nzz] - Vx[i + (j-2)*nzz]) +
                      8575.0f*(Vx[i + (j-1)*nzz] - Vx[i + (j+2)*nzz]) + 
                    128625.0f*(Vx[i + (j+1)*nzz] - Vx[i + j*nzz]))/(107520.0f*dx);

        float dVz_dz = (75.0f*(Vz[(i-3) + j*nzz] - Vz[(i+4) + j*nzz]) +   
                      1029.0f*(Vz[(i+3) + j*nzz] - Vz[(i-2) + j*nzz]) +
                      8575.0f*(Vz[(i-1) + j*nzz] - Vz[(i+2) + j*nzz]) +
                    128625.0f*(Vz[(i+1) + j*nzz] - Vz[i + j*nzz]))/(107520.0f*dz);     

        Txx[index] += dt*((L[index] + 2.0f*M[index])*dVx_dx + L[index]*dVz_dz);   
        Tzz[index] += dt*((L[index] + 2.0f*M[index])*dVz_dz + L[index]*dVx_dx);
    }

    if ((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3)) 
    {
        float dVz_dx = (75.0f*(Vz[i + (j-4)*nzz] - Vz[i + (j+3)*nzz]) +
                      1029.0f*(Vz[i + (j+2)*nzz] - Vz[i + (j-3)*nzz]) +
                      8575.0f*(Vz[i + (j-2)*nzz] - Vz[i + (j+1)*nzz]) +
                    128625.0f*(Vz[i + j*nzz]     - Vz[i + (j-1)*nzz]))/(107520.0f*dx);

        float dVx_dz = (75.0f*(Vx[(i-4) + j*nzz] - Vx[(i+3) + j*nzz]) +
                      1029.0f*(Vx[(i+2) + j*nzz] - Vx[(i-3) + j*nzz]) +
                      8575.0f*(Vx[(i-2) + j*nzz] - Vx[(i+1) + j*nzz]) +
                    128625.0f*(Vx[i + j*nzz]     - Vx[(i-1) + j*nzz]))/(107520.0f*dz);

        float Mxz = powf(0.25f*(1.0f/M[(i+1) + j*nzz]     + 1.0f/M[i + (j+1)*nzz] + 
                                1.0f/M[(i+1) + (j+1)*nzz] + 1.0f/M[i + j*nzz]), -1.0f); 

        Txz[index] += dt*Mxz*(dVx_dz + dVz_dx);            
    }          

    if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4))
    {    
        P[index] = 0.5f * (Txx[index] + Tzz[index]);    
    }
}

__global__ void compute_velocity(float * Vx, float * Vz, float * Txx, float * Tzz, float * Txz, float * B, float * d1D, float * d2D, int nb, int nxx, int nzz, float dx, float dz, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; 

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3)) 
    {
        float dTxx_dx = (75.0f*(Txx[i + (j-4)*nzz] - Txx[i + (j+3)*nzz]) +
                       1029.0f*(Txx[i + (j+2)*nzz] - Txx[i + (j-3)*nzz]) +
                       8575.0f*(Txx[i + (j-2)*nzz] - Txx[i + (j+1)*nzz]) +
                     128625.0f*(Txx[i + j*nzz]     - Txx[i + (j-1)*nzz]))/(107520.0f*dx);

        float dTxz_dz = (75.0f*(Txz[(i-3) + j*nzz] - Txz[(i+4) + j*nzz]) +
                       1029.0f*(Txz[(i+3) + j*nzz] - Txz[(i-2) + j*nzz]) + 
                       8575.0f*(Txz[(i-1) + j*nzz] - Txz[(i+2) + j*nzz]) +
                     128625.0f*(Txz[(i+1) + j*nzz] - Txz[i + j*nzz]))/(107520.0f*dz);

        float Bx = 0.5f*(B[i + (j+1)*nzz] + B[i + j*nzz]);

        Vx[index] += dt*Bx*(dTxx_dx + dTxz_dz);  
    }
    
    if ((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4)) 
    {
        float dTxz_dx = (75.0f*(Txz[i + (j-3)*nzz] - Txz[i + (j+4)*nzz]) +
                       1029.0f*(Txz[i + (j+3)*nzz] - Txz[i + (j-2)*nzz]) +
                       8575.0f*(Txz[i + (j-1)*nzz] - Txz[i + (j+2)*nzz]) +
                     128625.0f*(Txz[i + (j+1)*nzz] - Txz[i + j*nzz]))/(107520.0f*dx);

        float dTzz_dz = (75.0f*(Tzz[(i-4) + j*nzz] - Tzz[(i+3) + j*nzz]) + 
                       1029.0f*(Tzz[(i+2) + j*nzz] - Tzz[(i-3) + j*nzz]) +
                       8575.0f*(Tzz[(i-2) + j*nzz] - Tzz[(i+1) + j*nzz]) +
                     128625.0f*(Tzz[i + j*nzz]     - Tzz[(i-1) + j*nzz]))/(107520.0f*dz);

        float Bz = 0.5f*(B[(i+1) + j*nzz] + B[i + j*nzz]);

        Vz[index] += dt*Bz*(dTxz_dx + dTzz_dz); 
    }

    float damper = get_boundary_damper(d1D, d2D, i, j, nxx, nzz, nb);

    if (index < nxx*nzz)
    {    
        Vx[index] *= damper;    
        Vz[index] *= damper;    
        Txx[index] *= damper;    
        Tzz[index] *= damper;    
        Txz[index] *= damper;    
    }
}

__global__ void compute_seismogram(float * P, int * rIdx, int * rIdz, float * seismogram, int spread, int tId, int tlag, int nt, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if ((index < spread) && (tId >= tlag))
        seismogram[(tId - tlag) + index * nt] = P[rIdz[index] + rIdx[index]*nzz];
}

__device__ float get_boundary_damper(float * d1D, float * d2D, int i, int j, int nxx, int nzz, int nb)
{
    float damper;

    // global case
    if ((i >= nb) && (i < nzz - nb) && (j >= nb) && (j < nxx - nb))
    {
        damper = 1.0f;
    }

    // 1D damping
    else if ((i >= 0) && (i < nb) && (j >= nb) && (j < nxx - nb)) 
    {
        damper = d1D[i];
    }         
    else if ((i >= nzz - nb) && (i < nzz) && (j >= nb) && (j < nxx - nb)) 
    {
        damper = d1D[nb - (i - (nzz - nb)) - 1];
    }         
    else if ((i >= nb) && (i < nzz - nb) && (j >= 0) && (j < nb)) 
    {
        damper = d1D[j];
    }
    else if ((i >= nb) && (i < nzz - nb) && (j >= nxx - nb) && (j < nxx)) 
    {
        damper = d1D[nb - (j - (nxx - nb)) - 1];
    }

    // 2D damping 
    else if ((i >= 0) && (i < nb) && (j >= 0) && (j < nb))
    {
        damper = d2D[i + j*nb];
    }
    else if ((i >= nzz - nb) && (i < nzz) && (j >= 0) && (j < nb))
    {
        damper = d2D[nb - (i - (nzz - nb)) - 1 + j*nb];
    }
    else if((i >= 0) && (i < nb) && (j >= nxx - nb) && (j < nxx))
    {
        damper = d2D[i + (nb - (j - (nxx - nb)) - 1)*nb];
    }
    else if((i >= nzz - nb) && (i < nzz) && (j >= nxx - nb) && (j < nxx))
    {
        damper = d2D[nb - (i - (nzz - nb)) - 1 + (nb - (j - (nxx - nb)) - 1)*nb];
    }

    return damper;
}