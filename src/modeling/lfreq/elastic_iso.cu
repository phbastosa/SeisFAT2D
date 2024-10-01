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
        Vp[index] = 1500.0f;
        Vs[index] = 0.0f;
        Rho[index] = 1000.0f;

        M[index] = Rho[index]*Vs[index]*Vs[index];
        L[index] = Rho[index]*Vp[index]*Vp[index] - 2.0f*M[index];
        B[index] = 1.0f / Rho[index];
    }
            
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


    // snapshots
    // seismogram

    sIdx = (int)(geometry->xsrc[geometry->sInd[srcId]] / dx) + nb;
    sIdz = (int)(geometry->zsrc[geometry->sInd[srcId]] / dz) + nb;
}

void elastic_Iso::forward_solver()
{
    for (int tId = 0; tId < tlag; tId++)
    {
        compute_pressure<<<nBlocks, nThreads>>>(d_Vx, d_Vz, d_Txx, d_Tzz, d_Txz, d_P, d_M, d_L, wavelet, sIdx, sIdz, tId, nt, nxx, nzz, dx, dz, dt);
        cudaDeviceSynchronize();

        compute_velocity<<<nBlocks, nThreads>>>(d_Vx, d_Vz, d_Txx, d_Tzz, d_Txz, d_B, nxx, nzz, dx, dz, dt);
        cudaDeviceSynchronize();
    }

    std::cout << fmax << " " << dt << std::endl;
    std::cout << nxx << " " << nzz << std::endl;
    std::cout << sIdx << " " << sIdz << std::endl;
    std::cout << nBlocks << " " << nThreads << std::endl;

    cudaMemcpy(P, d_Vx, matsize*sizeof(float), cudaMemcpyDeviceToHost);
    export_binary_float("vx_shot_" + std::to_string(srcId+1) + ".bin", P, matsize);

    cudaMemcpy(P, d_Txx, matsize*sizeof(float), cudaMemcpyDeviceToHost);
    export_binary_float("txx_shot_" + std::to_string(srcId+1) + ".bin", P, matsize);

}

__global__ void compute_pressure(float * Vx, float * Vz, float * Txx, float * Tzz, float * Txz, float * P, float * M, float * L, float * wavelet, int sIdx, int sIdz, int tId, int nt, int nxx, int nzz, float dx, float dz, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index / nzz);
    int j = (int)(index % nzz);

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

__global__ void compute_velocity(float * Vx, float * Vz, float * Txx, float * Tzz, float * Txz, float * B, int nxx, int nzz, float dx, float dz, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; 

    int i = (int)(index / nzz);
    int j = (int)(index % nzz);

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
}