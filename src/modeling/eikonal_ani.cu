# include "eikonal_ani.cuh"

void Eikonal_ANI::set_conditions()
{
    modeling_type = "eikonal_ani";
    modeling_name = "Modeling type: Anisotropic eikonal solver";

    auto * Cij = new float[nPoints]();

    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    auto * C11 = new float[matsize]();
    auto * uC11 = new uintc[matsize]();
    import_binary_float(Cijkl_folder + "C11.bin", Cij, nPoints);
    expand_boundary(Cij, C11);
    compression(C11, uC11, matsize, maxC11, minC11);        
    cudaMalloc((void**)&(d_C11), matsize*sizeof(uintc));
    cudaMemcpy(d_C11, uC11, matsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C11;
    delete[] uC11;

    auto * C13 = new float[matsize]();
    auto * uC13 = new uintc[matsize]();
    import_binary_float(Cijkl_folder + "C13.bin", Cij, nPoints);
    expand_boundary(Cij, C13);
    compression(C13, uC13, matsize, maxC13, minC13);    
    cudaMalloc((void**)&(d_C13), matsize*sizeof(uintc));
    cudaMemcpy(d_C13, uC13, matsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C13;
    delete[] uC13;

    auto * C15 = new float[matsize]();
    auto * uC15 = new uintc[matsize]();
    import_binary_float(Cijkl_folder + "C15.bin", Cij, nPoints);
    expand_boundary(Cij, C15);
    compression(C15, uC15, matsize, maxC15, minC15);    
    cudaMalloc((void**)&(d_C15), matsize*sizeof(uintc));
    cudaMemcpy(d_C15, uC15, matsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C15;
    delete[] uC15;

    auto * C33 = new float[matsize]();
    auto * uC33 = new uintc[matsize]();
    import_binary_float(Cijkl_folder + "C33.bin", Cij, nPoints);
    expand_boundary(Cij, C33);
    compression(C33, uC33, matsize, maxC33, minC33);    
    cudaMalloc((void**)&(d_C33), matsize*sizeof(uintc));
    cudaMemcpy(d_C33, uC33, matsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C33;
    delete[] uC33;
    
    auto * C35 = new float[matsize]();
    auto * uC35 = new uintc[matsize]();
    import_binary_float(Cijkl_folder + "C35.bin", Cij, nPoints);
    expand_boundary(Cij, C35);
    compression(C35, uC35, matsize, maxC35, minC35);    
    cudaMalloc((void**)&(d_C35), matsize*sizeof(uintc));
    cudaMemcpy(d_C35, uC35, matsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C35;
    delete[] uC35;

    auto * C55 = new float[matsize]();
    auto * uC55 = new uintc[matsize]();
    import_binary_float(Cijkl_folder + "C55.bin", Cij, nPoints);
    expand_boundary(Cij, C55);
    compression(C55, uC55, matsize, maxC55, minC55);    
    cudaMalloc((void**)&(d_C55), matsize*sizeof(uintc));
    cudaMemcpy(d_C55, uC55, matsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C55;
    delete[] uC55;

    delete[] Cij;
}

void Eikonal_ANI::time_propagation()
{
    initialization();
    eikonal_solver();

    get_quasi_slowness<<<nBlocks,nThreads>>>(d_T,d_S,dx,dz,sIdx,sIdz,nxx,nzz,nb,d_C11,d_C13,d_C15,d_C33,d_C35,d_C55,minC11,
                                             maxC11,minC13,maxC13,minC15,maxC15,minC33,maxC33,minC35,maxC35,minC55,maxC55);
                                             
    initialization();
    eikonal_solver();

    cudaMemcpy(d_S, S, matsize * sizeof(float), cudaMemcpyHostToDevice);
}

__global__ void get_quasi_slowness(float * T, float * S, float dx, float dz, int sIdx, int sIdz, int nxx, int nzz, int nb,
                                   uintc * C11, uintc * C13, uintc * C15, uintc * C33, uintc * C35, uintc * C55, 
                                   float minC11, float maxC11, float minC13, float maxC13, float minC15, float maxC15, 
                                   float minC33, float maxC33, float minC35, float maxC35, float minC55, float maxC55)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;    
    int i = (int) (index - j*nzz - k*nxx*nzz);  

    const int n = 2;
    const int v = 3;

    float p[n];
    float C[v*v];
    float Gv[n];

    if ((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb))
    {
        if (!((i == sIdz) && (j == sIdx)))    
        {
            float dTz = 0.5f*(T[(i+1) + j*nzz] - T[(i-1) + j*nzz]) / dz;
            float dTx = 0.5f*(T[i + (j+1)*nzz] - T[i + (j-1)*nzz]) / dx;

            float norm = sqrtf(dTx*dTx + dTz*dTz);

            p[0] = dTx / norm;
            p[1] = dTz / norm;
            
            float c11 = (minC11 + (static_cast<float>(C11[index]) - 1.0f) * (maxC11 - minC11) / (COMPRESS - 1));
            float c13 = (minC13 + (static_cast<float>(C13[index]) - 1.0f) * (maxC13 - minC13) / (COMPRESS - 1));
            float c15 = (minC15 + (static_cast<float>(C15[index]) - 1.0f) * (maxC15 - minC15) / (COMPRESS - 1));

            float c33 = (minC33 + (static_cast<float>(C33[index]) - 1.0f) * (maxC33 - minC33) / (COMPRESS - 1));
            float c35 = (minC35 + (static_cast<float>(C35[index]) - 1.0f) * (maxC35 - minC35) / (COMPRESS - 1));

            float c55 = (minC55 + (static_cast<float>(C55[index]) - 1.0f) * (maxC55 - minC55) / (COMPRESS - 1));

            C[0+0*v] = c11; C[0+1*v] = c13; C[0+2*v] = c15;
            C[1+0*v] = c13; C[1+1*v] = c33; C[1+2*v] = c35;
            C[2+0*v] = c15; C[2+1*v] = c35; C[2+2*v] = c55;

            float Ro = c33*S[index]*S[index];    
            
            for (int indp = 0; indp < v*v; indp++)
                C[indp] = C[indp] / Ro / Ro;

            float Gxx = C[0+0*v]*p[0]*p[0] + C[2+2*v]*p[1]*p[1] + 2.0f*C[0+2*v]*p[0]*p[1];
            float Gzz = C[2+2*v]*p[0]*p[0] + C[1+1*v]*p[1]*p[1] + 2.0f*C[1+2*v]*p[0]*p[1];
            float Gxz = C[0+2*v]*p[0]*p[0] + C[1+2*v]*p[1]*p[1] + (C[0+1*v] + C[2+2*v])*p[0]*p[1]; 
            
            float coeff1 = Gxx + Gzz;
            float coeff2 = Gxx - Gzz;
            
            float det = std::sqrt((coeff2 * coeff2) / 4.0f + Gxz * Gxz);

            Gv[0] = coeff1 / 2.0 + det;
            Gv[1] = coeff1 / 2.0 - det;
            
            if (Gv[0] < Gv[1]) {float aux = Gv[0]; Gv[0] = Gv[1]; Gv[1] = aux;} 

            S[index] = 1.0f / sqrtf(Gv[0] * Ro);
        }
    }
}