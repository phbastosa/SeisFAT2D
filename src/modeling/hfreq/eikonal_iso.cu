# include "eikonal_iso.cuh"

void Eikonal_Iso::set_properties()
{
    Vp = new float[nPoints]();

    std::string model_file = catch_parameter("vp_model_file", parameters);

    import_binary_float(model_file, Vp, nPoints);

    for (int index = 0; index < nPoints; index++)
        Vp[index] = 1.0f / Vp[index];

    S = new float[matsize]();

    expand_boundary(Vp, S);
}

void Eikonal_Iso::set_conditions()
{
    modeling_type = "eikonal_iso_";
    modeling_name = "Eikonal isotropic time propagation";

    nSweeps = 4;
    meshDim = 2;

    nThreads = 32;      

    total_levels = (nxx - 1) + (nzz - 1);

    std::vector<std::vector<int>> sgnv = {{1, 1}, {0, 1},  {1, 0}, {0, 0}};
    std::vector<std::vector<int>> sgnt = {{1, 1}, {-1, 1}, {1, -1}, {-1, -1}};

    int * h_sgnv = new int [nSweeps*meshDim]();
    int * h_sgnt = new int [nSweeps*meshDim](); 

    for (int index = 0; index < nSweeps*meshDim; index++)
    {
        int j = index / nSweeps;
    	int i = index % nSweeps;				

	    h_sgnv[i + j*nSweeps] = sgnv[i][j];
	    h_sgnt[i + j*nSweeps] = sgnt[i][j];    
    }

    T = new float[matsize]();

    cudaMalloc((void**)&(d_T), matsize*sizeof(float));
    cudaMalloc((void**)&(d_S), matsize*sizeof(float));

    cudaMalloc((void**)&(d_sgnv), nSweeps*meshDim*sizeof(int));
    cudaMalloc((void**)&(d_sgnt), nSweeps*meshDim*sizeof(int));

    cudaMemcpy(d_sgnv, h_sgnv, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sgnt, h_sgnt, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);

    delete[] h_sgnt;
    delete[] h_sgnv;

    std::vector<std::vector<int>>().swap(sgnv);
    std::vector<std::vector<int>>().swap(sgnt);
}

void Eikonal_Iso::forward_solver()
{
    float dz2i = 1.0f / (dz * dz);
    float dx2i = 1.0f / (dx * dx);

    int min_level = std::min(nxx, nzz);
    int max_level = std::max(nxx, nzz);

    int z_offset, x_offset, n_elements;

    cudaMemcpy(d_T, T, matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_S, S, matsize*sizeof(float), cudaMemcpyHostToDevice);

    for (int sweep = 0; sweep < nSweeps; sweep++)
    { 
        int zd = (sweep == 2 || sweep == 3) ? -1 : 1; 
        int xd = (sweep == 0 || sweep == 2) ? -1 : 1;

        int sgni = sweep + 0*nSweeps;
        int sgnj = sweep + 1*nSweeps;

        for (int level = 0; level < total_levels; level++)
        {
            z_offset = (sweep == 0) ? ((level < nxx) ? 0 : level - nxx + 1) :
                    (sweep == 1) ? ((level < nzz) ? nzz - level - 1 : 0) :
                    (sweep == 2) ? ((level < nzz) ? level : nzz - 1) :
                                    ((level < nxx) ? nzz - 1 : nzz - 1 - (level - nxx + 1));

            x_offset = (sweep == 0) ? ((level < nxx) ? level : nxx - 1) :
                    (sweep == 1) ? ((level < nzz) ? 0 : level - nzz + 1) :
                    (sweep == 2) ? ((level < nzz) ? nxx - 1 : nxx - 1 - (level - nzz + 1)) :
                                    ((level < nxx) ? nxx - level - 1 : 0);

            n_elements = (level < min_level) ? level + 1 : 
                        (level >= max_level) ? total_levels - level : 
                        total_levels - min_level - max_level + level;

            nBlocks = (int)((n_elements + nThreads - 1) / nThreads);

            inner_sweep<<<nBlocks, nThreads>>>(d_T, d_S, d_sgnv, d_sgnt, sgni, sgnj, x_offset, z_offset, xd, zd, nxx, nzz, dx, dz, dx2i, dz2i);

            cudaDeviceSynchronize();    
        }
    }

    cudaMemcpy(T, d_T, matsize*sizeof(float), cudaMemcpyDeviceToHost);

    compute_seismogram();
}

__global__ void inner_sweep(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz, float dx2i, float dz2i)
{
    int element = blockIdx.x*blockDim.x + threadIdx.x;

    int i = z_offset + zd*element;
    int j = x_offset + xd*element;

    float Sref, t1, t2, t3;  

    if ((i > 0) && (i < nzz - 1) && (j > 0) && (j < nxx - 1))
    {
        int i1 = i - sgnv[sgni];
        int j1 = j - sgnv[sgnj];

        float tv = T[i - sgnt[sgni] + j*nzz];
        float te = T[i + (j - sgnt[sgnj])*nzz];
        float tev = T[(i - sgnt[sgni]) + (j - sgnt[sgnj])*nzz];

        Sref = min(S[i1 + max(j - 1, 1)*nzz], S[i1 + min(j, nxx - 1)*nzz]);
        
        float t1d1 = tv + dz*Sref; 

        Sref = min(S[max(i - 1, 1) + j1*nzz], S[min(i, nzz - 1) + j1*nzz]);

        float t1d2 = te + dx*Sref; 

        float t1D = min(t1d1, t1d2);

        t1 = t2 = t3 = 1e6f; 

        if ((tv <= te + dx*Sref) && (te <= tv + dz*Sref) && (te - tev >= 0.0f) && (tv - tev >= 0.0f))
        {
            float ta = tev + te - tv;
            float tb = tev - te + tv;

            t1 = ((tb*dz2i + ta*dx2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
        }
        else if ((te - tev <= Sref*dz*dz / sqrtf(dx*dx + dz*dz)) && (te - tev > 0.0f))
        {
            t2 = te + dx*sqrtf(Sref*Sref - ((te - tev) / dz)*((te - tev) / dz));
        }    
        else if ((tv - tev <= Sref*dx*dx / sqrt(dx*dx + dz*dz)) && (tv - tev > 0.0f))
        {
            t3 = tv + dz*sqrtf(Sref*Sref - ((tv - tev) / dx)*((tv - tev) / dx));
        }    

        float t2D = min(t1, min(t2, t3));

        T[i + j*nzz] = min(T[i + j*nzz], min(t1D, t2D));
    }
}