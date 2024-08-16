# include "eikonal.cuh"

void Eikonal::set_eikonal_parameters()
{
    nSweeps = 4;
    meshDim = 2;

    threadsPerBlock = 32;      

    total_levels = nxx + nzz - 1;
    
    cudaMalloc((void**)&(T), matsize*sizeof(float));
    cudaMalloc((void**)&(S), matsize*sizeof(float));

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

    cudaMalloc((void**)&(d_sgnv), nSweeps*meshDim*sizeof(int));
    cudaMalloc((void**)&(d_sgnt), nSweeps*meshDim*sizeof(int));

    cudaMemcpy(d_sgnv, h_sgnv, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sgnt, h_sgnt, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);

    delete[] h_sgnt;
    delete[] h_sgnv;

    std::vector<std::vector<int>>().swap(sgnv);
    std::vector<std::vector<int>>().swap(sgnt);
}

void Eikonal::initialization()
{
    sidx = (int)(geometry->xsrc[shotId] / dx) + nb;
    sidz = (int)(geometry->zsrc[shotId] / dz) + nb;

    for (int index = 0; index < matsize; index++) 
        eikonalT[index] = 1e6f;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int xi = sidx - (j - 1);
            int zi = sidz - (i - 1);

            eikonalT[zi + xi*nzz] = slowness[zi + xi*nzz] * sqrtf(powf(xi*dx - geometry->xsrc[shotId], 2.0f) + 
                                                                  powf(zi*dz - geometry->zsrc[shotId], 2.0f));
        }
    }

    cudaMemcpy(T, eikonalT, matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(S, slowness, matsize*sizeof(float), cudaMemcpyHostToDevice);
}

void Eikonal::forward_solver()
{
    int min_level = std::min(nx, nz);
    int max_level = std::max(nx, nz);

    int z_offset;
    int x_offset;
    int n_elements;

    float dz2i = 1.0f / (dz*dz);    
    float dx2i = 1.0f / (dx*dx);

    for (int sweep = 0; sweep < nSweeps; sweep++)
    { 
        int zd = (sweep == 2 || sweep == 3) ? -1 : 1; 
        int xd = (sweep == 0 || sweep == 2) ? -1 : 1;

        int sgni = sweep + 0*nSweeps;
        int sgnj = sweep + 1*nSweeps;

        for (int level = 0; level < total_levels; level++)
        {
            if (sweep == 0)
            {
                z_offset = (level < nxx) ? 0 : level - nxx + 1;
                x_offset = (level < nxx) ? level : nxx - 1;
            }
            else if (sweep == 1)
            {
                z_offset = (level < nzz) ? nzz - level - 1 : 0;
                x_offset = (level < nzz) ? 0 : level - nzz + 1;
            }
            else if (sweep == 2)
            {
                z_offset = (level < nzz) ? level : nzz - 1;
                x_offset = (level < nzz) ? nxx - 1 : nxx - 1 - (level - nzz + 1);
            }
            else if (sweep == 3)
            {
                z_offset = (level < nxx) ? nzz - 1 : nzz - 1 - (level - nxx + 1);
                x_offset = (level < nxx) ? nxx - level - 1 : 0;
            }

            if (level < min_level) 
                n_elements = level + 1;  
            else if (level >= max_level) 
                n_elements = total_levels - level;


            blocksPerGrid = (int)((n_elements - 1) / threadsPerBlock) + 1;

            fast_sweeping_method<<<blocksPerGrid, threadsPerBlock>>>(T, S, d_sgnv, d_sgnt, sgni, sgnj, x_offset, z_offset, xd, zd, nxx, nzz, dx, dz, dx2i, dz2i);

            cudaDeviceSynchronize();    
        }
    }

    cudaMemcpy(eikonalT, T, matsize*sizeof(float), cudaMemcpyDeviceToHost);
}

__global__ void fast_sweeping_method(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz, float dx2i, float dz2i)
{
    int element = (threadIdx.x + blockIdx.x*blockDim.x);

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
        else if ((te - tev <= Sref*dz*dz / sqrtf(dx*dx + dz*dz)) && (te - tev >= 0.0f))
        {
            t2 = te + dx*sqrtf(Sref*Sref - ((te - tev) / dz)*((te - tev) / dz));
        }    
        else if ((tv - tev <= Sref*dx*dx / sqrt(dx*dx + dz*dz)) && (tv - tev >= 0.0f))
        {
            t3 = tv + dz*sqrtf(Sref*Sref - ((tv - tev) / dx)*((tv - tev) / dx));
        }    

        float t2D = min(t1, min(t2, t3));

        T[i + j*nzz] = min(T[i + j*nzz], min(t1D, t2D));
    }
}
