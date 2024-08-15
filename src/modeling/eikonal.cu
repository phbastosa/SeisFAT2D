# include "eikonal.cuh"

void Eikonal::set_eikonal_parameters()
{
    nSweeps = 4;
    meshDim = 2;

    dz2i = 1.0f / (dz*dz);    
    dx2i = 1.0f / (dx*dx);

    total_levels = nxx + nzz - 1;
    
    cudaMalloc((void**)&(T), matsize*sizeof(float));
    cudaMalloc((void**)&(S), matsize*sizeof(float));

    std::vector<std::vector<int>> sgnv = {{1,1}, {0,1},  {1,0}, {0,0}};
    std::vector<std::vector<int>> sgnt = {{1,1}, {-1,1}, {1,-1}, {-1,-1}};

    int * h_sgnv = new int [nSweeps*meshDim]();
    int * h_sgnt = new int [nSweeps*meshDim](); 

    for (int index = 0; index < nSweeps * meshDim; index++)
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
    sidx = (int)(geometry->xsrc[shot_index] / dx) + nb;
    sidz = (int)(geometry->zsrc[shot_index] / dz) + nb;

    for (int index = 0; index < matsize; index++) 
        eikonalT[index] = 1e6f;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int xi = sidx - (j - 1);
            int zi = sidz - (i - 1);

            eikonalT[zi + xi*nzz] = slowness[zi + xi*nzz] * sqrtf(powf(xi*dx - geometry->xsrc[shot_index], 2.0f) + 
                                                                  powf(zi*dz - geometry->zsrc[shot_index], 2.0f));
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

    for (int sweep = 0; sweep < nSweeps; sweep++)
    { 
        int zd = (sweep == 2 || sweep == 3) ? -1 : 1; 
        int xd = (sweep == 0 || sweep == 2) ? -1 : 1;

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
            else 
                min_level;

            std::cout << n_elements << "\n";



        }
    }
}


__global__ void fast_sweeping_method(int z_offset, int zd, int x_offset, int xd, int nxx, int nzz)
{
    int element = (threadIdx.x + blockIdx.x * blockDim.x);

    int i = z_offset + zd*element;
    int j = x_offset + xd*element;

    if ((i > 0) && (i < nzz - 1) && (j > 0) && (j < nxx - 1))
    {


    }


}








    // i1 = i - sgnvz
    // j1 = j - sgnvx

    // tv = T[i - sgntz, j]
    // te = T[i, j - sgntx]
    // tev = T[i - sgntz, j - sgntx]

    // Sref = min(S[i1, max(j - 1, 1)], S[i1, min(j, nxx - 1)])
    // t1d1 = tv + dz * Sref 

    // Sref = min(S[max(i-1, 1), j1], S[min(i, nzz - 1), j1])
    // t1d2 = te + dx * Sref 

    // t1D = min(t1d1, t1d2)

    // t1 = t2 = t3 = 1e6

    // if (tv <= te + dx * Sref) and (te <= tv + dz * Sref) and (te - tev >= 0.0) and (tv - tev >= 0.0):
        
    //     ta = tev + te - tv
    //     tb = tev - te + tv
    //     t1 = ((tb * dz2i + ta * dx2i) + np.sqrt(4.0 * Sref**2 * (dz2i + dx2i) - dz2i * dx2i * (ta - tb)**2)) / (dz2i + dx2i)
    
    // elif (te - tev <= Sref*dz**2 / np.sqrt(dx**2 + dz**2)) and (te - tev >= 0.0):
        
    //     t2 = te + dx*np.sqrt(Sref**2 - ((te - tev) / dz)**2)

    // elif (tv - tev <= Sref*dx**2 / np.sqrt(dx**2 + dz**2)) and (tv - tev >= 0.0):
        
    //     t3 = tv + dz * np.sqrt(Sref**2 - ((tv - tev) / dx)**2)

    // t2D = min(t1, t2, t3)

    // T[i,j] = min(T[i,j], t1D, t2D)
