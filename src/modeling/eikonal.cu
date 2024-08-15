# include "eikonal.cuh"

void Eikonal::set_eikonal_parameters()
{
    nSweeps = 4;
    meshDim = 2;

    dz2i = 1.0f / (dz*dz);    
    dx2i = 1.0f / (dx*dx);

    totalLevels = (nxx - 1) + (nzz - 1);
    
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
    for (int sweep = 0; sweep < 1; sweep++)
    { 



    }
}