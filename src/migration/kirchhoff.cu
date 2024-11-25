# include "kirchhoff.cuh"

void Kirchhoff::set_specifications()
{
    cudaMalloc((void**)&(d_Tr), modeling->nPoints*sizeof(float));
    cudaMalloc((void**)&(d_Ts), modeling->nPoints*sizeof(float));
    cudaMalloc((void**)&(d_image), modeling->nPoints*sizeof(float));
    cudaMalloc((void**)&(d_gather), modeling->nz*modeling->max_spread*sizeof(float));
    cudaMalloc((void**)&(d_seismic), modeling->nt*modeling->max_spread*sizeof(float));

    nThreads = 256;
    nBlocks = (int)((modeling->nPoints + nThreads - 1) / nThreads);
}

void Kirchhoff::run_cross_correlation()
{
    cudaMemset(d_image, 0.0f, modeling->nPoints*sizeof(float));

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        read_seismic_data();

        modeling->show_information();

        std::cout << "\nKirchhoff depth migration: computing image matrix\n";

        modeling->initialization();
        modeling->forward_solver();

        modeling->reduce_boundary(modeling->T, Ts);

        cudaMemcpy(d_Ts, Ts, modeling->nPoints*sizeof(float), cudaMemcpyHostToDevice);
        
        cudaMemset(d_gather, 0.0f, modeling->nz*modeling->geometry->spread[modeling->srcId]*sizeof(float));
        cudaMemcpy(d_seismic, seismic, modeling->nt*modeling->geometry->spread[modeling->srcId]*sizeof(float), cudaMemcpyHostToDevice);
        
        int spread = 0;

        float sx = modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]];

        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            import_binary_float(output_table_folder + "traveltimes_receiver_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);

            float rx = modeling->geometry->xrec[modeling->recId];

            float cmp = sx + 0.5f*(rx - sx);

            cudaMemcpy(d_Tr, Tr, modeling->nPoints*sizeof(float), cudaMemcpyHostToDevice);

            cross_correlation<<<nBlocks, nThreads>>>(d_Ts, d_Tr, d_image, d_gather, d_seismic, aperture, cmp, modeling->nPoints, spread, modeling->nz, modeling->nt, modeling->dt, modeling->dx, modeling->dz);

            ++spread;
        }

        cudaMemcpy(gather, d_gather, modeling->nz*modeling->geometry->spread[modeling->srcId]*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(output_image_folder + "gather_" + std::to_string(modeling->nz) + "x" + std::to_string(spread) + "_shot_" + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin", gather, modeling->nz*modeling->geometry->spread[modeling->srcId]);
    }

    cudaMemcpy(image, d_image, modeling->nPoints*sizeof(float), cudaMemcpyDeviceToHost);

    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        image[index] *= 1.0f / modeling->geometry->nrel;
}

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * gather, float * seismic, float aperture, float cmp, int nPoints, int spread, int nz, int nt, float dt, float dx, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < nPoints)
    {
        int i = (int)(index % nz);
        int j = (int)(index / nz);

        float sigx = tanf(aperture*PI/180.0f)*i*dz;        
        float value = expf(-0.5*powf((j*dx - cmp)/(sigx + 1e-6f), 2.0f));

        float T = Ts[index] + Tr[index]; 
    
        int tId = (int)(T / dt);

        if (tId < nt) 
        {
            image[index] += value * seismic[tId + spread*nt];

            if (j == (int)(cmp/dx))
                gather[i + spread*nz] = value * seismic[tId + spread*nt];
        }    
    }
}
