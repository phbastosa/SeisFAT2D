# include "modeling.cuh"

void Modeling::set_parameters()
{
    nx = std::stoi(catch_parameter("x_samples", parameters));        
    nz = std::stoi(catch_parameter("z_samples", parameters));        

    dx = std::stof(catch_parameter("x_spacing", parameters));
    dz = std::stof(catch_parameter("z_spacing", parameters));

    data_folder = catch_parameter("modeling_output_folder", parameters);

    nPoints = nx*nz;

    geometry = new Geometry();
    geometry->parameters = parameters;
    geometry->set_parameters();

    max_spread = 0;
    for (int index = 0; index < geometry->nrel; index++)
    {   
        if (max_spread < geometry->spread[index])
            max_spread = geometry->spread[index]; 
    }

    seismogram = new float[max_spread]();

    nb = 3;
    
    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    matsize = nxx*nzz;

    nThreads = 256;
    nBlocks = (int)((matsize + nThreads - 1) / nThreads);

    set_properties();    
    set_conditions();    
    set_eikonal();
}

void Modeling::set_properties()
{
    float * vp = new float[nPoints]();

    std::string vp_file = catch_parameter("vp_model_file", parameters);

    import_binary_float(vp_file, vp, nPoints);

    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
        vp[index] = 1.0f / vp[index];

    S = new float[matsize]();

    expand_boundary(vp, S);

    delete[] vp;
}

void Modeling::set_eikonal()
{
    dz2i = 1.0f / (dz * dz);
    dx2i = 1.0f / (dx * dx);

    total_levels = (nxx - 1) + (nzz - 1);

    std::vector<std::vector<int>> sgnv = {{1, 1}, {0, 1},  {1, 0}, {0, 0}};
    std::vector<std::vector<int>> sgnt = {{1, 1}, {-1, 1}, {1, -1}, {-1, -1}};

    int * h_sgnv = new int [NSWEEPS*MESHDIM]();
    int * h_sgnt = new int [NSWEEPS*MESHDIM](); 

    for (int index = 0; index < NSWEEPS*MESHDIM; index++)
    {
        int j = index / NSWEEPS;
    	int i = index % NSWEEPS;				

	    h_sgnv[i + j*NSWEEPS] = sgnv[i][j];
	    h_sgnt[i + j*NSWEEPS] = sgnt[i][j];    
    }

    T = new float[matsize]();

    cudaMalloc((void**)&(d_T), matsize*sizeof(float));
    cudaMalloc((void**)&(d_S), matsize*sizeof(float));

    cudaMalloc((void**)&(d_sgnv), NSWEEPS*MESHDIM*sizeof(int));
    cudaMalloc((void**)&(d_sgnt), NSWEEPS*MESHDIM*sizeof(int));
    
    copy_slowness_to_device();

    cudaMemcpy(d_sgnv, h_sgnv, NSWEEPS*MESHDIM*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sgnt, h_sgnt, NSWEEPS*MESHDIM*sizeof(int), cudaMemcpyHostToDevice);

    std::vector<std::vector<int>>().swap(sgnv);
    std::vector<std::vector<int>>().swap(sgnt);

    delete[] h_sgnt;
    delete[] h_sgnv;
}

void Modeling::set_shot_point()
{
    sx = geometry->xsrc[geometry->sInd[srcId]]; 
    sz = geometry->zsrc[geometry->sInd[srcId]]; 
}

void Modeling::initialization()
{
    sIdx = (int)((sx + 0.5f*dx) / dx) + nb;
    sIdz = (int)((sz + 0.5f*dz) / dz) + nb;

    time_set<<<nBlocks,nThreads>>>(d_T, matsize);

    dim3 grid(1,1,1);
    dim3 block(MESHDIM+1,MESHDIM+1,1);

    time_init<<<grid,block>>>(d_T,d_S,sx,sz,dx,dz,sIdx,sIdz,nzz,nb);
}

void Modeling::eikonal_solver()
{
    int min_level = std::min(nxx, nzz);
    int max_level = std::max(nxx, nzz);

    int z_offset, x_offset, n_elements;

    for (int sweep = 0; sweep < NSWEEPS; sweep++)
    { 
        int zd = (sweep == 2 || sweep == 3) ? -1 : 1; 
        int xd = (sweep == 0 || sweep == 2) ? -1 : 1;

        int sgni = sweep + 0*NSWEEPS;
        int sgnj = sweep + 1*NSWEEPS;

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

            int nblk = (int)((n_elements + nThreads - 1) / nThreads);

            inner_sweep<<<nblk, nThreads>>>(d_T, d_S, d_sgnv, d_sgnt, sgni, sgnj, x_offset, z_offset, xd, zd, nxx, nzz, dx, dz, dx2i, dz2i); 
        }
    }
}

void Modeling::compute_seismogram()
{
    int spread = 0;

    float P[4][4];

    cudaMemcpy(T, d_T, matsize*sizeof(float), cudaMemcpyDeviceToHost);    

    for (recId = geometry->iRec[srcId]; recId < geometry->fRec[srcId]; recId++)
    {
        float x = geometry->xrec[recId];
        float z = geometry->zrec[recId];

        float x0 = floorf((x + 0.5f*dx) / dx) * dx;
        float z0 = floorf((z + 0.5f*dz) / dz) * dz;

        float x1 = floorf((x + 0.5f*dx) / dx) * dx + dx;
        float z1 = floorf((z + 0.5f*dz) / dz) * dz + dz;

        float xd = (x - x0) / (x1 - x0);
        float zd = (z - z0) / (z1 - z0);

        int i = (int)((z + 0.5f*dz) / dz) + nb; 
        int j = (int)((x + 0.5f*dx) / dx) + nb;   

        for (int pIdx = 0; pIdx < 4; pIdx++)
        {
            for (int pIdz = 0; pIdz < 4; pIdz++)
            {    
                P[pIdx][pIdz] = T[(i + pIdz - 1) + (j + pIdx - 1)*nzz];
            }
        }   

        seismogram[spread++] = cubic2d(P, xd, zd);
    }
}

float Modeling::cubic1d(float P[4], float dx)
{
    return P[1] + 0.5f*dx*(P[2] - P[0] + dx*(2.0f*P[0] - 5.0f*P[1] + 4.0f*P[2] - P[3] + dx*(3.0f*(P[1] - P[2]) + P[3] - P[0])));
}

float Modeling::cubic2d(float P[4][4], float dx, float dy)
{    
    float p[4];
    p[0] = cubic1d(P[0], dy);
    p[1] = cubic1d(P[1], dy);
    p[2] = cubic1d(P[2], dy);
    p[3] = cubic1d(P[3], dy);    
    return cubic1d(p, dx);
}

void Modeling::export_seismogram()
{    
    compute_seismogram();

    std::string data_file = data_folder + modeling_type + "_nStations" + std::to_string(geometry->spread[srcId]) + "_shot_" + std::to_string(geometry->sInd[srcId]+1) + ".bin";
    export_binary_float(data_file, seismogram, geometry->spread[srcId]);    
}

void Modeling::expand_boundary(float * input, float * output)
{
    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
    {
        int i = (int) (index % nz);  
        int j = (int) (index / nz);    

        output[(i + nb) + (j + nb)*nzz] = input[i + j*nz];     
    }

    for (int i = 0; i < nb; i++)
    {
        for (int j = nb; j < nxx - nb; j++)
        {
            output[i + j*nzz] = output[nb + j*nzz];
            output[(nzz - i - 1) + j*nzz] = output[(nzz - nb - 1) + j*nzz];
        }
    }

    for (int i = 0; i < nzz; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            output[i + j*nzz] = output[i + nb*nzz];
            output[i + (nxx - j - 1)*nzz] = output[i + (nxx - nb - 1)*nzz];
        }
    }
}

void Modeling::reduce_boundary(float * input, float * output)
{
    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
    {
        int i = (int) (index % nz);  
        int j = (int) (index / nz);    

        output[i + j*nz] = input[(i + nb) + (j + nb)*nzz];
    }
}

void Modeling::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------------------\n";
    std::cout << "                                 \033[34mSeisFAT2D\033[0;0m\n";
    std::cout << "-------------------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (nz - 1)*dz << ", x = " << (nx - 1) * dx <<") m\n\n";

    std::cout << "Running shot " << srcId + 1 << " of " << geometry->nrel << " in total\n\n";

    std::cout << "Current shot position: (z = " << geometry->zsrc[geometry->sInd[srcId]] << 
                                       ", x = " << geometry->xsrc[geometry->sInd[srcId]] << ") m\n\n";

    std::cout << "Modeling type: " << modeling_name << "\n";
}

void Modeling::compression(float * input, uintc * output, int matsize, float &max_value, float &min_value)
{
    max_value =-1e20f;
    min_value = 1e20f;
    
    # pragma omp parallel for
    for (int index = 0; index < matsize; index++)
    {
        min_value = std::min(input[index], min_value);
        max_value = std::max(input[index], max_value);        
    }

    # pragma omp parallel for
    for (int index = 0; index < matsize; index++)
        output[index] = static_cast<uintc>(1.0f + (COMPRESS - 1)*(input[index] - min_value) / (max_value - min_value));
}

void Modeling::copy_slowness_to_device()
{
    cudaMemcpy(d_S, S, matsize * sizeof(float), cudaMemcpyHostToDevice);
}

__global__ void time_set(float * T, int matsize)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    if (index < matsize) T[index] = 1e6f;
}

__global__ void time_init(float * T, float * S, float sx, float sz, float dx, 
                          float dz, int sIdx, int sIdz, int nzz, int nb)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    int zi = sIdz + (i - 1);
    int xi = sIdx + (j - 1);

    int index = zi + xi*nzz;

    T[index] = S[index] * sqrtf(powf((xi - nb)*dx - sx, 2.0f) + 
                                powf((zi - nb)*dz - sz, 2.0f));
}

__global__ void inner_sweep(float * T, float * S, int * sgnv, int * sgnt, int sgni, int sgnj, 
                            int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, 
                            float dx, float dz, float dx2i, float dz2i)
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

        float t1d1 = tv + dz*min(S[i1 + max(j - 1, 1)*nzz], S[i1 + min(j, nxx - 1)*nzz]); 
        float t1d2 = te + dx*min(S[max(i - 1, 1) + j1*nzz], S[min(i, nzz - 1) + j1*nzz]); 

        float t1D = min(t1d1, t1d2);

        t1 = t2 = t3 = 1e6f; 

        Sref = S[i1 + j1*nzz];

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
