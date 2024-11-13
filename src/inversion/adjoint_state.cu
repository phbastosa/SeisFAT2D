# include "adjoint_state.cuh"

void Adjoint_State::set_specifications()
{
    inversion_name = "adjoint_state_";
    inversion_method = "Adjoint-State First-Arrival Tomography";

    nSweeps = 4;
    meshDim = 2;
    nThreads = 32;      

    cell_area = modeling->dx*modeling->dz;

    total_levels = (modeling->nxx - 1) + (modeling->nzz - 1);

    m = new float[modeling->nPoints]();
    v = new float[modeling->nPoints]();

    source = new float[modeling->matsize]();
    adjoint = new float[modeling->matsize]();
    gradient = new float[modeling->nPoints]();

    cudaMalloc((void**)&(d_T), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_source), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_adjoint), modeling->matsize*sizeof(float));
}

void Adjoint_State::apply_inversion_technique()
{
    initialization();

    cudaMemcpy(d_T, modeling->T, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_source, source, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_adjoint, adjoint, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    int min_level = std::min(modeling->nxx, modeling->nzz);
    int max_level = std::max(modeling->nxx, modeling->nzz);

    int z_offset, x_offset, n_elements;

    for (int sweep = 0; sweep < nSweeps; sweep++)
    { 
        int zd = (sweep == 2 || sweep == 3) ? -1 : 1; 
        int xd = (sweep == 0 || sweep == 2) ? -1 : 1;

        for (int level = 0; level < total_levels; level++)
        {
            z_offset = (sweep == 0) ? ((level < modeling->nxx) ? 0 : level - modeling->nxx + 1) :
                       (sweep == 1) ? ((level < modeling->nzz) ? modeling->nzz - level - 1 : 0) :
                       (sweep == 2) ? ((level < modeling->nzz) ? level : modeling->nzz - 1) :
                                      ((level < modeling->nxx) ? modeling->nzz - 1 : modeling->nzz - 1 - (level - modeling->nxx + 1));

            x_offset = (sweep == 0) ? ((level < modeling->nxx) ? level : modeling->nxx - 1) :
                       (sweep == 1) ? ((level < modeling->nzz) ? 0 : level - modeling->nzz + 1) :
                       (sweep == 2) ? ((level < modeling->nzz) ? modeling->nxx - 1 : modeling->nxx - 1 - (level - modeling->nzz + 1)) :
                                      ((level < modeling->nxx) ? modeling->nxx - level - 1 : 0);

            n_elements = (level < min_level) ? level + 1 : 
                         (level >= max_level) ? total_levels - level : 
                         total_levels - min_level - max_level + level;

            nBlocks = (int)((n_elements + nThreads - 1) / nThreads);

            inner_sweep<<<nBlocks, nThreads>>>(d_T, d_adjoint, d_source, x_offset, z_offset, xd, zd, modeling->nxx, modeling->nzz, modeling->dx, modeling->dz);

            cudaDeviceSynchronize();    
        }
    }

    cudaMemcpy(adjoint, d_adjoint, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz);  

        int indp = i + j*modeling->nz; 
        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

        gradient[indp] += cell_area*adjoint[indb] / modeling->geometry->nrel;
    }
}

void Adjoint_State::initialization()
{   
    # pragma omp parallel for
    for (int index = 0; index < modeling->matsize; index++) 
    {
        source[index] = 0.0f;    
        adjoint[index] = 1e6f;

        int i = (int) (index % modeling->nzz);    
        int j = (int) (index / modeling->nzz);  

        if ((i == 0) || (i == modeling->nzz - 1) || 
            (j == 0) || (j == modeling->nxx - 1)) 
            adjoint[index] = 0.0f;        
    }

    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    int sIdx = (int)(modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dx) + modeling->nb;
    int sIdz = (int)(modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dz) + modeling->nb;

    int spread = 0;

    for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    {
        int rIdx = (int)(modeling->geometry->xrec[modeling->recId] / modeling->dx) + modeling->nb;
        int rIdz = (int)(modeling->geometry->zrec[modeling->recId] / modeling->dz) + modeling->nb;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int xi = rIdx + (j - 1);
                int zi = rIdz + (i - 1);

                int index = zi + xi*modeling->nzz;
                
                source[index] += (dobs[spread + skipped] - modeling->T[index]) / cell_area;
            }
        }

        ++spread;
    }   
} 

void Adjoint_State::optimization()
{   
    float gdot = 0.0f;
    #pragma omp parallel for reduction(+:gdot)
    for (int index = 0; index < modeling->nPoints; index++)
        gdot += gradient[index]*gradient[index];
    
    float beta1 = 0.5f;
    float beta2 = 0.9f;

    float epsilon = 1e-8f;

    for (int index = 0; index < modeling->nPoints; index++)
    {
        gradient[index] *= 1.0f / gdot;

        m[index] = beta1*m[index] + (1.0f - beta1)*gradient[index];
        
        v[index] = beta2*v[index] + (1.0f - beta2)*gradient[index]*gradient[index];

        float m_hat = m[index] / (1.0f - powf(beta1, iteration));
        
        float v_hat = v[index] / (1.0f - powf(beta2, iteration));

        perturbation[index] = powf(0.5, iteration)*max_slowness_variation*m_hat/(sqrtf(v_hat) + epsilon);
    }

    memset(gradient, 0.0f, modeling->nPoints);
}

__global__ void inner_sweep(float * T, float * adjoint, float * source, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz)
{
    int element = blockIdx.x*blockDim.x + threadIdx.x;

    int i = z_offset + zd*element;
    int j = x_offset + xd*element;    

    if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1))
    {
        float a1 = -1.0f*(T[i + j*nzz] - T[i + (j-1)*nzz]) / dx;
        float ap1 = (a1 + fabsf(a1))/2.0f;
        float am1 = (a1 - fabsf(a1))/2.0f;

        float a2 = -1.0f*(T[i + (j+1)*nzz] - T[i + j*nzz]) / dx;
        float ap2 = (a2 + fabsf(a2))/2.0f;
        float am2 = (a2 - fabsf(a2))/2.0f;

        float c1 = -1.0f*(T[i + j*nzz] - T[(i-1) + j*nzz]) / dz;
        float cp1 = (c1 + fabsf(c1))/2.0f;
        float cm1 = (c1 - fabsf(c1))/2.0f;

        float c2 = -1.0f*(T[(i+1) + j*nzz] - T[i + j*nzz]) / dz;
        float cp2 = (c2 + fabsf(c2))/2.0f;
        float cm2 = (c2 - fabsf(c2))/2.0f;

        float d = (ap2 - am1) / dx + (cp2 - cm1) / dz;

        if (fabsf(d) < 1e-6f)
        {
            adjoint[i + j*nzz] = 0.0f;    
        }
        else
        {
            float e = (am1*adjoint[i + (j-1)*nzz] - ap2*adjoint[i + (j+1)*nzz]) / dx +
                      (cm1*adjoint[(i-1) + j*nzz] - cp2*adjoint[(i+1) + j*nzz]) / dz;


            float f = (e + source[i + j*nzz]) / d;

            if (adjoint[i + j*nzz] > f) adjoint[i + j*nzz] = f;
        }
    }
}

