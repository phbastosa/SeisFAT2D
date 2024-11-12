# include "adjoint_state.cuh"

void Adjoint_State::set_specifications()
{
    inversion_name = "adjoint_state_";
    inversion_method = "Adjoint-State First-Arrival Tomography";

    nSweeps = 4;
    meshDim = 2;
    nThreads = 32;      

    cell_area = modeling->dx*modeling->dz;

    total_levels = modeling->nxx + modeling->nzz - 1;

    m = new float[modeling->nPoints]();
    v = new float[modeling->nPoints]();

    source_grad = new float[modeling->matsize]();
    source_comp = new float[modeling->matsize]();
    
    adjoint_grad = new float[modeling->matsize]();
    adjoint_comp = new float[modeling->matsize]();

    gradient = new float[modeling->nPoints]();
    illumination = new float[modeling->nPoints]();     

    cudaMalloc((void**)&(d_T), modeling->matsize*sizeof(float));

    cudaMalloc((void**)&(d_source_grad), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_source_comp), modeling->matsize*sizeof(float));

    cudaMalloc((void**)&(d_adjoint_grad), modeling->matsize*sizeof(float));
    cudaMalloc((void**)&(d_adjoint_comp), modeling->matsize*sizeof(float));
}

void Adjoint_State::apply_inversion_technique()
{
    initialization();

    cudaMemcpy(d_T, modeling->T, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(d_source_grad, source_grad, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_source_comp, source_comp, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);    

    cudaMemcpy(d_adjoint_grad, adjoint_grad, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_adjoint_comp, adjoint_comp, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);    

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

            inner_sweep<<<nBlocks, nThreads>>>(d_T, d_adjoint_grad, d_adjoint_comp, d_source_grad, d_source_comp, x_offset, z_offset, xd, zd, modeling->nxx, modeling->nzz, modeling->dx, modeling->dz);

            cudaDeviceSynchronize();    
        }
    }

    cudaMemcpy(adjoint_grad, d_adjoint_grad, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(adjoint_comp, d_adjoint_comp, modeling->matsize*sizeof(float), cudaMemcpyDeviceToHost);

    # pragma omp parallel for reduction(+:gradient[:modeling->nPoints])
    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz);  

        int indp = i + j*modeling->nz; 
        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

        gradient[indp] += (adjoint_grad[indb] / (adjoint_comp[indb] + 1e-6f))*cell_area / modeling->geometry->nrel;
    }
}

void Adjoint_State::initialization()
{   
    # pragma omp parallel for
    for (int index = 0; index < modeling->matsize; index++) 
    {
        source_grad[index] = 0.0f;    
        source_comp[index] = 0.0f;    
        
        adjoint_grad[index] = 1e6f;
        adjoint_comp[index] = 1e6f;

        int i = (int) (index % modeling->nzz);    
        int j = (int) (index / modeling->nzz);  

        if ((i == 0) || (i == modeling->nzz - 1) || 
            (j == 0) || (j == modeling->nxx - 1)) 
        {    
            adjoint_grad[index] = 0.0f;        
            adjoint_comp[index] = 0.0f;        
        }
    }

    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    int sIdx = (int)(modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dx) + modeling->nb;
    int sIdz = (int)(modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dz) + modeling->nb;

    float So = modeling->S[sIdz + sIdx*modeling->nzz];

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

                float X = sqrtf(powf((sIdx - xi)*modeling->dx, 2.0f) + powf((sIdz - zi)*modeling->dz, 2.0f));

                int index = zi + xi*modeling->nzz;
                
                source_grad[index] += (dobs[spread + skipped] - modeling->T[index]) / cell_area;
                source_comp[index] += 1.0f / (X*X*So);
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

        perturbation[index] = max_slowness_variation*m_hat / (sqrtf(v_hat) + epsilon);
    }

    memset(gradient, 0.0f, modeling->nPoints);
}

__global__ void inner_sweep(float * T, float * adjoint_grad, float * adjoint_comp, float * source_grad, float * source_comp, int x_offset, int z_offset, int xd, int zd, int nxx, int nzz, float dx, float dz)
{
    int element = blockIdx.x*blockDim.x + threadIdx.x;

    int i = z_offset + zd*element;
    int j = x_offset + xd*element;    

    if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1))
    {
        float a1 = -1.0f*(T[i + j*nzz] - T[i + (j-1)*nzz]) / dx;
        float ap1 = 0.5f*(a1 + fabsf(a1));
        float am1 = 0.5f*(a1 - fabsf(a1));

        float a2 = -1.0f*(T[i + (j+1)*nzz] - T[i + j*nzz]) / dx;
        float ap2 = 0.5f*(a2 + fabsf(a2));
        float am2 = 0.5f*(a2 - fabsf(a2));

        float c1 = -1.0f*(T[i + j*nzz] - T[(i-1) + j*nzz]) / dz;
        float cp1 = 0.5f*(c1 + fabsf(c1));
        float cm1 = 0.5f*(c1 - fabsf(c1));

        float c2 = -1.0f*(T[(i+1) + j*nzz] - T[i + j*nzz]) / dz;
        float cp2 = 0.5f*(c2 + fabsf(c2));
        float cm2 = 0.5f*(c2 - fabsf(c2));

        float d = (ap2 - am1) / dx + (cp2 - cm1) / dz;

        if (fabsf(d) < 1e-6f)
        {
            adjoint_grad[i + j*nzz] = 0.0f;    
            adjoint_comp[i + j*nzz] = 0.0f;    
        }
        else
        {
            float eg = (ap1*adjoint_grad[i + (j-1)*nzz] - am2*adjoint_grad[i + (j+1)*nzz]) / dx +
                       (cp1*adjoint_grad[(i-1) + j*nzz] - cm2*adjoint_grad[(i+1) + j*nzz]) / dz;

            float ec = (ap1*adjoint_comp[i + (j-1)*nzz] - am2*adjoint_comp[i + (j+1)*nzz]) / dx +
                       (cp1*adjoint_comp[(i-1) + j*nzz] - cm2*adjoint_comp[(i+1) + j*nzz]) / dz;

            float fg = (eg + source_grad[i + j*nzz]) / d;
            float fc = (ec + source_comp[i + j*nzz]) / d;

            if (adjoint_grad[i + j*nzz] > fg) adjoint_grad[i + j*nzz] = fg;
            if (adjoint_comp[i + j*nzz] > fc) adjoint_comp[i + j*nzz] = fc;
        }
    }
}

