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

    m_hat = new float[modeling->nPoints]();
    v_hat = new float[modeling->nPoints]();

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

            nBlocks = (int)(n_elements / nThreads) + 1;

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

    int sId = modeling->geometry->sInd[modeling->srcId];

    int skipped = modeling->srcId * modeling->geometry->spread[sId];

    int sIdx = (int)(modeling->geometry->xsrc[sId] / modeling->dx);
    int sIdz = (int)(modeling->geometry->zsrc[sId] / modeling->dz);

    float So = modeling->S[sIdz + sIdx*modeling->nzz];

    int spreadId = 0;

    for (modeling->recId = modeling->geometry->iRec[sId]; modeling->recId < modeling->geometry->fRec[sId]; modeling->recId++)
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

                source_grad[zi + xi*modeling->nzz] += (dobs[spreadId + skipped] - modeling->T[zi + xi*modeling->nzz]) / cell_area;
                source_comp[zi + xi*modeling->nzz] += 1.0f / (X*X*So);
            }
        }

        ++spreadId;
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

        m_hat[index] = m[index] / (1.0f - powf(beta1, iteration));
        
        v_hat[index] = v[index] / (1.0f - powf(beta2, iteration));

        perturbation[index] = max_slowness_variation*m_hat[index] / (sqrtf(v_hat[index]) + epsilon);
    }

    memset(gradient, 0.0f, modeling->nPoints);
}

void Adjoint_State::gradient_preconditioning()
{
    float sigx = 0.005 * M_PI / modeling->dx;
    float sigz = 0.005 * M_PI / modeling->dz;

    float * kx = new float[modeling->nx]();
    float * kz = new float[modeling->nz]();

    for (int i = 0; i < modeling->nz; i++) 
    {
        float direction = (i <= modeling->nz / 2) ? (float) i : (float)(modeling->nz - i);

        kz[i] = 2.0f*direction*M_PI/(modeling->nz*modeling->dz);
    }

    for (int j = 0; j < modeling->nx; j++) 
    {
        float direction = (j <= modeling->nx / 2) ? (float) j : (float)(modeling->nx - j);

        kx[j] = 2.0f*direction*M_PI/(modeling->nx*modeling->dx);
    }

    fftw_complex * input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*modeling->nPoints);
    fftw_complex * output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*modeling->nPoints);

    for (int index = 0; index < modeling->nPoints; index++) 
        input[index][0] = static_cast<double>(gradient[index]); 

    fftw_plan forward_plan = fftw_plan_dft_2d(modeling->nz, modeling->nx, input, output, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan inverse_plan = fftw_plan_dft_2d(modeling->nz, modeling->nx, output, input, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(forward_plan);

    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz);  

        double gaussian_weight = 1.0 - exp(-(0.5*pow(kx[j] / sigx, 2.0) + 0.5*pow(kz[i] / sigz, 2.0)));

        output[i + j*modeling->nz][0] *= gaussian_weight; 
        output[i + j*modeling->nz][1] *= gaussian_weight; 
    }

    fftw_execute(inverse_plan);

    for (int index = 0; index < modeling->nPoints; index++) 
    {
        float preconditioning = static_cast<float>(input[index][0]) / modeling->nPoints / modeling->nPoints;        

        gradient[index] *= fabsf(preconditioning);        
    }    

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(inverse_plan);
    
    fftw_free(input);
    fftw_free(output);
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

