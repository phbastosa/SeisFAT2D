# include "adjoint_state.hpp"

void Adjoint_State::set_specifications()
{
    inversion_name = "adjoint_state_";
    inversion_method = "Adjoint-State First-Arrival Tomography";
    
    source_grad = new float[modeling->matsize]();
    source_comp = new float[modeling->matsize]();
    
    adjoint_grad = new float[modeling->matsize]();
    adjoint_comp = new float[modeling->matsize]();

    gradient = new float[modeling->nPoints]();
    illumination = new float[modeling->nPoints]();     

    cell_area = modeling->dx*modeling->dz;
}

void Adjoint_State::apply_inversion_technique()
{
    initialization();

    for (i = 1; i < modeling->nzz; i++)
    {
        for (j = 1; j < modeling->nxx; j++) 
            inner_sweep();

        for (j = modeling->nxx-2; j >= 0; j--) 
            inner_sweep();
    }

    for (i = modeling->nzz-2; i >= 0; i--)
    {
        for (j = 1; j < modeling->nxx; j++)
            inner_sweep();

        for (j = modeling->nxx-2; j >= 0; j--)
            inner_sweep();
    }

    for (j = 1; j < modeling->nxx; j++)
    {        
        for (i = 1; i < modeling->nzz; i++)
            inner_sweep();

        for (i = modeling->nzz-2; i >= 0; i--)
            inner_sweep();
    }

    for (j = modeling->nxx - 2; j >= 0; j--)
    {
        for (i = 1; i < modeling->nzz; i++)
            inner_sweep();

        for (i = modeling->nzz-2; i >= 0; i--)
            inner_sweep();
    }

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

    int sIdx = (int)(modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dx);
    int sIdz = (int)(modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dz);

    float So = modeling->S[sIdz + sIdx*modeling->nzz];

    int spreadId = 0;

    for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    {
        int rIdx = (int)(modeling->geometry->xrec[modeling->recId] / modeling->dx) + modeling->nb;
        int rIdz = (int)(modeling->geometry->zrec[modeling->recId] / modeling->dz) + modeling->nb;

        float X = sqrtf(powf(rIdx - sIdx, 2.0f) + powf(rIdz - sIdz, 2.0f));

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int xi = rIdx + (j - 1);
                int zi = rIdz + (i - 1);

                source_grad[zi + xi*modeling->nzz] += (dobs[spreadId + skipped] - modeling->T[zi + xi*modeling->nzz]) / cell_area;
                source_comp[zi + xi*modeling->nzz] += 1.0f / (X*X*So);
            }
        }

        ++spreadId;
    }   
} 

void Adjoint_State::inner_sweep()
{
    if ((i > 0) && (i < modeling->nzz-1) && (j > 0) && (j < modeling->nxx-1))
    {
        float a1 = -1.0f*(modeling->T[i + j*modeling->nzz] - modeling->T[i + (j-1)*modeling->nzz]) / modeling->dx;
        float ap1 = (a1 + fabsf(a1)) / 2.0f;
        float am1 = (a1 - fabsf(a1)) / 2.0f;

        float a2 = -1.0f*(modeling->T[i + (j+1)*modeling->nzz] - modeling->T[i + j*modeling->nzz]) / modeling->dx;
        float ap2 = (a2 + fabsf(a2)) / 2.0f;
        float am2 = (a2 - fabsf(a2)) / 2.0f;

        float c1 = -1.0f*(modeling->T[i + j*modeling->nzz] - modeling->T[(i-1) + j*modeling->nzz]) / modeling->dz;
        float cp1 = (c1 + fabsf(c1)) / 2.0f;
        float cm1 = (c1 - fabsf(c1)) / 2.0f;

        float c2 = -1.0f*(modeling->T[(i+1) + j*modeling->nzz] - modeling->T[i + j*modeling->nzz]) / modeling->dz;
        float cp2 = (c2 + fabsf(c2)) / 2.0f;
        float cm2 = (c2 - fabsf(c2)) / 2.0f;

        float d = (ap2 - am1)/modeling->dx + (cp2 - cm1)/modeling->dz;

        if (fabsf(d) < 1e-6f)
        {
            adjoint_grad[i + j*modeling->nzz] = 0.0f;    
            adjoint_comp[i + j*modeling->nzz] = 0.0f;    
        }
        else
        {
            float eg = (ap1*adjoint_grad[i + (j-1)*modeling->nzz] - am2*adjoint_grad[i + (j+1)*modeling->nzz]) / modeling->dx +
                       (cp1*adjoint_grad[(i-1) + j*modeling->nzz] - cm2*adjoint_grad[(i+1) + j*modeling->nzz]) / modeling->dz;

            float ec = (ap1*adjoint_comp[i + (j-1)*modeling->nzz] - am2*adjoint_comp[i + (j+1)*modeling->nzz]) / modeling->dx +
                       (cp1*adjoint_comp[(i-1) + j*modeling->nzz] - cm2*adjoint_comp[(i+1) + j*modeling->nzz]) / modeling->dz;

            float fg = (eg + source_grad[i + j*modeling->nzz]) / d;
            float fc = (ec + source_comp[i + j*modeling->nzz]) / d;

            if (adjoint_grad[i + j*modeling->nzz] > fg) 
                adjoint_grad[i + j*modeling->nzz] = fg;

            if (adjoint_comp[i + j*modeling->nzz] > fc) 
                adjoint_comp[i + j*modeling->nzz] = fc;
        }
    }
}

void Adjoint_State::optimization()
{   
    gradient_preconditioning();

    float gmax = 0.0f;
    float gdot = 0.0f;

    for (int index = 0; index < modeling->nPoints; index++)
    {
        if (gmax < fabsf(gradient[index]))
            gmax = fabsf(gradient[index]);

        gdot += gradient[index]*gradient[index];
    }

    float gamma = max_slowness_variation;

    float lambda = 0.5f * residuo.back() / gdot;     
    
    float alpha = (lambda*gmax > gamma) ? (gamma / (lambda*gmax)) : 1.0f; 

    for (int index = 0; index < modeling->nPoints; index++)
        perturbation[index] = alpha*lambda*gradient[index];        
}

void Adjoint_State::gradient_preconditioning()
{
    float sigx = 0.01 * M_PI / modeling->dx;
    float sigz = 0.01 * M_PI / modeling->dz;

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

    for (int i = 0; i < modeling->nz; i++) 
    {
        for (int j = 0; j < modeling->nx; j++) 
        {
            double gaussian_weight = 1.0 - exp(-(0.5*pow(kx[j] / sigx, 2.0) + 0.5*pow(kz[i] / sigz, 2.0)));

            output[i + j*modeling->nz][0] *= gaussian_weight; 
            output[i + j*modeling->nz][1] *= gaussian_weight; 
        }
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


