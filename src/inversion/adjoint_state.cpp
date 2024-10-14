# include "adjoint_state.hpp"

void Adjoint_State::set_specifications()
{
    inversion_method = "Adjoint-State First-Arrival Tomography";

    source = new float[modeling->matsize]();
    adjoint = new float[modeling->matsize]();
    gradient = new float[modeling->nPoints]();

    cell_area = modeling->dx * modeling->dz;
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

    gradient_preconditioning();
}

void Adjoint_State::initialization()
{
    for (int index = 0; index < modeling->matsize; index++) 
    {
        source[index] = 0.0f;    
        adjoint[index] = 1e6f;

        int i = (int) (index % modeling->nzz);    
        int j = (int) (index / modeling->nzz);  

        if ((i == 0) || (i == modeling->nzz - 1) || 
            (j == 0) || (j == modeling->nxx - 1)) 
        {    
            adjoint[index] = 0.0f;        
        }
    }

    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    int spreadId = 0;

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

                source[zi + xi*modeling->nzz] += (dobs[spreadId + skipped] - modeling->T[zi + xi*modeling->nzz]) / cell_area;
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
            adjoint[i + j*modeling->nzz] = 0.0f;    
        }
        else
        {
            float e = (ap1*adjoint[i + (j-1)*modeling->nzz] - am2*adjoint[i + (j+1)*modeling->nzz]) / modeling->dx +
                    (cp1*adjoint[(i-1) + j*modeling->nzz] - cm2*adjoint[(i+1) + j*modeling->nzz]) / modeling->dz;

            float f = (e + source[i + j*modeling->nzz]) / d;

            if (adjoint[i + j*modeling->nzz] > f) 
                adjoint[i + j*modeling->nzz] = f;
        }
    }
}

void Adjoint_State::gradient_preconditioning()
{
    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz);  

        int indp = i + j*modeling->nz; 
        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

        gradient[indp] += adjoint[indb]*modeling->T[indb]*modeling->S[indb]*modeling->S[indb]*cell_area / modeling->geometry->nrel;
    }
}

void Adjoint_State::optimization()
{



    
}
