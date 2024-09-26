# include "serial_aFSM.hpp"

void Serial_aFSM::set_specifications()
{
    dz2i = 1.0f / (dz * dz);
    dx2i = 1.0f / (dx * dx);
}

void Serial_aFSM::forward_solver()
{
    sgntz = 1; sgnvz = 1;
    for (i = 1; i < nzz; i++)
    {
        sgntx = 1; sgnvx = 1;
        for (j = 1; j < nxx; j++) 
            inner_sweep();

        sgntx = -1; sgnvx = 0;
        for (j = nxx-2; j >= 0; j--) 
            inner_sweep();
    }

    sgntz = -1; sgnvz = 0;
    for (i = nzz-2; i >= 0; i--)
    {
        sgntx = 1; sgnvx = 1;
        for (j = 1; j < nxx; j++)
            inner_sweep();

        sgntx = -1; sgnvx = 0;
        for (j = nxx-2; j >= 0; j--)
            inner_sweep();
    }

    sgntx = 1; sgnvx = 1;
    for (j = 1; j < nxx; j++)
    {        
        sgntz = 1; sgnvz = 1;
        for (i = 1; i < nzz; i++)
            inner_sweep();

        sgntz = -1; sgnvz = 0;
        for (i = nzz-2; i >= 0; i--)
            inner_sweep();
    }

    sgntx = -1; sgnvx = 0;
    for (j = nxx - 2; j >= 0; j--)
    {
        sgntz = 1; sgnvz = 1;
        for (i = 1; i < nzz; i++)
            inner_sweep();

        sgntz = -1; sgnvz = 0;
        for (i = nzz-2; i >= 0; i--)
            inner_sweep();
    }
}

void Serial_aFSM::inner_sweep()
{
    i1 = i - sgnvz;
    j1 = j - sgnvx;

    tv = eikonalT[(i - sgntz) + j*nzz];
    te = eikonalT[i + (j - sgntx)*nzz];
    tev = eikonalT[(i - sgntz) + (j - sgntx)*nzz];

    Sref = std::min(slowness[i1 + std::max(j - 1, 1)*nzz], slowness[i1 + std::min(j, nxx - 1)*nzz]);
    t1d1 = tv + dz*Sref;

    Sref = std::min(slowness[std::max(i - 1, 1) + j1*nzz], slowness[std::min(i, nzz - 1) + j1*nzz]);
    t1d2 = te + dx*Sref;

    t1D = std::min(t1d1, t1d2);

    t1 = t2 = t3 = 1e6f;

    if ((tv <= te + dx*Sref) && (te <= tv + dz*Sref) && (te - tev >= 0.0f) && (tv - tev >= 0.0f))
    {
        float ta = tev + te - tv;
        float tb = tev - te + tv;

        t1 = ((tb * dz2i + ta * dx2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
    }        
    else if ((te - tev <= dz*dz*Sref / sqrtf(dx*dx + dz*dz)) && (te - tev >= 0.0f))
    {    
        t2 = te + dx*sqrtf(Sref*Sref - ((te - tev) / dz)*((te - tev) / dz));
    }
    else if ((tv - tev <= dx*dx*Sref / sqrtf(dx*dx + dz*dz)) && (tv - tev >= 0.0f))
    {
        t3 = tv + dz*sqrtf(Sref*Sref - ((tv - tev) / dx)*((tv - tev) / dx));
    }

    t2D = std::min(t1, std::min(t2, t3));

    eikonalT[i + j*nzz] = std::min(eikonalT[i + j*nzz], std::min(t1D, t2D));
}