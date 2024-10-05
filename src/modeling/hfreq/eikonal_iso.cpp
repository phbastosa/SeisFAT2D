# include "eikonal_iso.hpp"

void Eikonal_Iso::set_properties()
{
    float * vp = new float[nPoints]();

    std::string model_file = catch_parameter("vp_model_file", parameters);

    import_binary_float(model_file, vp, nPoints);

    for (int index = 0; index < nPoints; index++)
        vp[index] = 1.0f / vp[index];

    S = new float[matsize]();

    expand_boundary(vp, S);

    delete[] vp;
}

void Eikonal_Iso::set_conditions()
{
    dz2i = 1.0f / (dz * dz);
    dx2i = 1.0f / (dx * dx);

    T = new float[matsize]();
}

void Eikonal_Iso::forward_solver()
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

    compute_seismogram();
}

void Eikonal_Iso::inner_sweep()
{
    i1 = i - sgnvz;
    j1 = j - sgnvx;

    tv = T[(i - sgntz) + j*nzz];
    te = T[i + (j - sgntx)*nzz];
    tev = T[(i - sgntz) + (j - sgntx)*nzz];

    Sref = std::min(S[i1 + std::max(j - 1, 1)*nzz], S[i1 + std::min(j, nxx - 1)*nzz]);
    t1d1 = tv + dz*Sref;

    Sref = std::min(S[std::max(i - 1, 1) + j1*nzz], S[std::min(i, nzz - 1) + j1*nzz]);
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

    T[i + j*nzz] = std::min(T[i + j*nzz], std::min(t1D, t2D));
}


