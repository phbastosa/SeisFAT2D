# include "wavefield.cuh"

void Wavefield::set_specifications()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));

    fmax = std::stof(catch_parameter("max_frequency", parameters));

    set_wavelet();
    set_boundaries();
    set_properties();    
    set_conditions();    

    nThreads = 256;
    nBlocks = (int)(matsize / nThreads) + 1;

    current_xrec = new int[max_spread]();
    current_zrec = new int[max_spread]();

    cudaMalloc((void**)&(rIdx), max_spread*sizeof(int));
    cudaMalloc((void**)&(rIdz), max_spread*sizeof(int));
}

void Wavefield::set_wavelet()
{
    float * aux_s = new float[nt]();
    float * signal = new float[nt]();

    float pi = 4.0f*atanf(1.0f);
    float t0 = 2.0f*sqrtf(pi) / fmax;
    float fc = fmax / (3.0f * sqrtf(pi));

    tlag = (int)(t0 / dt) + 1;

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*pi*fc*fc*td*td;

        aux_s[n] = 1e5f*(1.0f - 2.0f*arg)*expf(-arg);
    }

    for (int n = 0; n < nt; n++)
    {
        float summation = 0;
        for (int i = 0; i < n; i++)
            summation += aux_s[i];    
        
        signal[n] = summation;
    }

    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] aux_s;
    delete[] signal;
}

void Wavefield::set_boundaries()
{
    nb = std::stoi(catch_parameter("boundary_samples", parameters));

    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    matsize = nxx*nzz;
}
