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

    define_cerjan_dampers();

    cudaMalloc((void**)&(rIdx), max_spread*sizeof(int));
    cudaMalloc((void**)&(rIdz), max_spread*sizeof(int));
}

void Wavefield::set_wavelet()
{
    float * signal_aux1 = new float[nt]();
    float * signal_aux2 = new float[nt]();

    float pi = 4.0f*atanf(1.0f);
    float t0 = 2.0f*sqrtf(pi) / fmax;
    float fc = fmax / (3.0f * sqrtf(pi));

    tlag = (int)(t0 / dt) + 1;

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*pi*fc*fc*td*td;

        signal_aux1[n] = 1e5f*(1.0f - 2.0f*arg)*expf(-arg);
    }

    for (int n = 0; n < nt; n++)
    {
        float summation = 0;
        for (int i = 0; i < n; i++)
            summation += signal_aux1[i];    
        
        signal_aux2[n] = summation;
    }

    


    // export_binary_float("wavelet_original.bin", signal_aux2, nt);
    // export_binary_float("wavelet_modified.bin", signal_aux1, nt);



    // cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    // cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);

    // delete[] aux_s;
    // delete[] signal;
}

void Wavefield::set_boundaries()
{
    nb = std::stoi(catch_parameter("boundary_samples", parameters));

    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    matsize = nxx*nzz;
}

void Wavefield::define_cerjan_dampers()
{
    float * damp1D = new float[nb]();
    float * damp2D = new float[nb*nb]();

    float factor = std::stof(catch_parameter("boundary_damper", parameters));

    for (int i = 0; i < nb; i++) 
    {
        damp1D[i] = expf(-powf(factor * (nb - i), 2.0f));
    }

    for(int i = 0; i < nb; i++) 
    {
        for (int j = 0; j < nb; j++)
        {   
            damp2D[j + i*nb] += damp1D[i]; // up to bottom
            damp2D[i + j*nb] += damp1D[i]; // left to right
        }
    }

    for (int index = 0; index < nb*nb; index++)
        damp2D[index] -= 1.0f;

	cudaMalloc((void**)&(d1D), nb*sizeof(float));
	cudaMalloc((void**)&(d2D), nb*nb*sizeof(float));

	cudaMemcpy(d1D, damp1D, nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d2D, damp2D, nb*nb*sizeof(float), cudaMemcpyHostToDevice);

    delete[] damp1D;
    delete[] damp2D;
}

