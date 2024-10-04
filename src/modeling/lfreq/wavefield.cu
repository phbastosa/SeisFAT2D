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

    double * time_domain = (double *) fftw_malloc(nt*sizeof(double));

    fftw_complex * freq_domain = (fftw_complex *) fftw_malloc(nt*sizeof(fftw_complex));

    fftw_plan forward_plan = fftw_plan_dft_r2c_1d(nt, time_domain, freq_domain, FFTW_ESTIMATE);
    fftw_plan inverse_plan = fftw_plan_dft_c2r_1d(nt, freq_domain, time_domain, FFTW_ESTIMATE);

    double delta_f = 1.0 / (nt * dt);  
    
    std::complex<double> I(0.0, 1.0);  

    for (int k = 0; k < nt; k++)
    {
        time_domain[k] = (double) signal_aux2[k];
    }

    fftw_execute(forward_plan);

    for (int k = 0; k < nt; ++k) 
    {
        double freq = (k <= nt / 2) ? k * delta_f : (k - nt) * delta_f;
        
        std::complex<double> factor = std::pow(I * freq, 0.5);  

        std::complex<double> complex_freq(freq_domain[k][0], freq_domain[k][1]);
        std::complex<double> filtered_freq = complex_freq * factor;

        freq_domain[k][0] = filtered_freq.real();
        freq_domain[k][1] = filtered_freq.imag();
    }

    fftw_execute(inverse_plan);    

    for (int k = 0; k < nt; k++)
    {
        signal_aux1[k] = (float) time_domain[k];
    }

    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal_aux1, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] signal_aux1;
    delete[] signal_aux2;
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

