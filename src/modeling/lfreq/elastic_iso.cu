# include "elastic_iso.cuh"

void elastic_Iso::set_properties()
{
    std::string vp_file = catch_parameter("vp_model_file", parameters);
    std::string vs_file = catch_parameter("vs_model_file", parameters);
    std::string rho_file = catch_parameter("rho_model_file", parameters);

    float * vp = new float[nPoints]();
    float * vs = new float[nPoints]();
    float * rho = new float[nPoints]();

    Vp = new float[matsize]();
    Vs = new float[matsize]();
    Rho = new float[matsize]();

    import_binary_float(vp_file, vp, nPoints);
    import_binary_float(vs_file, vs, nPoints);
    import_binary_float(rho_file, rho, nPoints);

    expand_boundary(vp, Vp);
    expand_boundary(vs, Vs);
    expand_boundary(rho, Rho);

    delete[] vp;
    delete[] vs;
    delete[] rho;
}

void elastic_Iso::set_conditions()
{
    M = new float[matsize]();
    L = new float[matsize]();
    B = new float[matsize]();
    P = new float[matsize]();

    for (int index = 0; index < matsize; index++)
    {
        M[index] = Rho[index]*Vs[index]*Vs[index];
        L[index] = Rho[index]*Vp[index]*Vp[index] - 2.0f*M[index];
        B[index] = 1.0f / Rho[index];
    }

    cudaMalloc((void**)&(d_M), matsize*sizeof(float));
    cudaMalloc((void**)&(d_L), matsize*sizeof(float));
    cudaMalloc((void**)&(d_B), matsize*sizeof(float));
    cudaMalloc((void**)&(d_P), matsize*sizeof(float));
    
    cudaMalloc((void**)&(d_Vx), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Vz), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Txx), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Tzz), matsize*sizeof(float));
    cudaMalloc((void**)&(d_Txz), matsize*sizeof(float));

    cudaMemcpy(d_M, M, matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L, L, matsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, matsize*sizeof(float), cudaMemcpyHostToDevice);
}

void elastic_Iso::forward_solver()
{


}