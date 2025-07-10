# include "tomography_iso.hpp"

void Tomography_ISO::set_modeling_type()
{
    modeling = new Eikonal_ISO();
    modeling->parameters = parameters;
    modeling->set_parameters();

    inversion_name = "tomography_iso";
    inversion_method = "Isotropic First-Arrival Tomography";

    dS = new float[modeling->nPoints]();
}

void Tomography_ISO::set_sensitivity_matrix()
{
    int nnz = (tk_order + 1) * (n_model - tk_order);                    

    M = n_model;                                  
    N = n_data;
    NNZ = vG.size();

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();
    x = new float[M]();

    for (int index = 0; index < n_data; index++) 
        B[index] = dobs[index] - dcal[index];

    for (int index = 0; index < vG.size(); index++)
    {
        iA[index] = iG[index];
        jA[index] = jG[index];
        vA[index] = vG[index];
    }

    for (int index = vG.size(); index < NNZ; index++)
    {
        iA[index] = n_data + iR[index - (NNZ - nnz)];
        jA[index] = jR[index - (NNZ - nnz)];
        vA[index] = tk_param * vR[index - (NNZ - nnz)];  
    }

    std::vector< int >().swap(iG);
    std::vector< int >().swap(jG);
    std::vector<float>().swap(vG);        

    delete[] iR;
    delete[] jR;
    delete[] vR;
}

void Tomography_ISO::get_parameter_variation()
{
    #pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        dS[index] = x[index];
}

void Tomography_ISO::model_update()
{
    model_smoothing(dS);

    for (int index = 0; index < modeling->nPoints; index++)
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz); 

        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

        modeling->S[indb] += dS[index];
    }

    modeling->copy_slowness_to_device();
}

void Tomography_ISO::export_estimated_models()
{
    float * Vp = new float[modeling->nPoints]();
    modeling->reduce_boundary(modeling->S, Vp);
    
    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        Vp[index] = 1.0f / Vp[index];

    std::string estimated_vp_path = estimated_model_folder + inversion_name + "_final_model_vp_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
    export_binary_float(estimated_vp_path, Vp, modeling->nPoints);
}