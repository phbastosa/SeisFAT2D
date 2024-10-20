# include "tomography.hpp"

void Tomography::set_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));

    obs_data_folder = catch_parameter("obs_data_folder", parameters);
    obs_data_prefix = catch_parameter("obs_data_prefix", parameters);

    smooth_model_per_iteration = str2bool(catch_parameter("smooth_per_iteration", parameters));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", parameters));
    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", parameters));
    
    convergence_map_folder = catch_parameter("convergence_folder", parameters);
    estimated_model_folder = catch_parameter("inversion_output_folder", parameters);

    write_model_per_iteration = str2bool(catch_parameter("export_model_per_iteration", parameters));

    set_forward_modeling();
    set_inversion_elements();
    set_specifications();
}

void Tomography::set_forward_modeling()
{
    modeling = new Eikonal_Iso();

    modeling->parameters = parameters;

    modeling->set_parameters();
}

void Tomography::set_inversion_elements()
{    
    ndata = 0;
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
        ndata += modeling->geometry->spread[shot];

    dcal = new float[ndata]();
    dobs = new float[ndata]();

    perturbation = new float[modeling->nPoints]();
}

void Tomography::import_obsData()
{
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        float * data = new float[modeling->geometry->spread[shot]]();

        std::string path = obs_data_folder + obs_data_prefix + std::to_string(modeling->geometry->sInd[shot]+1) + ".bin";

        import_binary_float(path, data, modeling->geometry->spread[shot]);

        int skipped = shot * modeling->geometry->spread[shot];    
        for (int i = 0; i < modeling->geometry->spread[shot]; i++) 
            dobs[i + skipped] = data[i];

        delete[] data;
    }
}

void Tomography::forward_modeling()
{
    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        modeling->srcId = shot;
    
        modeling->show_information();

        show_information();

        modeling->initialization();
        modeling->forward_solver();

        concatenate_data();
        
        if (iteration != max_iteration)
            apply_inversion_technique();
    }
}

void Tomography::show_information()
{
    std::cout << "\nInversion type: ";
    std::cout << inversion_method << "\n\n";

    if (iteration == max_iteration) 
        std::cout << "-------- Checking final residuo --------\n\n";
    else
    {    
        std::cout << "-------- Computing iteration " << iteration + 1 << " of " << max_iteration << " --------\n\n";

        if (iteration > 0) std::cout << "Previous residuo: " << residuo.back() << "\n\n";   
    }
}

void Tomography::concatenate_data()
{
    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    for (int i = 0; i < modeling->geometry->spread[modeling->srcId]; i++) 
        dcal[i + skipped] = modeling->synthetic_data[i];    
}

void Tomography::check_convergence()
{
    float square_difference = 0.0f;

    for (int i = 0; i < ndata; i++)
        square_difference += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(square_difference));

    if ((iteration >= max_iteration))
    {
        std::cout << "Final residuo: "<< residuo.back() <<"\n";
        converged = true;
    }
    else
    {
        iteration += 1;
        converged = false;
    }
}

void Tomography::smooth_matrix(float * input, float * output, int nx, int nz)
{
    int init = smoother_samples / 2;
    int nPoints = nx * nz;
    int nKernel = smoother_samples * smoother_samples;

    float pi = 4.0f * atanf(1.0f); 

    float * kernel = new float[nKernel]();

    for (int i = 0; i < nPoints; i++) 
        output[i] = input[i];

    int mid = (int)(smoother_samples / 2); 

    kernel[mid + mid*smoother_samples] = 1.0f;

    if (smoother_stdv != 0.0f)
    {
        float sum = 0.0f;

        for (int x = -init; x <= init; x++)
        {
            for (int z = -init; z <= init; z++)
            {          
                int index = (z + init) + (x + init)*smoother_samples; 
                
                float r = sqrtf(x*x + z*z);

                kernel[index] = 1.0f / (pi*smoother_stdv) * expf(-((r*r)/(2.0f*smoother_stdv*smoother_stdv)));
    
                sum += kernel[index]; 
            }
        }

        for (int i = 0; i < nKernel; i++) 
            kernel[i] /= sum;
    }
        
    for (int j = init; j < nx - init; j++)
    {
        for (int i = init; i < nz - init; i++)
        {       
            float accum = 0.0f;
                
            for (int xk = 0; xk < smoother_samples; xk++)
            {      
                for (int zk = 0; zk < smoother_samples; zk++)
                {   
                    int index = zk + xk*smoother_samples;   
                    int partial = (i - init + zk) + (j - init + xk)*nz; 

                    accum += input[partial] * kernel[index];
                }        
            }
                
            output[i + j*nz] = accum;
        }   
    }

    delete[] kernel;
}

void Tomography::model_update()
{
    if (smooth_model_per_iteration)
    {
        int aux_nx = modeling->nx + 2*smoother_samples;
        int aux_nz = modeling->nz + 2*smoother_samples;

        int aux_nPoints = aux_nx*aux_nz;

        float * dm_aux = new float[aux_nPoints]();
        float * dm_smooth = new float[aux_nPoints]();

        for (int index = 0; index < modeling->nPoints; index++)
        {
            int i = (int) (index % modeling->nz);    
            int j = (int) (index / modeling->nz);  

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz;

            dm_aux[ind_filt] = perturbation[i + j*modeling->nz];
        }

        smooth_matrix(dm_aux, dm_smooth, aux_nx, aux_nz);

        for (int index = 0; index < modeling->nPoints; index++)
        {
            int i = (int) (index % modeling->nz);    
            int j = (int) (index / modeling->nz);  

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz;

            perturbation[i + j*modeling->nz] = dm_smooth[ind_filt];
        }
    
        delete[] dm_aux;
        delete[] dm_smooth;
    }   
    
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz); 

        modeling->S[(i + modeling->nb) + (j + modeling->nb)*modeling->nzz] += perturbation[index];
        modeling->Vp[index] = 1.0f / modeling->S[(i + modeling->nb) + (j + modeling->nb)*modeling->nzz];
    }

    if (write_model_per_iteration)
    {
        std::string model_iteration_path = estimated_model_folder + "model_iteration_" + std::to_string(iteration) + "_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";

        export_binary_float(model_iteration_path, modeling->Vp, modeling->nPoints);
    }
}

void Tomography::export_results()
{    
    std::string estimated_model_path = estimated_model_folder + "final_model_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
    std::string convergence_map_path = convergence_map_folder + "convergence_" + std::to_string(iteration) + "_iterations.txt"; 

    export_binary_float(estimated_model_path, modeling->Vp, modeling->nPoints);

    std::ofstream resFile(convergence_map_path, std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) 
        resFile << residuo[r] << "\n";

    resFile.close();

    std::cout << "Text file \033[34m" << convergence_map_path << "\033[0;0m was successfully written." << std::endl;
}