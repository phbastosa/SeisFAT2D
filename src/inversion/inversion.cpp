# include "inversion.hpp"

void Inversion::set_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));

    obs_data_folder = catch_parameter("obs_data_folder", parameters);
    obs_data_prefix = catch_parameter("obs_data_prefix", parameters);

    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", parameters));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", parameters));
    
    convergence_map_folder = catch_parameter("convergence_folder", parameters);
    estimated_model_folder = catch_parameter("inversion_output_folder", parameters);

    write_model_per_iteration = str2bool(catch_parameter("export_model_per_iteration", parameters));
    smooth_model_per_iteration = str2bool(catch_parameter("smooth_per_iteration", parameters));

    set_forward_modeling();
}

void Inversion::import_obsData()
{
    n_data = modeling->max_spread * modeling->geometry->nrel;

    dcal = new float[n_data]();
    dobs = new float[n_data]();

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        float * data = new float[modeling->geometry->spread[modeling->srcId]]();

        std::string path = obs_data_folder + obs_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

        import_binary_float(path, data, modeling->geometry->spread[modeling->srcId]);

        int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];    
        
        for (int i = 0; i < modeling->geometry->spread[modeling->srcId]; i++) 
            dobs[i + skipped] = data[i];

        delete[] data;
    }
}

// void Tomography::forward_modeling()
// {
//     for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
//     {
//         show_information();

//         modeling->initialization();
//         modeling->forward_solver();

//         concatenate_data();
        
//         if (iteration != max_iteration)
//             apply_inversion_technique();
//     }
// }

// void Tomography::show_information()
// {
//     modeling->show_information();    
    
//     std::cout << "\nInversion type: " << inversion_method << "\n\n";

//     if (iteration == max_iteration) 
//         std::cout << "-------- Checking final residuo --------\n\n";
//     else
//     {    
//         std::cout << "-------- Computing iteration " << iteration + 1 << " of " << max_iteration << " --------\n\n";

//         if (iteration > 0) std::cout << "Previous residuo: " << residuo.back() << "\n\n";   
//     }
// }

// void Tomography::concatenate_data()
// {
//     int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

//     for (int i = 0; i < modeling->geometry->spread[modeling->srcId]; i++) 
//         dcal[i + skipped] = modeling->synthetic_data[i];    
// }

// void Tomography::check_convergence()
// {
//     float square_difference = 0.0f;

//     for (int i = 0; i < n_data; i++)
//         square_difference += powf(dobs[i] - dcal[i], 2.0f);

//     residuo.push_back(sqrtf(square_difference));

//     if ((iteration >= max_iteration))
//     {
//         std::cout << "Final residuo: "<< residuo.back() <<"\n";
//         converged = true;
//     }
//     else
//     {
//         iteration += 1;
//         converged = false;
//     }
// }

// void Tomography::model_update()
// {
//     if (smooth_model_per_iteration)
//     {
//         int aux_nx = modeling->nx + 2*smoother_samples;
//         int aux_nz = modeling->nz + 2*smoother_samples;

//         int aux_nPoints = aux_nx*aux_nz;

//         float * dm_aux = new float[aux_nPoints]();
//         float * dm_smooth = new float[aux_nPoints]();

//         # pragma omp parallel for
//         for (int index = 0; index < modeling->nPoints; index++)
//         {
//             int i = (int) (index % modeling->nz);    
//             int j = (int) (index / modeling->nz);  

//             int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz;

//             dm_aux[ind_filt] = perturbation[i + j*modeling->nz];
//         }

//         smooth_matrix(dm_aux, dm_smooth, aux_nx, aux_nz);

//         # pragma omp parallel for    
//         for (int index = 0; index < modeling->nPoints; index++)
//         {
//             int i = (int) (index % modeling->nz);    
//             int j = (int) (index / modeling->nz);  

//             int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz;

//             perturbation[i + j*modeling->nz] = dm_smooth[ind_filt];
//         }
    
//         delete[] dm_aux;
//         delete[] dm_smooth;
//     }   
    
//     for (int index = 0; index < modeling->nPoints; index++)
//     {
//         int i = (int) (index % modeling->nz);    
//         int j = (int) (index / modeling->nz); 

//         int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

//         if ((i > 0) && (i < modeling->nz - 1) && (j > 0) && (j < modeling->nx - 1))
//             modeling->S[indb] += perturbation[index];

//         modeling->Vp[index] = 1.0f / modeling->S[indb];
//     }

//     if (write_model_per_iteration)
//     {
//         std::string model_iteration_path = estimated_model_folder + inversion_name + "model_iteration_" + std::to_string(iteration) + "_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";

//         export_binary_float(model_iteration_path, modeling->Vp, modeling->nPoints);
//     }
// }

// void Tomography::smooth_matrix(float * input, float * output, int nx, int nz)
// {
//     int nPoints = nx * nz;
//     int nKernel = smoother_samples * smoother_samples;

//     float pi = 4.0f * atanf(1.0f); 

//     float * kernel = new float[nKernel]();

//     # pragma omp parallel for
//     for (int i = 0; i < nPoints; i++) 
//         output[i] = input[i];

//     int mid = (int)(smoother_samples / 2); 

//     kernel[mid + mid*smoother_samples] = 1.0f;

//     if (smoother_stdv != 0.0f)
//     {
//         float sum = 0.0f;

//         for (int x = -mid; x <= mid; x++)
//         {
//             for (int z = -mid; z <= mid; z++)
//             {          
//                 int index = (z + mid) + (x + mid)*smoother_samples; 
                
//                 float r = sqrtf(x*x + z*z);

//                 kernel[index] = 1.0f / (pi*smoother_stdv) * expf(-((r*r)/(2.0f*smoother_stdv*smoother_stdv)));
    
//                 sum += kernel[index]; 
//             }
//         }

//         for (int i = 0; i < nKernel; i++) 
//             kernel[i] /= sum;
//     }
        
//     for (int j = mid; j < nx - mid; j++)
//     {
//         for (int i = mid; i < nz - mid; i++)
//         {       
//             float accum = 0.0f;
                
//             for (int xk = 0; xk < smoother_samples; xk++)
//             {      
//                 for (int zk = 0; zk < smoother_samples; zk++)
//                 {   
//                     int index = zk + xk*smoother_samples;   
//                     int partial = (i - mid + zk) + (j - mid + xk)*nz; 

//                     accum += input[partial] * kernel[index];
//                 }        
//             }
                
//             output[i + j*nz] = accum;
//         }   
//     }

//     delete[] kernel;
// }

// void Tomography::export_results()
// {    
//     std::string estimated_model_path = estimated_model_folder + inversion_name + "final_model_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
//     std::string convergence_map_path = convergence_map_folder + inversion_name + "convergence_" + std::to_string(iteration) + "_iterations.txt"; 

//     export_binary_float(estimated_model_path, modeling->Vp, modeling->nPoints);

//     std::ofstream resFile(convergence_map_path, std::ios::out);
    
//     for (int r = 0; r < residuo.size(); r++) 
//         resFile << residuo[r] << "\n";

//     resFile.close();

//     std::cout << "Text file \033[34m" << convergence_map_path << "\033[0;0m was successfully written." << std::endl;
// }

// void Least_Squares::set_specifications()
// {
//     inversion_name = "least_squares_";
//     inversion_method = "Least-Squares First-Arrival Tomography";

//     tk_order = std::stoi(catch_parameter("tk_order", parameters));
//     tk_param = std::stof(catch_parameter("tk_param", parameters));

//     n_model = modeling->nPoints;

//     ray_path_max_samples = 0;

//     for (int shot = 0; shot < modeling->geometry->nrel; shot++)
//     {
//         for (int node = modeling->geometry->iRec[shot]; node < modeling->geometry->fRec[shot]; node++)
//         {
//             float dx = (modeling->geometry->xsrc[modeling->geometry->sInd[shot]] - modeling->geometry->xrec[node]) / modeling->dx;
//             float dz = (modeling->geometry->zsrc[modeling->geometry->sInd[shot]] - modeling->geometry->zrec[node]) / modeling->dz;
            
//             ray_path_max_samples += (size_t)(sqrtf(dx*dx + dz*dz));
//         }
//     }

//     iG.reserve(ray_path_max_samples);
//     jG.reserve(ray_path_max_samples);
//     vG.reserve(ray_path_max_samples);
// }

// void Least_Squares::apply_inversion_technique()
// {
//     int sIdx = (int)(modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dx);
//     int sIdz = (int)(modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] / modeling->dz);

//     int sId = sIdz + sIdx*modeling->nz; 

//     float rayStep = 0.2f * modeling->dz;

//     std::vector < int > ray_index; 

//     for (int ray_id = modeling->geometry->iRec[modeling->srcId]; ray_id < modeling->geometry->fRec[modeling->srcId]; ray_id++)
//     {
//         float xi = modeling->geometry->xrec[ray_id];        
//         float zi = modeling->geometry->zrec[ray_id];

//         if ((modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] == zi) && 
//             (modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] == xi))
//             continue;        

//         while (true)
//         {
//             int j = (int)(xi / modeling->dx) + modeling->nb;
//             int i = (int)(zi / modeling->dz) + modeling->nb;

//             float dTx = (modeling->T[i + (j+1)*modeling->nzz] - modeling->T[i + (j-1)*modeling->nzz]) / (2.0f*modeling->dx);    
//             float dTz = (modeling->T[(i+1) + j*modeling->nzz] - modeling->T[(i-1) + j*modeling->nzz]) / (2.0f*modeling->dz);    

//             float norm = sqrtf(dTx*dTx + dTz*dTz);

//             xi -= rayStep*dTx / norm;   
//             zi -= rayStep*dTz / norm;    

//             int jm = (int)(xi / modeling->dx); 
//             int im = (int)(zi / modeling->dz); 

//             int index = im + jm*modeling->nz;
            
//             ray_index.push_back(index);

//             if (ray_index.back() == sId) break;
//         }
   
//         float final_distance = sqrtf(powf(zi - modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]], 2.0f) + 
//                                      powf(xi - modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]], 2.0f));

//         std::sort(ray_index.begin(), ray_index.end());

//         int current_voxel_index = ray_index[0];
//         float distance_per_voxel = rayStep;

//         for (int index = 0; index < ray_index.size(); index++)
//         {
//             if (ray_index[index] == current_voxel_index)
//             {
//                 distance_per_voxel += rayStep;
//             }
//             else
//             {
//                 vG.emplace_back(distance_per_voxel);
//                 jG.emplace_back(current_voxel_index);
//                 iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);

//                 if (current_voxel_index == sId) vG.back() = final_distance;

//                 distance_per_voxel = rayStep;
//                 current_voxel_index = ray_index[index];    
//             }
//         }

//         if (current_voxel_index == sId)
//         {
//             vG.emplace_back(final_distance);
//             jG.emplace_back(current_voxel_index);
//             iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);
//         }
//         else 
//         {
//             vG.emplace_back(distance_per_voxel);
//             jG.emplace_back(current_voxel_index);
//             iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);
//         }

//         std::vector<int>().swap(ray_index);
//     }
// }

// void Least_Squares::optimization()
// {
//     M = n_model;                                  
//     N = n_data + n_model - tk_order;                    
//     NNZ = vG.size() + (tk_order + 1) * (n_model - tk_order);

//     iA = new int[NNZ]();
//     jA = new int[NNZ]();
//     vA = new float[NNZ]();

//     B = new float[N]();
//     x = new float[M]();

//     for (int index = 0; index < n_data; index++) 
//         B[index] = dobs[index] - dcal[index];

//     for (int index = 0; index < vG.size(); index++)
//     {
//         iA[index] = iG[index];
//         jA[index] = jG[index];
//         vA[index] = vG[index];
//     }

//     std::vector< int >().swap(iG);
//     std::vector< int >().swap(jG);
//     std::vector<float>().swap(vG);

//     apply_regularization();
//     solve_linear_system_lscg();

//     # pragma omp parallel for	
//     for (int index = 0; index < n_model; index++)
//         perturbation[index] = x[index];

//     delete[] x;
//     delete[] B;
//     delete[] iA;
//     delete[] jA;
//     delete[] vA;    
// }

// void Least_Squares::apply_regularization()
// {
//     int elements = tk_order + 1;
		
//     int n = n_model - tk_order;
//     int nnz = elements * n;	
    
//     int * iL = new int[nnz]();
//     int * jL = new int[nnz]();
//     float * vL = new float[nnz]();

//     if (tk_order <= 0)
//     {
// 	for (int index = 0; index < nnz; index++)
//         {
//             iL[index] = index;
// 	    jL[index] = index;
// 	    vL[index] = 1.0f;
// 	}
//     } 
//     else
//     {
//         int * df = new int[elements]();	
//         int * df1 = new int[elements + 1]();
//         int * df2 = new int[elements + 1]();
        
//         df[0] = -1; df[1] = 1;
        
//         for (int index = 1; index < tk_order; index++)
//         {
//             for (int k = 0; k < elements; k++)
//             {
//                 df2[k] = df[k];
//                 df1[k + 1] = df[k];

//                 df[k] = df1[k] - df2[k]; 
//             }		 
//         }
        
//         for (int index = 0; index < n; index++)
//         {
//             for (int k = 0; k < elements; k++)
//             {
//                 iL[elements*index + k] = index;	
//                 jL[elements*index + k] = index + k;
//                 vL[elements*index + k] = df[k];
//             }	
//         }

//         delete[] df;
//         delete[] df1;
//         delete[] df2;
//     }

//     for (int index = NNZ - nnz; index < NNZ; index++) 
//     {
//         iA[index] = n_data + iL[index - (NNZ - nnz)];
//         jA[index] = jL[index - (NNZ - nnz)];
//         vA[index] = tk_param * vL[index - (NNZ - nnz)];        
//     }

//     delete[] iL;
//     delete[] jL;
//     delete[] vL;
// }

// void Least_Squares::solve_linear_system_lscg()
// {
//     float a, b, qTq, rTr, rd;
//     int cg_max_iteration = 10;

//     float * s = new float[N]();
//     float * q = new float[N]();
//     float * r = new float[M]();
//     float * p = new float[M]();

//     // s = d - G * x, where d = dobs - dcal and x = slowness variation
//     for (int i = 0; i < N; i++) 
//         s[i] = B[i]; 

//     // r = G' * s    
//     for (int i = 0; i < NNZ; i++) 
//         r[jA[i]] += vA[i] * s[iA[i]];        

//     // p = r and x = 0;
//     for (int i = 0; i < M; i++) 
//     {
//         p[i] = r[i]; 
//         x[i] = 0.0f;
//     }

//     // q = G * p
//     for (int i = 0; i < NNZ; i++) 
//         q[iA[i]] += vA[i] * p[jA[i]];        

//     for (int i = 0; i < cg_max_iteration; i++)
//     {
//         qTq = 0.0f;
//         for (int k = 0; k < N; k++)           // q inner product
//             qTq += q[k] * q[k];               // qTq = q' * q

//         rTr = 0.0f;
//         for (int k = 0; k < M; k++)           // r inner product
//             rTr += r[k] * r[k];               // rTr = r' * r 

//         a = rTr / qTq;                        // a = (r' * r) / (q' * q)                    

//         for (int k = 0; k < M; k++)           // model atualization
//             x[k] += a * p[k];                 // x = x + a * p

//         for (int k = 0; k < N; k++)           // s atualization  
//             s[k] -= a * q[k];                 // s = s - a * q 

//         rd = 0.0f;
//         for (int k = 0; k < M; k++)           // r inner product for division 
//             rd += r[k] * r[k];                // rd = r' * r

//         for (int k = 0; k < M; k++)           // Zeroing r 
//             r[k] = 0.0f;                      // r = 0, for multiplication
        
//         for (int k = 0; k < NNZ; k++)         // r atualization 
//             r[jA[k]] += vA[k] * s[iA[k]];     // r = G' * s    

//         rTr = 0.0f;                
//         for (int k = 0; k < M; k++)           // r inner product
//             rTr += r[k] * r[k];               // rTr = r' * r

//         b = rTr / rd;                         // b = (r' * r) / rd

//         for (int k = 0; k < M; k++)          
//             p[k] = r[k] + b * p[k];           // p = r + b * p 

//         for (int k = 0; k < N; k++) 
//             q[k] = 0.0f;                      // q = 0, for multiplication

//         for (int k = 0; k < NNZ; k++) 
//             q[iA[k]] += vA[k] * p[jA[k]];     // q = G * p   
//     }
// }
