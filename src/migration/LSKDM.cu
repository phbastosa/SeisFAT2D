# include "LSKDM.cuh"

void LSKDM::kirchhoff_depth_migration()
{
    // set_src_travel_times();
    // set_rec_travel_times();
    // prepare_convolution();

    // set_src_domain();

    // initialization();

    // while (true)
    // {
    //     compute_gradient();        
    //     compute_residuals();

    //     if (converged) break;

    //     regularization();

    //     compute_direction();
    //     compute_stepLength();

    //     update_model();
    // }
}

void LSKDM::export_outputs()
{


}









// void IDLSKDM::image_building()
// {
//     set_src_travel_times();
//     set_rec_travel_times();
//     prepare_convolution();
//     prepare_components();
//     set_seismic_data();

//     set_initial_model();

//     while (true)
//     {                             
//         check_convergence();      

//         if (converged) break;     
                                     
//         compute_gradient();      
//         compute_direction();    
//         compute_stepLength();   
//         update_reflectivity();  
//     }
// }

// void IDLSKDM::prepare_components()
// {
//     dobs = new float[nt*nTraces]();
//     dcal = new float[nt*nTraces]();
//     dres = new float[nt*nTraces]();
    
//     model = new float[modeling->matsize]();
    
//     gradient_p = new float[modeling->matsize]();
//     gradient_m = new float[modeling->matsize]();
    
//     direction = new float[modeling->matsize]();

//     IMAGE = new float[modeling->nPoints]();

//     iteration = 0;
//     converged = false;
// }

// void IDLSKDM::set_initial_model()
// {
//     current_operation = "IDLSKDM: adjoint operator -> Generating initial model";

//     adjoint(dobs, model);
// }

// void IDLSKDM::check_convergence()
// {
//     current_operation = "IDLSKDM: forward operator -> Checking convergence";

//     if (iteration == max_it) ++iteration;
    
//     forward(model, dcal);        

//     if (iteration > max_it) --iteration;

//     float squared_difference_sum = 0.0f;

//     for (int index = 0; index < nt*nTraces; index++)
//     {
//         dres[index] = dobs[index] - dcal[index];

//         squared_difference_sum += dres[index]*dres[index];
//     }

//     residuo.push_back(sqrtf(squared_difference_sum));

//     converged = (++iteration > max_it) ? true : false;

//     if (converged) std::cout << "Final residuo: "<< residuo.back() <<"\n";
// }

// void IDLSKDM::compute_gradient()
// {   
//     current_operation = "IDLSKDM: adjoint operator -> Computing gradient";

//     adjoint(dres, gradient_p);
    
//     float grad_norm = 0.0f;
//     for (int index = 0; index < modeling->matsize; index++)
//         grad_norm += gradient_p[index]*gradient_p[index];

//     grad_norm = sqrtf(grad_norm);

//     for (int index = 0; index < modeling->matsize; index++)
//         gradient_p[index] /= grad_norm;
// }

// void IDLSKDM::compute_direction()
// {
//     float sum1 = 0.0f;
//     float sum2 = 0.0f;

//     for (int index = 0; index < modeling->matsize; index++)
//     {
//         sum1 += gradient_p[index] * (gradient_p[index] - gradient_m[index]);
//         sum2 += gradient_m[index] * gradient_m[index];
//     }

//     beta = sum1 / (sum2 + EPS); 

//     for (int index = 0; index < modeling->matsize; index++)
//         direction[index] = beta*direction[index] - gradient_p[index];
// }

// void IDLSKDM::compute_stepLength()
// {
//     current_operation = "IDLSKDM: forward operator -> Computing step length";

//     forward(direction, dcal);

//     float sum1 = 0.0f; 
//     float sum2 = 0.0f;
    
//     for (int index = 0; index < nt*nTraces; index++)
//     {
//         sum1 += dres[index]*dcal[index];
//         sum2 += dcal[index]*dcal[index];
//     }
    
//     alpha = sum1 / (sum2 + EPS);
// }

// void IDLSKDM::update_reflectivity()
// {
//     for (int index = 0; index < modeling->matsize; index++)
//     {
//         model[index] = model[index] + alpha*direction[index];
        
//         gradient_m[index] = gradient_p[index];
//     }   
// }


// void IDLSKDM::show_iteration_info()
// {
//     if (iteration > max_it) 
//         std::cout << "\n-------- Checking final residuo --------\n\n";
//     else
//     {    
//         if (iteration == 0) 
//             std::cout << "\n-------- Computing first residuo --------\n";        
//         else
//         {
//             std::cout << "\n-------- Computing iteration " << iteration << " of " << max_it << " --------\n\n";
            
//             std::cout << "Previous residuo: " << residuo.back() << "\n";   
//         }
//     }
// }

// void IDLSKDM::export_outputs()
// {
//     std::string path = residuo_folder + "IDLSKDM_convergence_" + std::to_string(iteration-1) + "_iterations.txt"; 

//     std::ofstream resFile(path, std::ios::out);
    
//     for (int r = 0; r < residuo.size(); r++) 
//         resFile << residuo[r] << "\n";

//     resFile.close();

//     std::cout << "Text file \033[34m" << path << "\033[0;0m was successfully written." << std::endl;

//     modeling->reduce_boundary(model, IMAGE);

//     export_binary_float(gathers_folder + "IDLSKDM_dcal_" + std::to_string(nt) + "x" + std::to_string(nTraces) + "_iteration_" + std::to_string(iteration-1) + ".bin", dcal, nt*nTraces);
//     export_binary_float(images_folder + "IDLSKDM_model_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "_iteration_" + std::to_string(iteration-1) + ".bin", IMAGE, modeling->nPoints);
// }

// __global__ void forward_kernel(float * S, float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
// {
//     int index = blockIdx.x * blockDim.x + threadIdx.x;

//     int i = (int)(index % nzz);
//     int j = (int)(index / nzz);

//     if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
//     {
//         float T = Ts[index] + Tr[index]; 
//         int tId = __float2int_rd(T / dt);
        
//         if (tId < nt) 
//         {
//             float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
//             float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

//             float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
//             float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
//             float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

//             float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
//             float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

//             float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
//             float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
//             float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

//             float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
//             float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
//             float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
//             float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

//             float nx_norm = 0.0f, nz_norm = -1.0f; 
//             float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
//             float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

//             float a = d2Ts_dx2  + d2Tr_dx2;
//             float b = d2Ts_dxdz + d2Tr_dxdz;
//             float c = d2Ts_dz2  + d2Tr_dz2;

//             float detH = a*c - b*b;

//             float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

//             float R_s = max(Ts[index] / S[index], EPS);
//             float R_r = max(Tr[index] / S[index], EPS);

//             float G = 1.0f / sqrt(R_s * R_r);

//             float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

//             float R = 1.0f + 0.2f*cos(theta);

//             float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

//             atomicAdd(&data[tId], weights * m[index]);                    
//         }            
//     }    
// }

// __global__ void adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * m, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz)
// {
//     int index = blockIdx.x * blockDim.x + threadIdx.x;

//     int i = (int)(index % nzz);
//     int j = (int)(index / nzz);

//     if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
//     {
//         float T = Ts[index] + Tr[index]; 
//         int tId = __float2int_rd(T / dt);

//         if (tId < nt) 
//         {
//             float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
//             float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

//             float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
//             float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
//             float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

//             float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
//             float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

//             float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
//             float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
//             float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

//             float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
//             float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
//             float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
//             float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

//             float nx_norm = 0.0f, nz_norm = -1.0f; 
//             float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
//             float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

//             float a = d2Ts_dx2  + d2Tr_dx2;
//             float b = d2Ts_dxdz + d2Tr_dxdz;
//             float c = d2Ts_dz2  + d2Tr_dz2;

//             float detH = a*c - b*b;

//             float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

//             float R_s = max(Ts[index] / S[index], EPS);
//             float R_r = max(Tr[index] / S[index], EPS);

//             float G = 1.0f / sqrt(R_s * R_r);

//             float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

//             float R = 1.0f + 0.2f*cos(theta);

//             float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            
                
//             atomicAdd(&m[index], weights * data[tId]);
//         }
//     }    
// }
