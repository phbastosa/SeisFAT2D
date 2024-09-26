# include "least_squares.hpp"

void Least_Squares::set_specifications()
{
    dx_tomo = std::stof(catch_parameter("dx_tomo", parameters));
    dz_tomo = std::stof(catch_parameter("dz_tomo", parameters));

    if ((dx_tomo < modeling->dx) || (dz_tomo < modeling->dz))
        throw std::invalid_argument("\033[31mDownsampling with smaller spacing than original!\033[0;0m\n");

    nx_tomo = (int)((modeling->nx-1) * modeling->dx / dx_tomo) + 1;    
    nz_tomo = (int)((modeling->nz-1) * modeling->dz / dz_tomo) + 1;    

    tk_order = std::stoi(catch_parameter("tk_order", parameters));
    tk_param = std::stof(catch_parameter("tk_param", parameters));

    n_model = nx_tomo * nz_tomo;

    inversion_method = "Least-Squares First-Arrival Tomography";

    ray_path_max_samples = 0;

    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        std::cout << shot << " " << modeling->geometry->sInd[shot] << " " << modeling->geometry->iRec[shot] << " " << modeling->geometry->fRec[shot] << std::endl;  

        for (int node = modeling->geometry->iRec[shot]; node < modeling->geometry->fRec[shot]; node++)
        {
            float dx = (modeling->geometry->xsrc[modeling->geometry->sInd[shot]] - modeling->geometry->xrec[node]) / modeling->dx;
            float dz = (modeling->geometry->zsrc[modeling->geometry->sInd[shot]] - modeling->geometry->zrec[node]) / modeling->dz;
            
            ray_path_max_samples += (size_t)(2.0f*sqrtf(dx*dx + dz*dz));
        }
    }

    iG.reserve(ray_path_max_samples);
    jG.reserve(ray_path_max_samples);
    vG.reserve(ray_path_max_samples);
}

void Least_Squares::apply_inversion_technique()
{
    int sIdx = (int)(modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] / dx_tomo);
    int sIdz = (int)(modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] / dz_tomo);

    int sId = sIdz + sIdx*nz_tomo; 

    float rayStep = 0.2f * modeling->dz;

    std::vector < int > ray_index; 

    for (int ray_id = modeling->geometry->iRec[modeling->srcId]; ray_id < modeling->geometry->fRec[modeling->srcId]; ray_id++)
    {
        float xi = modeling->geometry->xrec[ray_id];        
        float zi = modeling->geometry->zrec[ray_id];

        if ((modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] == zi) && 
            (modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] == xi))
            continue;        

        while (true)
        {
            int j = (int)(xi / modeling->dx) + modeling->nb;
            int i = (int)(zi / modeling->dz) + modeling->nb;

            float dTx = (modeling->eikonalT[i + (j+1)*modeling->nzz] - modeling->eikonalT[i + (j-1)*modeling->nzz]) / (2.0f*modeling->dx);    
            float dTz = (modeling->eikonalT[(i+1) + j*modeling->nzz] - modeling->eikonalT[(i-1) + j*modeling->nzz]) / (2.0f*modeling->dz);    

            float norm = sqrtf(dTx*dTx + dTz*dTz);

            xi -= rayStep*dTx / norm;   
            zi -= rayStep*dTz / norm;    

            int jm = (int)(xi / dx_tomo); 
            int im = (int)(zi / dz_tomo); 

            int index = im + jm*nz_tomo;
            
            ray_index.push_back(index);

            if (ray_index.back() == sId) break;
        }
   
        float final_distance = sqrtf(powf(zi - modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]], 2.0f) + 
                                     powf(xi - modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]], 2.0f));

        std::sort(ray_index.begin(), ray_index.end());

        int current_voxel_index = ray_index[0];
        float distance_per_voxel = rayStep;

        for (int index = 0; index < ray_index.size(); index++)
        {
            if (ray_index[index] == current_voxel_index)
            {
                distance_per_voxel += rayStep;
            }
            else
            {
                vG.emplace_back(distance_per_voxel);
                jG.emplace_back(current_voxel_index);
                iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);

                if (current_voxel_index == sId) vG.back() = final_distance;

                distance_per_voxel = rayStep;
                current_voxel_index = ray_index[index];    
            }
        }

        if (current_voxel_index == sId)
        {
            vG.emplace_back(final_distance);
            jG.emplace_back(current_voxel_index);
            iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);
        }
        else 
        {
            vG.emplace_back(distance_per_voxel);
            jG.emplace_back(current_voxel_index);
            iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);
        }

        std::vector < int >().swap(ray_index);
    }
}

void Least_Squares::optimization()
{
    std::cout<<"\nSolving linear system using Tikhonov regularization with order " + std::to_string(tk_order) + "\n\n";

    M = n_model;                                  
    N = ndata + n_model - tk_order;                    
    NNZ = vG.size() + (tk_order + 1) * (n_model - tk_order);

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();
    x = new float[M]();

    for (int index = 0; index < ndata; index++) 
        B[index] = dobs[index] - dcal[index];

    for (int index = 0; index < vG.size(); index++)
    {
        iA[index] = iG[index];
        jA[index] = jG[index];
        vA[index] = vG[index];
    }

    std::vector< int >().swap(iG);
    std::vector< int >().swap(jG);
    std::vector<float>().swap(vG);

    apply_regularization();
    solve_linear_system_lscg();
    slowness_variation_rescaling();

    delete[] B;
    delete[] iA;
    delete[] jA;
    delete[] vA;    
}

void Least_Squares::apply_regularization()
{
    int elements = tk_order + 1;
		
    int n = n_model - tk_order;
    int nnz = elements * n;	
    
    int * iL = new int[nnz]();
    int * jL = new int[nnz]();
    float * vL = new float[nnz]();

    if (tk_order <= 0)
    {
	for (int index = 0; index < nnz; index++)
	{
        iL[index] = index;
	    jL[index] = index;
	    vL[index] = 1.0f;
	}
    } 
    else
    {
        int * df = new int[elements]();	
        int * df1 = new int[elements + 1]();
        int * df2 = new int[elements + 1]();
        
        df[0] = -1; df[1] = 1;
        
        for (int index = 1; index < tk_order; index++)
        {
            for (int k = 0; k < elements; k++)
            {
                df2[k] = df[k];
                df1[k + 1] = df[k];

                df[k] = df1[k] - df2[k]; 
            }		 
        }
        
        for (int index = 0; index < n; index++)
        {
            for (int k = 0; k < elements; k++)
            {
                iL[elements*index + k] = index;	
                jL[elements*index + k] = index + k;
                vL[elements*index + k] = df[k];
            }	
        }

        delete[] df;
        delete[] df1;
        delete[] df2;
    }

    for (int index = NNZ - nnz; index < NNZ; index++) 
    {
        iA[index] = ndata + iL[index - (NNZ - nnz)];
        jA[index] = jL[index - (NNZ - nnz)];
        vA[index] = tk_param * vL[index - (NNZ - nnz)];        
    }

    delete[] iL;
    delete[] jL;
    delete[] vL;
}

void Least_Squares::solve_linear_system_lscg()
{
    float a, b, qTq, rTr, rd;
    int cg_max_iteration = 10;

    float * s = new float[N]();
    float * q = new float[N]();
    float * r = new float[M]();
    float * p = new float[M]();

    // s = d - G * x, where d = dobs - dcal and x = slowness variation
    for (int i = 0; i < N; i++) 
        s[i] = B[i]; 

    // r = G' * s    
    for (int i = 0; i < NNZ; i++) 
        r[jA[i]] += vA[i] * s[iA[i]];        

    // p = r and x = 0;
    for (int i = 0; i < M; i++) 
    {
        p[i] = r[i]; 
        x[i] = 0.0f;
    }

    // q = G * p
    for (int i = 0; i < NNZ; i++) 
        q[iA[i]] += vA[i] * p[jA[i]];        

    for (int i = 0; i < cg_max_iteration; i++)
    {
        qTq = 0.0f;
        for (int k = 0; k < N; k++)           // q inner product
            qTq += q[k] * q[k];               // qTq = q' * q

        rTr = 0.0f;
        for (int k = 0; k < M; k++)           // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r 

        a = rTr / qTq;                        // a = (r' * r) / (q' * q)                    

        for (int k = 0; k < M; k++)           // model atualization
            x[k] += a * p[k];                 // x = x + a * p

        for (int k = 0; k < N; k++)           // s atualization  
            s[k] -= a * q[k];                 // s = s - a * q 

        rd = 0.0f;
        for (int k = 0; k < M; k++)           // r inner product for division 
            rd += r[k] * r[k];                // rd = r' * r

        for (int k = 0; k < M; k++)           // Zeroing r 
            r[k] = 0.0f;                      // r = 0, for multiplication
        
        for (int k = 0; k < NNZ; k++)         // r atualization 
            r[jA[k]] += vA[k] * s[iA[k]];     // r = G' * s    

        rTr = 0.0f;                
        for (int k = 0; k < M; k++)           // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r

        b = rTr / rd;                         // b = (r' * r) / rd

        for (int k = 0; k < M; k++)          
            p[k] = r[k] + b * p[k];           // p = r + b * p 

        for (int k = 0; k < N; k++) 
            q[k] = 0.0f;                      // q = 0, for multiplication

        for (int k = 0; k < NNZ; k++) 
            q[iA[k]] += vA[k] * p[jA[k]];     // q = G * p   
    }
}

void Least_Squares::slowness_variation_rescaling()
{
    // for (int index = 0; index < modeling->nPoints; index++)
    // {
    //     int k = (int) (index / (modeling->nx*modeling->nz));        
    //     int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
    //     int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  
        
    //     if ((i >= (int)(0.5f*dz_tomo/modeling->dz)) && (i < modeling->nz - (int)(0.5f*dz_tomo/modeling->dz)) &&
    //         (j >= (int)(0.5f*dx_tomo/modeling->dx)) && (j < modeling->nx - (int)(0.5f*dx_tomo/modeling->dx)) &&
    //         (k >= (int)(0.5f*dy_tomo/modeling->dy)) && (k < modeling->ny - (int)(0.5f*dy_tomo/modeling->dy)))
    //     {
    //         float zp = (i - (int)(0.5f*dz_tomo/modeling->dz))*modeling->dz; 
    //         float xp = (j - (int)(0.5f*dx_tomo/modeling->dx))*modeling->dx; 
    //         float yp = (k - (int)(0.5f*dy_tomo/modeling->dy))*modeling->dy; 
            
    //         float x0 = floorf(xp/dx_tomo)*dx_tomo;
    //         float y0 = floorf(yp/dy_tomo)*dy_tomo;
    //         float z0 = floorf(zp/dz_tomo)*dz_tomo;

    //         float x1 = floorf(xp/dx_tomo)*dx_tomo + dx_tomo;
    //         float y1 = floorf(yp/dy_tomo)*dy_tomo + dy_tomo;
    //         float z1 = floorf(zp/dz_tomo)*dz_tomo + dz_tomo;

    //         int idz = (int)(zp/dz_tomo);
    //         int idx = (int)(xp/dx_tomo);
    //         int idy = (int)(yp/dy_tomo);

    //         int ind_m = (int)(idz + idx*nz_tomo + idy*nx_tomo*nz_tomo);

    //         float c000 = x[ind_m];                  
    //         float c001 = x[ind_m + 1];
    //         float c100 = x[ind_m + nz_tomo];
    //         float c101 = x[ind_m + 1 + nz_tomo];
    //         float c010 = x[ind_m + nx_tomo*nz_tomo];
    //         float c011 = x[ind_m + 1 + nx_tomo*nz_tomo];
    //         float c110 = x[ind_m + nz_tomo + nx_tomo*nz_tomo];
    //         float c111 = x[ind_m + 1 + nz_tomo + nx_tomo*nz_tomo];  

    //         float xd = (xp - x0) / (x1 - x0);
    //         float yd = (yp - y0) / (y1 - y0);
    //         float zd = (zp - z0) / (z1 - z0);

    //         float c00 = c000*(1 - xd) + c100*xd;    
    //         float c01 = c001*(1 - xd) + c101*xd;    
    //         float c10 = c010*(1 - xd) + c110*xd;    
    //         float c11 = c011*(1 - xd) + c111*xd;    

    //         float c0 = c00*(1 - yd) + c10*yd;
    //         float c1 = c01*(1 - yd) + c11*yd;

    //         float dm_ijk = (c0*(1 - zd) + c1*zd);

    //         dm[index] = dm_ijk;            
    //     } 
    // }

    // if (smooth)
    // {
    //     int aux_nx = modeling->nx + 2*smoother_samples;
    //     int aux_ny = modeling->ny + 2*smoother_samples;
    //     int aux_nz = modeling->nz + 2*smoother_samples;

    //     int aux_nPoints = aux_nx*aux_ny*aux_nz;

    //     float * dm_aux = new float[aux_nPoints]();
    //     float * dm_smooth = new float[aux_nPoints]();

    //     for (int index = 0; index < modeling->nPoints; index++)
    //     {
    //         int k = (int) (index / (modeling->nx*modeling->nz));        
    //         int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
    //         int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

    //         int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

    //         dm_aux[ind_filt] = dm[i + j*modeling->nz + k*modeling->nx*modeling->nz];
    //     }

    //     smooth_volume(dm_aux, dm_smooth, aux_nx, aux_ny, aux_nz);

    //     for (int index = 0; index < modeling->nPoints; index++)
    //     {
    //         int k = (int) (index / (modeling->nx*modeling->nz));        
    //         int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
    //         int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

    //         int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

    //         dm[i + j*modeling->nz + k*modeling->nx*modeling->nz] = dm_smooth[ind_filt];
    //     }
    
    //     delete[] dm_aux;
    //     delete[] dm_smooth;
    // }    
}
