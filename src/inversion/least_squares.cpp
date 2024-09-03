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




    }
}

void Least_Squares::optimization()
{

    
}