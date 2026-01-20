# include "ADLSKDM.cuh"

void ADLSKDM::set_migration()
{
    domain = "Angle Domain";
    migType = "ADLSKDM";
    m_samples = modeling->nz*nang*nCMP;

    output_path = seismic_folder + migType + "_result_" + std::to_string(modeling->nz) + "x" + std::to_string(nCMP) + "x" + std::to_string(nang) + "_iteration_" + std::to_string(max_it) + ".bin";
}

void ADLSKDM::regularization()
{
    for (int cmpId = 0; cmpId < nCMP; cmpId++)
    {
        for (int zId = 0; zId < modeling->nz; zId++)
        {
            for (int aId = 0; aId < nang; aId++)
            {
                int targetId = zId + aId*modeling->nz + cmpId*nang*modeling->nz;
                int ang_forw = zId + (aId+1)*modeling->nz + cmpId*nang*modeling->nz;
                int ang_back = zId + (aId-1)*modeling->nz + cmpId*nang*modeling->nz;

                if (aId == 0)
                {
                    h_direction[targetId] += smooth*(h_model[ang_forw] - 2.0f*h_model[targetId] + h_model[ang_forw]) / (da*da);
                }
                else if (aId == nang-1)
                {
                    h_direction[targetId] += smooth*(h_model[ang_back] - 2.0f*h_model[targetId] + h_model[ang_back]) / (da*da);
                }
                else 
                {
                    h_direction[targetId] += smooth*(h_model[ang_forw] - 2.0f*h_model[targetId] + h_model[ang_back]) / (da*da);
                }
            }
        }
    }
}

void ADLSKDM::perform_forward()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, aperture, CMP, cmpId);
}

void ADLSKDM::perform_adjoint()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, aperture, CMP, cmpId);
}

void ADLSKDM::perform_adjoint_gradient()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_gradient, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, aperture, CMP, cmpId);
}

void ADLSKDM::perform_forward_direction()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_direction, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, aperture, CMP, cmpId);    
}
