# include "IDLSKDM.cuh"

void IDLSKDM::set_migration()
{
    domain = "Image Domain";
    migType = "IDLSKDM";
    m_samples = modeling->nPoints;

    output_path = seismic_folder + migType + "_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "_iteration_" + std::to_string(max_it) + ".bin";
}

void IDLSKDM::regularization()
{

    

}

void IDLSKDM::perform_adjoint()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, dt, nt, modeling->dx, modeling->dz, modeling->nxx, modeling->nzz, modeling->nb, aperture, CMP);
}

void IDLSKDM::perform_forward()
{
    image_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, dt, nt, modeling->dx, modeling->dz, modeling->nxx, modeling->nzz, modeling->nb, aperture, CMP);
}

void IDLSKDM::perform_adjoint_gradient()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_gradient, dt, nt, modeling->dx, modeling->dz, modeling->nxx, modeling->nzz, modeling->nb, aperture, CMP);
}

void IDLSKDM::perform_forward_direction()
{
    image_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_direction, dt, nt, modeling->dx, modeling->dz, modeling->nxx, modeling->nzz, modeling->nb, aperture, CMP);    
}
