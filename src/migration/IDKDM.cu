# include "IDKDM.cuh"

void IDKDM::set_migration()
{
    domain = "Image Domain";
    migType = "IDKDM";
    m_samples = modeling->nPoints;
    d_samples = nt*modeling->max_spread*modeling->geometry->nsrc;
}

void IDKDM::perform_forward()
{
    image_domain_forward_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, modeling->nxx, modeling->nzz, nt, modeling->nb);
}

void IDKDM::perform_adjoint()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, modeling->nxx, modeling->nzz, nt, modeling->nb);
}
