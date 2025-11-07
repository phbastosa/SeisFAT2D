# include "ADKDM.cuh"

void ADKDM::set_migration()
{
    domain = "Angle Domain";
    migType = "ADKDM";
    m_samples = modeling->nz*nang*nCMP;
    d_samples = nt*modeling->max_spread*modeling->geometry->nsrc;
}

void ADKDM::perform_forward()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADKDM::perform_adjoint()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}
