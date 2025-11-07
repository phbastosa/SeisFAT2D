# include "ADLSKDM.cuh"

void ADLSKDM::set_migration()
{
    domain = "Angle Domain";
    migType = "ADLSKDM";
    m_samples = modeling->nz*nang*nCMP;
    d_samples = nt*modeling->max_spread*modeling->geometry->nsrc;
}

void ADLSKDM::perform_forward()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_model, modeling->nxx, modeling->nzz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADLSKDM::perform_adjoint()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(d_Ts, d_Tr, d_data, d_model, modeling->nxx, modeling->nzz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}
