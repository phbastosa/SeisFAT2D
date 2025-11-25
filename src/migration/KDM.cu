# include "KDM.cuh"

void KDM::kirchhoff_depth_migration()
{
    set_src_travel_times();
    set_rec_travel_times();
    prepare_convolution();

    set_src_domain();

    cudaMemset(d_model, 0.0f, m_samples*sizeof(float));    

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
        import_binary_float(data_path, seismic, nt*modeling->geometry->nrec);
     
        import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
        cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
    
        set_current_src();
    
        current_operation = domain + " Kirchhoff Depth Migration";

        show_information();

        float sx = modeling->geometry->xsrc[modeling->srcId];
        
        for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
        {            
            float rx = modeling->geometry->xrec[modeling->recId];

            float offset = fabsf((sx - rx));

            if (offset < max_offset) 
            {
                CMP = 0.5f*(sx + rx);

                import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
                cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

                for (int tId = 0; tId < nt; tId++)
                    h_data[tId] = seismic[tId + modeling->recId*nt];

                adjoint_convolution();

                cudaMemcpy(d_data, h_data, nt * sizeof(float), cudaMemcpyHostToDevice);
                
                perform_adjoint();
            }
        }
    }
}

void KDM::export_outputs()
{
    cudaMemcpy(h_model, d_model, m_samples*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float(output_path, h_model, m_samples);
}
