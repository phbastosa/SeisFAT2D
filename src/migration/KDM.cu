# include "KDM.cuh"

void KDM::kirchhoff_depth_migration()
{
    set_src_travel_times();
    set_rec_travel_times();
    prepare_convolution();

    set_src_domain();

    initialization();

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
        
        import_binary_float(data_path, seismic, nt*modeling->max_spread);
     
        import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);

        cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
    
        set_current_src();

        show_information();

        spreadId = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {            
            cmp = 0.5f*(modeling->sx + modeling->geometry->xrec[modeling->recId]);
            cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;               

            import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);

            cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

            adjoint_convolution();

            cudaMemcpy(d_data, h_data, nt*sizeof(float), cudaMemcpyHostToDevice);
                
            perform_migration();

            ++spreadId;
        }
    }
}
