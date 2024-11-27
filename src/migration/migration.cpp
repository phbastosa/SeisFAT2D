# include "migration.hpp"

void Migration::set_parameters()
{
    modeling = new Eikonal_Iso();
    modeling->parameters = parameters;
    modeling->set_parameters();

    aperture_x = std::stof(catch_parameter("mig_aperture_x", parameters));

    input_data_folder = catch_parameter("input_data_folder", parameters);
    input_data_prefix = catch_parameter("input_data_prefix", parameters);

    output_image_folder = catch_parameter("output_image_folder", parameters);
    output_table_folder = catch_parameter("output_table_folder", parameters);

    Tr = new float[modeling->nPoints]();
    Ts = new float[modeling->nPoints]();

    image = new float[modeling->nPoints]();

    seismic = new float[modeling->nt*modeling->max_spread]();

    set_specifications();
}

void Migration::read_seismic_data()
{
    std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

    import_binary_float(data_path, seismic, modeling->nt*modeling->geometry->spread[modeling->srcId]);
}

void Migration::image_building()
{
    get_receiver_traveltimes();

    run_cross_correlation();
}

void Migration::get_receiver_traveltimes()
{
    for (int recId = 0; recId < modeling->geometry->nrec; recId++)
    {
        modeling->recId = recId;

        initialization();
        
        show_information();

        modeling->forward_solver();

        modeling->recId = recId;

        export_receiver_traveltimes();
    }
}

void Migration::initialization()
{
    float rIdx = (int)(modeling->geometry->xrec[modeling->recId] / modeling->dx) + modeling->nb;
    float rIdz = (int)(modeling->geometry->zrec[modeling->recId] / modeling->dz) + modeling->nb;

    # pragma omp parallel for
    for (int index = 0; index < modeling->matsize; index++) 
        modeling->T[index] = 1e6f;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int xi = rIdx + (j - 1);
            int zi = rIdz + (i - 1);

            int index = zi + xi*modeling->nzz;

            modeling->T[index] = modeling->S[index] * 
                sqrtf(powf((xi - modeling->nb)*modeling->dx - modeling->geometry->xrec[modeling->recId], 2.0f) + 
                      powf((zi - modeling->nb)*modeling->dz - modeling->geometry->zrec[modeling->recId], 2.0f));
        }
    }
}

void Migration::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------------------\n";
    std::cout << "                                 \033[34mSeisFAT2D\033[0;0m\n";
    std::cout << "-------------------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (modeling->nz - 1)*modeling->dz << 
                                  ", x = " << (modeling->nx - 1)*modeling->dx << ") m\n\n";

    std::cout << "Running receiver " << modeling->recId + 1 << " of " << modeling->geometry->nrec << " in total\n\n";

    std::cout << "Current receiver position: (z = " << modeling->geometry->zrec[modeling->recId] << 
                                           ", x = " << modeling->geometry->xrec[modeling->recId] << ") m\n\n";

    std::cout << "Kirchhoff Depth Migration: computing receiver travel time matrices\n";
}

void Migration::export_receiver_traveltimes()
{
    modeling->reduce_boundary(modeling->T, Tr);
    export_binary_float(output_table_folder + "traveltimes_receiver_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);    
}

void Migration::export_outputs()
{
    export_binary_float(output_image_folder + "kirchhoff_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin", image, modeling->nPoints);
}
