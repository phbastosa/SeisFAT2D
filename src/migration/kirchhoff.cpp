# include "kirchhoff.hpp"

void Kirchhoff::set_parameters()
{
    modeling = new Eikonal_Iso();
    modeling->parameters = parameters;
    modeling->set_parameters();

    input_data_folder = catch_parameter("input_data_folder", parameters);
    input_data_prefix = catch_parameter("input_data_prefix", parameters);

    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));

    Tr = new float[modeling->nPoints]();
    Ts = new float[modeling->nPoints]();
    Im = new float[modeling->nPoints]();

    image = new float[modeling->nPoints]();

    seismic = new float[nt*modeling->max_spread]();
}

void Kirchhoff::image_building()
{
    get_receiver_traveltimes();

    run_cross_correlation();
}

void Kirchhoff::get_receiver_traveltimes()
{
    for (int recId = 0; recId < modeling->geometry->nrec; recId++)
    {
        modeling->recId = recId;

        initialization();

        modeling->forward_solver();

        modeling->recId = recId;

        export_receiver_traveltimes();
    }
}

void Kirchhoff::initialization()
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

            modeling->T[zi + xi*modeling->nzz] = modeling->S[zi + xi*modeling->nzz] * 
                sqrtf(powf((xi - modeling->nb)*modeling->dx - modeling->geometry->xrec[modeling->recId], 2.0f) + 
                      powf((zi - modeling->nb)*modeling->dz - modeling->geometry->zrec[modeling->recId], 2.0f));
        }
    }
}

void Kirchhoff::export_receiver_traveltimes()
{
    auto clear = system("clear");    
    modeling->reduce_boundary(modeling->T, Tr);
    export_binary_float("../outputs/travelTimeTables/traveltimes_receiver_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);    
}

void Kirchhoff::run_cross_correlation()
{
    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        read_seismic_data();
        
        modeling->show_information();

        std::cout << "\nKirchhoff depth migration \n\n";

        modeling->initialization();
        modeling->forward_solver();

        modeling->reduce_boundary(modeling->T, Ts);

        int sId = modeling->geometry->sInd[modeling->srcId];

        int spread = 0;

        for (modeling->recId = modeling->geometry->iRec[sId]; modeling->recId < modeling->geometry->fRec[sId]; modeling->recId++)
        {
            import_binary_float("../outputs/travelTimeTables/traveltimes_receiver_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);

            # pragma omp parallel for
            for (int index = 0; index < modeling->nPoints; index++)
            {
                int i = (int)(index % modeling->nz);
                int j = (int)(index / modeling->nz);

                Im[index] = Ts[index] + Tr[index] - modeling->synthetic_data[spread];

                int seisId = (int)(Im[index] / dt);

                if (seisId < nt) image[index] += seismic[seisId + spread*nt];
            }

            ++spread;
        }
    }
}

void Kirchhoff::read_seismic_data()
{
    std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

    import_binary_float(data_path, seismic, nt*modeling->geometry->spread[modeling->geometry->sInd[modeling->srcId]]);
}

void Kirchhoff::export_outputs()
{
    export_binary_float("../outputs/migratedImages/kirchhoff_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin", image, modeling->nPoints);
}