# include "migration.hpp"

void Migration::set_parameters()
{
    modeling = new Eikonal_Iso();
    modeling->parameters = parameters;
    modeling->set_parameters();

    input_data_folder = catch_parameter("input_data_folder", parameters);
    input_data_prefix = catch_parameter("input_data_prefix", parameters);

    Tr = new float[modeling->nPoints]();
    Ts = new float[modeling->nPoints]();

    image = new float[modeling->nPoints]();

    seismic = new float[modeling->nt*modeling->max_spread]();

    set_components();
}

void Migration::read_seismic_data()
{
    std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

    import_binary_float(data_path, seismic, modeling->nt*modeling->geometry->spread[modeling->geometry->sInd[modeling->srcId]]);
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

            modeling->T[zi + xi*modeling->nzz] = modeling->S[zi + xi*modeling->nzz] * 
                sqrtf(powf((xi - modeling->nb)*modeling->dx - modeling->geometry->xrec[modeling->recId], 2.0f) + 
                      powf((zi - modeling->nb)*modeling->dz - modeling->geometry->zrec[modeling->recId], 2.0f));
        }
    }
}

void Migration::export_receiver_traveltimes()
{
    auto clear = system("clear");    
    modeling->reduce_boundary(modeling->T, Tr);
    export_binary_float("../outputs/travelTimeTables/traveltimes_receiver_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);    
}

void Migration::export_outputs()
{
    export_binary_float("../outputs/migratedImages/kirchhoff_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin", image, modeling->nPoints);
}