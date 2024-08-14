# include "modeling/eikonal.cuh"

int main(int argc, char **argv)
{
    Modeling * modeling = new Eikonal();

    modeling->parameters = std::string(argv[1]);

    modeling->set_parameters();

    // set_runtime()

    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        modeling->shot_index = modeling->geometry->sInd[shot];

        modeling->initialization();
        modeling->forward_solver();


    }


    // get_runtime()

    return 0;
}