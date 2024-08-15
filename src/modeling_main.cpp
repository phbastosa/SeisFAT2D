# include "modeling/eikonal.cuh"

int main(int argc, char **argv)
{
    Modeling * modeling = new Eikonal();

    modeling->parameters = std::string(argv[1]);

    modeling->set_parameters();

    auto ti = std::chrono::system_clock::now();

    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        modeling->shot_index = modeling->geometry->sInd[shot];

        modeling->initialization();
        modeling->forward_solver();


    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;
    
    return 0;
}