# include "modeling/serial_aFSM.hpp"
# include "modeling/parallel_aFSM.cuh"

int main(int argc, char **argv)
{
    std::vector<Eikonal *> modeling = 
    {
        new Serial_aFSM(),
        new Parallel_aFSM()
    };

    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("modeling_type", file));    

    modeling[type]->parameters = file;

    modeling[type]->set_parameters();

    auto ti = std::chrono::system_clock::now();

    for (int shot = 0; shot < modeling[type]->geometry->nrel; shot++)
    {
        modeling[type]->srcId = shot;

        modeling[type]->show_information();

        modeling[type]->initialization();
        modeling[type]->forward_solver();

        modeling[type]->get_synthetic_data();

        modeling[type]->export_synthetic_data();
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;
    
    return 0;
}