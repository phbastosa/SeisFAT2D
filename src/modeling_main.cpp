# include "modeling/lfreq/isotropic/scalar.cuh"
# include "modeling/lfreq/isotropic/elastic.cuh"
# include "modeling/lfreq/isotropic/acoustic.cuh"

# include "modeling/hfreq/isotropic/serial_aFSM.hpp"
# include "modeling/hfreq/isotropic/parallel_aFSM.cuh"
# include "modeling/hfreq/isotropic/podvin_lecomte.cuh" 

int main(int argc, char **argv)
{
    std::vector<Modeling *> modeling = 
    {
        new Podvin_Lecomte(),
        new Serial_aFSM(),
        new Parallel_aFSM(),

        new Scalar(),
        new Acoustic(),
        new Elastic(),
    };

    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("modeling_type", file));

    modeling[type]->set_name();
    modeling[type]->print_name();

//     modeling->parameters = file;

//     modeling->set_parameters();

//     auto ti = std::chrono::system_clock::now();

//     for (int shot = 0; shot < modeling->geometry->nrel; shot++)
//     {
//         modeling->srcId = shot;

//         modeling->show_information();

//         modeling->initialization();
//         modeling->forward_solver();

//         modeling->get_synthetic_data();

//         modeling->export_synthetic_data();
//     }

//     auto tf = std::chrono::system_clock::now();

//     std::chrono::duration<double> elapsed_seconds = tf - ti;
//     std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;
    
    return 0;
}