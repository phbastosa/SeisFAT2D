# include "../src/inversion/tomography_iso.hpp"
# include "../src/inversion/tomography_vti.hpp"

int main(int argc, char **argv)
{
    std::vector<Inversion *> inversion = 
    {
        new Tomography_ISO(), 
        new Tomography_VTI() 
    }; 
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("inversion_type", file));

    inversion[type]->parameters = file;

    inversion[type]->set_parameters();
    inversion[type]->import_obsData();

    auto ti = std::chrono::system_clock::now();

    while (true)
    {
        inversion[type]->forward_modeling();
        inversion[type]->check_convergence();

        if (inversion[type]->converged) break; 

        inversion[type]->optimization();
        inversion[type]->model_update();
    }

    auto tf = std::chrono::system_clock::now();

    inversion[type]->export_results();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;

    return 0;
}
