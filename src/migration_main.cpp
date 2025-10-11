# include "migration/KDM.cuh"
# include "migration/IDLSKDM.cuh"
# include "migration/ODLSKDM.cuh"
# include "migration/ADLSKDM.cuh"

int main(int argc, char **argv)
{
    std::vector<Migration *> migration = 
    {
        new KDM(),
        new IDLSKDM(),
        new ODLSKDM(),
        new ADLSKDM()
    };

    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("mig_type", file));    

    migration[type]->parameters = file;

    migration[type]->set_parameters();

    auto ti = std::chrono::system_clock::now();

    migration[type]->image_building();

    auto tf = std::chrono::system_clock::now();

    migration[type]->export_outputs();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;
    
    return 0;
}