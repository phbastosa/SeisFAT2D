# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

class Migration
{
private:

    int nt; 
    int nBlocks; 
    int nThreads;

    float dt; 
    float aperture;
    float max_offset;

    float * d_Tr = nullptr;

    float * f_image = nullptr;
    float * h_image = nullptr;
    float * d_image = nullptr;

    float * h_seismic = nullptr;
    float * d_seismic = nullptr;

    int nTraces;
    int nang, ncmp;
    float ds, dr, da;

    int * partial_cmp_sum = nullptr;

    float * h_ODCIG = nullptr;
    float * d_ODCIG = nullptr;

    float * h_ADCIG = nullptr;
    float * d_ADCIG = nullptr;

    float * ODCIG = nullptr;
    float * ADCIG = nullptr;

    std::string input_data_folder;
    std::string input_data_prefix;

    std::string output_image_folder;
    std::string output_table_folder;

    void show_information();
    void read_seismic_data();
    void set_common_gathers();
    void set_receiver_point();
    void get_receiver_eikonal();
    void run_cross_correlation();
    void export_receiver_eikonal();

protected:

    Modeling * modeling = nullptr;

    virtual void set_modeling_type() = 0;
    
public:
    
    std::string parameters;

    void set_parameters();
    void image_building();
    void export_outputs();
};

__global__ void cross_correlation(float * S, float * Ts, float * Tr, float * image, float * seismic, float * ODCIG, float * ADCIG, float aperture, float cmp, int spread, int nxx, int nzz, int nb, int nt, float dt, float dx, float dz);

# endif