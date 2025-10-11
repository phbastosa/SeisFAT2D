# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

# define NTHREADS 256

class Migration
{
protected:

    int nBlocks, max_it; 
    int nt, nang, nw, nfft;
    int nTraces, nCMP, *pSUM;

    float aperture, max_angle;
    float dt, ds, dr, da, fmax;

    bool anisotropy;

    float * seismic = nullptr; 
    float * wavelet = nullptr;

    float * trace_in = nullptr;
    float * trace_out = nullptr;

    double * time_trace = nullptr;
    double * time_result = nullptr;
    double * time_wavelet = nullptr;

    fftw_complex * freq_trace = nullptr;
    fftw_complex * freq_result = nullptr;
    fftw_complex * freq_wavelet = nullptr;

    fftw_plan trace_forward_plan;
    fftw_plan result_adjoint_plan;
    fftw_plan wavelet_forward_plan;

    float * ODCIG = nullptr;
    float * ADCIG = nullptr;
    float * IMAGE = nullptr;

    float * d_Tr = nullptr;

    float * h_angle = nullptr;
    float * h_trace = nullptr;
    float * h_image = nullptr;

    float * d_data = nullptr;    
    float * d_angle = nullptr;
    float * d_trace = nullptr;
    float * d_image = nullptr;

    Modeling * modeling = nullptr;

    std::string input_data_folder;
    std::string input_data_prefix;
    
    std::string tables_folder;
    std::string outputs_folder;

    void set_wavelet();
    void set_gathers();

    void show_information();    
    void set_rec_travel_times();
    void prepare_convolution();

    virtual void adjoint() = 0;
    virtual void forward() = 0;

public:
    
    std::string parameters;

    void set_parameters();
    
    virtual void image_building() = 0;
    virtual void export_outputs() = 0;
};

# endif