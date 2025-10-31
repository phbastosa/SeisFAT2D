# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

# define EPS 1e-6f

# define NTHREADS 256

class Migration
{
protected:

    int nBlocks, spreadId; 
    int nt, nang, nw, nfft;
    int nTraces, nCMP, *pSUM;

    float aperture, max_angle;
    float dt, ds, dr, da, fmax;

    float cmp;
    int cmpId, traceId;
    bool anisotropy, converged;

    float * seismic = nullptr; 
    float * wavelet = nullptr;

    double * time_trace = nullptr;
    double * time_wavelet = nullptr;

    fftw_complex * freq_trace = nullptr;
    fftw_complex * freq_wavelet = nullptr;

    fftw_plan trace_forward_plan;
    fftw_plan trace_inverse_plan;
    fftw_plan wavelet_forward_plan;

    float * h_Ts = nullptr;
    float * h_Tr = nullptr;

    float * d_Ts = nullptr;
    float * d_Tr = nullptr;

    float * d_data = nullptr;
    float * h_data = nullptr;    

    Modeling * modeling = nullptr;
    
    std::string current, keyword;
    std::string xpos, zpos, total;
    std::string current_operation;

    std::string input_data_folder;
    std::string input_data_prefix;
    
    std::string tables_folder;
    std::string images_folder;
    std::string gathers_folder;
    std::string residuo_folder;

    void set_wavelet();
    void set_gathers();

    void set_src_domain();
    void set_current_src();

    void set_rec_domain();
    void set_current_rec();

    void show_information();    
    void set_src_travel_times();
    void set_rec_travel_times();
    void prepare_convolution();

    void adjoint_convolution();
    void forward_convolution();

    virtual void set_migration() = 0;

public:
    
    std::string parameters;

    void set_parameters();

    virtual void kirchhoff_depth_migration() = 0;

    virtual void export_outputs() = 0;
};

# endif