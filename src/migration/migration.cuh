# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

# define EPS 1e-6f

class Migration
{
protected:

    int old_nx, old_nz, old_nPoints;
    int new_nx, new_nz, new_nPoints;

    float old_dx, old_dz;
    float new_dx, new_dz;

    int nfreq, cmpId, nCMP;
    int nt, nw, nfft, max_it, nang;
    int nBlocks, d_samples, m_samples; 

    float ds, dr, dt, da, dCMP;
    float minCMP, maxCMP;
    float fmax, max_angle;
    float max_offset, CMP; 
    float aperture;

    bool anisotropy;
    
    float * seismic = nullptr; 
    float * wavelet = nullptr;

    double * time_trace = nullptr;
    double * time_wavelet = nullptr;

    fftw_complex * freq_trace = nullptr;
    fftw_complex * freq_wavelet = nullptr;

    fftw_plan trace_forward_plan;
    fftw_plan trace_inverse_plan;
    fftw_plan wavelet_forward_plan;

    float * m1 = nullptr;
    float * d1 = nullptr;

    float * m2 = nullptr;
    float * d2 = nullptr;

    float * h_Ts = nullptr;
    float * h_Tr = nullptr;

    float * d_Ts = nullptr;
    float * d_Tr = nullptr;

    float * d_data = nullptr;
    float * h_data = nullptr;    

    float * h_model = nullptr;
    float * d_model = nullptr;

    std::string domain;
    std::string migType;
    std::string output_path;

    Modeling * modeling = nullptr;
    
    std::string current, keyword;
    std::string xpos, zpos, total;
    std::string current_operation;

    std::string input_data_folder;
    std::string input_data_prefix;
    
    std::string tables_folder;
    std::string seismic_folder;
    std::string residuo_folder;

    void set_interpolation();
    void set_anisotropy();
    void set_slowness();
    void set_wavelet();
    void set_gathers();

    void perform_cubic(float * input, float * output);

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
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;

public:
    
    std::string parameters;

    void set_parameters();

    //void dot_product_test();

    virtual void kirchhoff_depth_migration() = 0;

    virtual void export_outputs() = 0;
};

__global__ void image_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dt, int nt, 
                                            float old_dx, float old_dz, float new_dx, float new_dz, int old_nx, int old_nz, 
                                            int new_nxx, int new_nzz, int nb, float aperture, float cmp);


__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dt, int nt, 
                                            float old_dx, float old_dz, float new_dx, float new_dz, int old_nx, int old_nz, 
                                            int new_nxx, int new_nzz, int nb, float aperture, float cmp);

__global__ void angle_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);
__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);

__device__ float d_cubic1d(float P[4], float dx);
__device__ float d_cubic2d(float P[4][4], float dx, float dy);

# endif