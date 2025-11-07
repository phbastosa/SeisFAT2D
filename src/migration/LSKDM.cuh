# ifndef LSKDM_CUH
# define LSKDM_CUH

# include "migration.cuh"

class LSKDM : public Migration
{
private:

    bool converged;

    std::vector<float> residuo;

    float alpha, beta;
    
    float * gradient_p = nullptr;
    float * gradient_m = nullptr;
    
    float * direction = nullptr;

    void update_model();
    void initialization();
    void regularization();
    void compute_gradient();
    void compute_residuals();
    void compute_direction();
    void compute_stepLength();

protected:
    
    virtual void set_migration() = 0;
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;

public:

    void kirchhoff_depth_migration();

    void export_outputs();
};

# endif
