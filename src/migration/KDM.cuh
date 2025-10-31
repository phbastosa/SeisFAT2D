# ifndef KDM_CUH
# define KDM_CUH

# include "migration.cuh"

class KDM : public Migration
{
protected:

    float * h_model = nullptr;
    float * d_model = nullptr;
    
    virtual void set_migration() = 0;
    virtual void initialization() = 0;
    virtual void perform_migration() = 0;
    
public:

    void kirchhoff_depth_migration();

    virtual void export_outputs() = 0;
};

# endif
