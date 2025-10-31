# ifndef LSKDM_CUH
# define LSKDM_CUH

# include "migration.cuh"

class LSKDM : public Migration
{
    virtual void set_migration() = 0;

    void image_building();

    virtual void export_outputs() = 0;
};

# endif
