# ifndef MODELING_HPP
# define MODELING_HPP

# include "../geometry/geometry.hpp"

class Modeling
{
private:

protected:

    std::string name;

public:

    virtual void set_name() = 0;

    void print_name();
};

# endif