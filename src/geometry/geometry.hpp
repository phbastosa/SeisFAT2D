# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../ioFunctions/ioFunctions.hpp"

class Coord
{
public:

    int total;

    float * x = nullptr;
    float * z = nullptr;    
};

class Geometry
{
private:

    std::vector<float> linspace(float xi, float xf, int n);

protected:

    bool import_geometry;

    std::string rps_file;
    std::string sps_file;
    std::string xps_file;

    std::vector<std::string> splitted;

    void import_geometry();
    void export_geometry();
    
    void set_regular(Coord &obj);

public:

    Coord src;
    Coord rec;

    std::string file;

    virtual void set_geometry() = 0;     
};

# endif