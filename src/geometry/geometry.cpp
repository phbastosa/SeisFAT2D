# include "geometry.hpp"

Geometry::Geometry(std::string parameters)
{    
    read_geometry = str2bool(catch_parameter("import_geometry", parameters));

    if (read_geometry)
    {
        std::vector<std::string> SPS;
        std::vector<std::string> RPS;
        std::vector<std::string> XPS;
        std::vector<std::string> AUX;

        import_text_file(sps_file, SPS); 
        import_text_file(rps_file, RPS); 
        import_text_file(xps_file, XPS); 

        nsrc = SPS.size();
        nrec = RPS.size();
        nrel = XPS.size();

        sInd = new int[nrel]();
        iRec = new int[nrel]();
        fRec = new int[nrel]();

        xsrc = new float[nsrc]();
        zsrc = new float[nsrc]();

        xrec = new float[nrec]();
        zrec = new float[nrec]();

        for (int i = 0; i < nrel; i++)
        {
            AUX = split(XPS[i], ',');

            sInd[i] = std::stoi(AUX[0]);
            iRec[i] = std::stoi(AUX[1]);
            fRec[i] = std::stoi(AUX[2]);
        }    

        std::vector<std::string>().swap(AUX);

        for (int i = 0; i < nsrc; i++)
        {
            AUX = split(SPS[i], ',');

            xsrc[i] = std::stof(AUX[0]);
            zsrc[i] = std::stof(AUX[1]);
        }    

        std::vector<std::string>().swap(AUX);

        for (int i = 0; i < nrec; i++)
        {
            AUX = split(RPS[i], ',');

            xrec[i] = std::stof(AUX[0]);
            zrec[i] = std::stof(AUX[1]);
        }    
    }
    else
    {
        sps_file = "../inputs/geometry/SeisFAT2D_SPS.txt";
        rps_file = "../inputs/geometry/SeisFAT2D_RPS.txt";
        xps_file = "../inputs/geometry/SeisFAT2D_XPS.txt";

        nsrc = std::stoi(catch_parameter("nsrc", parameters));
        nrec = std::stoi(catch_parameter("nrec", parameters));

        nrel = nsrc;

        sInd = new int[nrel]();
        iRec = new int[nrel]();
        fRec = new int[nrel]();

        xsrc = new float[nsrc]();
        zsrc = new float[nsrc]();

        xrec = new float[nrec]();
        zrec = new float[nrec]();

        float xsrc_beg = std::stof(catch_parameter("xsrc_beg", parameters));
        float xsrc_end = std::stof(catch_parameter("xsrc_end", parameters));

        float xrec_beg = std::stof(catch_parameter("xrec_beg", parameters));
        float xrec_end = std::stof(catch_parameter("xrec_end", parameters));

        float sz = std::stof(catch_parameter("zsrc", parameters));
        float rz = std::stof(catch_parameter("zrec", parameters));

        auto sx = linspace(xsrc_beg, xsrc_end, nsrc);
        auto rx = linspace(xsrc_beg, xsrc_end, nsrc);
        
        for (int i = 0; i < nrec; i++)
        {
            xrec[i] = rx[i];
            zrec[i] = rz;
        }

        for (int i = 0; i < nsrc; i++)
        {
            xsrc[i] = sx[i];
            zsrc[i] = sz;
        }
    
        for (int i = 0; i < nrel; i++)
        {
            sInd[i] = i;
            iRec[i] = 0;
            fRec[i] = nrec;
        }    
    
        std::ofstream SPS(sps_file, std::ios::out);        
        std::ofstream RPS(rps_file, std::ios::out);        
        std::ofstream XPS(xps_file, std::ios::out);        

        for (int i = 0; i < nsrc; i++)        
            SPS << xsrc[i] << ", " << zsrc[i] << "\n";    

        for (int i = 0; i < nrec; i++)        
            RPS << xrec[i] << ", " << zrec[i] << "\n";    
            
        for (int i = 0; i < nrel; i++)        
            XPS << sInd[i] << ", " << iRec[i] << ", " << fRec[i] << "\n";    
    }
}

std::vector<float> Geometry::linspace(float xi, float xf, int n)
{
    std::vector<float> linspaced;
    
    if (n == 0) return linspaced;
    if (n == 1)
    {
        linspaced.push_back(xi);
        return linspaced;
    } 

    linspaced.reserve(n);

    float delta = (xf - xi) / (n - 1);

    for (int i = 0; i < n; i++)
    {
        linspaced.emplace_back(xi + (float)(delta*i));
    }

    return linspaced;
}
