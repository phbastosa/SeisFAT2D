# include "eikonal.hpp"

void Eikonal::set_specifications()
{
    set_properties();    
    set_conditions();    

    synthetic_data = new float[max_spread]();
}

void Eikonal::set_boundaries()
{
    nb = 1;    
}

void Eikonal::initialization()
{
    sIdx = (int)(geometry->xsrc[geometry->sInd[srcId]] / dx) + nb;
    sIdz = (int)(geometry->zsrc[geometry->sInd[srcId]] / dz) + nb;

    for (int index = 0; index < matsize; index++) 
        T[index] = 1e6f;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int xi = sIdx + (j - 1);
            int zi = sIdz + (i - 1);

            T[zi + xi*nzz] = S[zi + xi*nzz] * sqrtf(powf((xi - nb)*dx - geometry->xsrc[geometry->sInd[srcId]], 2.0f) + 
                                                    powf((zi - nb)*dz - geometry->zsrc[geometry->sInd[srcId]], 2.0f));
        }
    }
}

void Eikonal::compute_seismogram()
{
    int spread = 0;

    for (recId = geometry->iRec[geometry->sInd[srcId]]; recId < geometry->fRec[geometry->sInd[srcId]]; recId++)
    {
        float x = geometry->xrec[recId];
        float z = geometry->zrec[recId];
 
        float x1 = floorf(x / dx)*dx;  
        float x2 = floorf(x / dx)*dx + dx;

        float z1 = floorf(z / dz)*dz;  
        float z2 = floorf(z / dz)*dz + dz;

        int i1 = (int)(z1 / dz) + nb;
        int i2 = (int)(z2 / dz) + nb;

        int j1 = (int)(x1 / dx) + nb;
        int j2 = (int)(x2 / dx) + nb;

        float q11 = T[i1 + j1*nzz];
        float q12 = T[i2 + j1*nzz];
        float q21 = T[i1 + j2*nzz];
        float q22 = T[i2 + j2*nzz];

        float p0 = 1.0 / ((x2 - x1) * (z2 - z1));

        float p1 = q11 * (x2 - x) * (z2 - z);
        float p2 = q21 * (x - x1) * (z2 - z);
        float p3 = q12 * (x2 - x) * (z - z1);
        float p4 = q22 * (x - x1) * (z - z1);

        synthetic_data[spread++] = p0*(p1 + p2 + p3 + p4);
    }
}

void Eikonal::export_synthetic_data()
{
    std::string data_file = data_folder + "eikonal_iso_nStations" + std::to_string(geometry->spread[srcId]) + "_shot_" + std::to_string(geometry->sInd[srcId]+1) + ".bin";
    export_binary_float(data_file, synthetic_data, geometry->spread[srcId]);    
}

