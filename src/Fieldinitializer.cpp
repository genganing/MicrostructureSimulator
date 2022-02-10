#include "../include/Fieldinitializer.h"
#include <cmath>

using namespace std;

Fieldinitializer::Fieldinitializer(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, 
                     Orientationfieldsolver& local_orientationfield, Settings& local_setting)
{
    
}

void Fieldinitializer::Initialize(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, 
                    Orientationfieldsolver& local_orientationfield, Settings& local_setting, int nucleation_type)
{

}

void Fieldinitializer::Nucleation_fixed_loc(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, 
                    Orientationfieldsolver& local_orientationfield, Settings& local_setting, int process_rank, int x_index, int y_index, double radius, double orientation)
{
    int local_index_x = x_index;
    int local_index_y = y_index - process_rank * (local_setting._Ny - 2); 

    for (int i = 0; i != local_setting._Nx; ++i)
        for (int j = 0; j != local_setting._Ny; ++j)
        {
            double distance = sqrt((local_index_x - i) * (local_index_x - i) + (local_index_y - j) * (local_index_y - j));
            if (distance <= radius)
            {
                local_phasefield.phi[i][j] = 1.0;
                local_orientationfield.orientation_field[i][j] = orientation;
            }
        }                   
}