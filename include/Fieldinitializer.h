#ifndef FIELDINITIALIZER_H
#define FIELDINITIALIZER_H

#include "Orientationfieldsolver.h"
#include "Phasefieldsolver.h"
#include "Solutefieldsolver.h"
#include "Settings.h"

class Fieldinitializer
{
    public:
    Fieldinitializer(){};    
    Fieldinitializer(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, 
                     Orientationfieldsolver& local_orientationfield, Settings& local_setting);

    void Initialize(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, 
                    Orientationfieldsolver& local_orientationfield, Settings& local_setting, int nucleation_type);
    void Nucleation_fixed_loc(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, 
                    Orientationfieldsolver& local_orientationfield, Settings& local_setting, int process_rank, int x_loc, int y_loc, double radius, double orientation);
    ~Fieldinitializer() {};
};
#endif