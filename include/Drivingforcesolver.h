#ifndef DRIVINGFORCESOLVER_H
#define DRIVINGFORCESOLVER_H

#include <iostream>
#include "Vector2D.h"
#include "Settings.h"
#include "Phasefieldsolver.h"
#include "Solutefieldsolver.h"
#include "../../libtorch/include/torch/script.h"
#include <memory>
#include "string.h"

using namespace std;

class Settings;
class Orientationfieldsolver;
class Phasefieldsolver;
class Temperaturefieldsolver;
class Solutefieldsover;

class Drivingforcesolver
{
    public:
    Vector2D<double> drivingforce_field; //Orientation field
    torch::jit::script::Module module;
    string Pytorch_model_path = "./model_1026_1.pt";
    
    Drivingforcesolver() {};
    Drivingforcesolver(const Settings& localsettings);

    void Create_fields(const Settings& localsettings); //create field to store phi-related variales
    void Load_pytorch_module(const Settings& localsettings);
    
    void Calculate_driving_force(const Settings& localsettings, const Phasefieldsolver& local_phasefield, 
                                 const Solutefieldsolver& localsolutefield, const Temperaturefieldsolver& local_temp);
    void Calculate_c_phase(const Settings& localsettings, const Phasefieldsolver& local_phasefield, 
                           Solutefieldsolver& localsolutefield, const Temperaturefieldsolver& local_temp);
    ~Drivingforcesolver() {};
};
#endif