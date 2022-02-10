#ifndef PHASEFIELDSOLVER_H
#define PHASEFIELDSOLVER_H

#include <iostream>
#include "Vector2D.h"
#include "Settings.h"
#include "Drivingforcesolver.h"
#include "Interrupt.h"


using namespace std;

class Settings;
class Orientationfieldsolver;
class Drivingforcesolver;
class Interruptsolver;

class Phasefieldsolver
{
    public:
    Vector2D<double> phi; //phase field of next step
    Vector2D<double> phi_old; //phase field of current step
    Vector2D<double> dphi_dt;
    Vector2D<double> orientation;
    Vector2D<double> interface_flag;
    Vector2D<double> interfacial_energy;
    Vector2D<double> interfacial_mobility;
    Vector2D<double> phi_laplacian;
    
    Phasefieldsolver(){};
    Phasefieldsolver(const Settings& localsettings);

    void Create_fields(const Settings& localsettings); //create field to store phi-related variales
    void Initialize_fields(const Settings& localsettings);

    void Solve_derivatives(const Settings& localsettings, const Orientationfieldsolver& local_orientation);
    void Solve_equations(const Settings& localsettings, const Drivingforcesolver& local_drivingforce, const Interruptsolver& local_interrupt);
    void Update_fields(const Settings& localsettings); 
    void Update_interface_flag(const Settings& localsettings);
    void Solve(const Settings& localsettings, const Orientationfieldsolver& local_orientation,
               const Drivingforcesolver& local_drivingforce, const Interruptsolver& local_interrupt);
    void Communicate(const Settings& localsetting, int process_rank);
    ~Phasefieldsolver() {};
};
#endif