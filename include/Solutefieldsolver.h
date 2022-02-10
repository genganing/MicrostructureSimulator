#ifndef SOLUTEFIELDSOLVER_H
#define SOLUTEFIELDSOLVER_H

#include <iostream>
#include "Vector2D.h"
#include "Settings.h"
#include "Phasefieldsolver.h"
#include "Drivingforcesolver.h"

using namespace std;

class Settings;
class Orientationfieldsolver;
class Drivingforcesolver;
class Phasefieldsolver;
class Temperaturefieldsolver;

class Solutefieldsolver
{
    public:
    Vector2D<double> c_overall; //solute field of next step
    Vector2D<double> c_overall_old; //solute field of previous step

    Vector2D<double> c_phase_1;  //liquid
    Vector2D<double> c_phase_2;  //solid
    Vector2D<double> dc_dt;

    Solutefieldsolver(){};
    Solutefieldsolver(const Settings& localsettings);

    void Create_fields(const Settings& localsettings);
    void Update_fields(const Settings& localsettings);
    void Solutefield_initialize(const Settings& localsettings); 
    void solve(const Settings& localsettings, const Phasefieldsolver& local_phasefield);
    ~Solutefieldsolver() {};
};
#endif