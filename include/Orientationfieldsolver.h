#ifndef ORIENTATIONFIELDSOLVER_H
#define ORIENTATIONFIELDSOLVER_H

#include <iostream>
#include "Vector2D.h"
#include "Settings.h"
#include "MPIcommunication.h"

using namespace std;

class Settings;
class MPIcommunication;

class Orientationfieldsolver
{
    friend class MPIcommunication;

public:
    Vector2D<double> orientation_field; //Orientation field

    Orientationfieldsolver() {};
    Orientationfieldsolver(const Settings& localsettings);
    Orientationfieldsolver(const Settings& localsettings, double _orientation);

    void Create_fields(const Settings& localsettings); //create field to store phi-related variales
    void Initialize(const Settings& localsettings, double _orientation);
    void solve(const Settings& localsettings, double _orientation);
    void communicate(const Settings& local_setting, int process_rank);
    ~Orientationfieldsolver() {};
};
#endif