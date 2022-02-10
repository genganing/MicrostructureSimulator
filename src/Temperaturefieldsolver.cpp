#include "../include/Temperaturefieldsolver.h"
#include "../include/Settings.h"
#include <cmath>

using namespace std;

Temperaturefieldsolver::Temperaturefieldsolver(const Settings& localsettings)
{
    this->Create_fields(localsettings);
}

Temperaturefieldsolver::Temperaturefieldsolver(const Settings& localsettings, double _temperature)
{
    this->Create_fields(localsettings);
    this->Initialize(localsettings, _temperature); 
}

void Temperaturefieldsolver::Create_fields(const Settings& localsettings)
{
    this->temperature_field.Resize(localsettings._Nx,localsettings._Ny);
}

void Temperaturefieldsolver::Initialize(const Settings& localsettings, double _temperature)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            this->temperature_field[i][j] = _temperature;
        }
}