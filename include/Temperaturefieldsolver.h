#ifndef TEMPERATUREFIELDSOLVER_H
#define TEMPERATUREFIELDSOLVER_H

#include <iostream>
#include "Vector2D.h"
#include "Settings.h"

using namespace std;

class Temperaturefieldsolver
{
    public:
    Vector2D<double> temperature_field;
    
    Temperaturefieldsolver(){};
    Temperaturefieldsolver(const Settings& localsettings);
    Temperaturefieldsolver(const Settings& localsettings, double _temperature);

    void Create_fields(const Settings& localsettings); //create field to store phi-related variales
    void Initialize(const Settings& localsettings, double _temperature);
    ~Temperaturefieldsolver() {};
};
#endif