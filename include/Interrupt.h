#ifndef INTERRUPT_H
#define INTERRUPT_H

#include <iostream>
#include "Vector2D.h"
#include "Settings.h"
#include "Phasefieldsolver.h"

using namespace std;

class Settings;
class Phasefieldsolver;

class Interruptsolver
{
public:
    Vector2D<double> term_interrupt;

    Interruptsolver() {};
    Interruptsolver(const Settings& localsettings);
    void Create_fields(const Settings& localsettings);
    void Term_interrupt(const Settings& localsettings, const Phasefieldsolver& localphasefield, int interrupt_coe, double intensity_coe);
    ~Interruptsolver() {};
};
#endif