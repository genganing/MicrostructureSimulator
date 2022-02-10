#include "../include/Interrupt.h"
#include "../include/Vector2D.h"
#include "../include/Settings.h"
#include "../include/Phasefieldsolver.h"
#include <cmath>
#include <iostream>

using namespace std;

Interruptsolver::Interruptsolver(const Settings& localsettings)
{
    this->Create_fields(localsettings);
}

void Interruptsolver::Create_fields(const Settings& localsettings)
{
    term_interrupt.Resize(localsettings._Nx, localsettings._Ny);
}

void Interruptsolver::Term_interrupt(const Settings& localsettings, const Phasefieldsolver& localphasefield, int interrupt_coe, double intensity_coe)
{
    srand((unsigned)time(NULL));//Set random number
    if (interrupt_coe != 0) 
    {
        for (int i = 1; i != localsettings._Nx - 1; ++i)
            for (int j = 1; j != localsettings._Ny - 1; ++j)
            {
                if (localphasefield.phi_old[i][j] > 0.01 && localphasefield.phi_old[i][j] < 0.99)
                {                    
                    term_interrupt[i][j] = intensity_coe * (rand() % 10001 / 10000.0 - 0.5);
                }
                else
                {
                    term_interrupt[i][j] = 0;
                }
            }
    }
    else
    {
        for (int i = 1; i != localsettings.Nx - 1; ++i)
            for (int j = 1; j != localsettings._Ny - 1; ++j)
            {
                term_interrupt[i][j] = 0;
            }
    }
}