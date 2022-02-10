#include <iostream>
#include <cmath>
#include "../include/Vector2D.h"
#include "../include/Settings.h"
#include "../include/Phasefieldsolver.h"
#include "../include/Orientationfieldsolver.h"
#include "../include/Drivingforcesolver.h"
#include "../include/Interrupt.h"

using namespace std;

Phasefieldsolver::Phasefieldsolver(const Settings& localsettings)
{
    this->Create_fields(localsettings);
    this->Initialize_fields(localsettings);
}

void Phasefieldsolver::Create_fields(const Settings& localsettings)
{
    phi.Resize(localsettings._Nx,localsettings._Ny);
    phi_old.Resize(localsettings._Nx,localsettings._Ny);
    dphi_dt.Resize(localsettings._Nx,localsettings._Ny);
    orientation.Resize(localsettings._Nx, localsettings._Ny);
    interface_flag.Resize(localsettings._Nx, localsettings._Ny);
    interfacial_energy.Resize(localsettings._Nx, localsettings._Ny);
    interfacial_mobility.Resize(localsettings._Nx, localsettings._Ny);
    phi_laplacian.Resize(localsettings._Nx, localsettings._Ny);
}

void Phasefieldsolver::Initialize_fields(const Settings& localsettings)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            phi[i][j] = 0.0;
            phi_old[i][j] = 0.0;
        }
}

void Phasefieldsolver::Update_interface_flag(const Settings& localsettings)
{
    for (int i = 1; i != localsettings._Nx - 1; ++i)
        for (int j = 1; j != localsettings._Ny - 1; ++j)
        {
            if (this->phi_old[i][j] < 1.0e-2)
            {
                this->interface_flag[i][j] = 0;//liquid
            }
            else if (this->phi_old[i][j] > 1 - 1.0e-2)
            {
                this->interface_flag[i][j] = 1;//solid
            }
            else
            {
                this->interface_flag[i][j] = 0.5;//interface
            }
            for (int ii = -1; ii <= 1; ++ii)
            {
                for (int jj = -1 + abs(ii); jj <= 1 - abs(ii); ++jj)
                {
                    if (this->phi_old[i + ii][j + jj] > 1.0e-2 && this->phi_old[i + ii][j + jj] < 1 - 1.0e-2)
                    {
                        this->interface_flag[i][j] = 0.5;
                    }
                }
            }
        }
}


void Phasefieldsolver::Solve_derivatives(const Settings& localsettings, const Orientationfieldsolver& local_orientation)
{
    for (int i = 1; i != localsettings._Nx - 1; ++i)
        for (int j = 1; j != localsettings._Ny - 1; ++j)
        {
            double dphi_dx = (phi_old[i + 1][j] - phi_old[i - 1][j]) / 2.0 / localsettings.grid_size;
            double dphi_dy = (phi_old[i][j + 1] - phi_old[i][j - 1]) / 2.0 / localsettings.grid_size;          
     
            phi_laplacian[i][j] = (2.0 * (phi_old[i + 1][j] + phi_old[i - 1][j] + phi_old[i][j + 1] + phi_old[i][j - 1])
                + (phi_old[i + 1][j + 1] + phi_old[i - 1][j + 1] + phi_old[i + 1][j - 1] + phi_old[i - 1][j - 1]) - 12.0 * phi_old[i][j])
                / 4.0 / localsettings.grid_size / localsettings.grid_size;

            if (local_orientation.orientation_field[i][j] != 0)
            {
                orientation[i][j] = local_orientation.orientation_field[i][j];
            }

            if (phi_laplacian[i][j] != 0) 
            {
                if (orientation[i][j] == 0) 
                {
                    for (int ii = -1;ii <= 1;ii++) 
                    {
                        for (int jj = -1;jj <= 1;jj++) 
                        {
                            if (phi_old[i + ii][j + jj] > phi[i][j]) 
                            {
                                if (orientation[i + ii][j + jj] != 0)  orientation[i][j] = orientation[i + ii][j + jj];
                                if (local_orientation.orientation_field[i + ii][j + jj] != 0)  orientation[i][j] = local_orientation.orientation_field[i + ii][j + jj];
                            }
                        }
                    }
                }
            }    
                        
            interfacial_energy[i][j] = 0.0;
            interfacial_mobility[i][j] = 0.0;
            if (phi_laplacian[i][j] != 0.0) 
            {                      
                double beta = atan2(dphi_dy, dphi_dx)+ orientation[i][j];
                interfacial_energy[i][j] = localsettings.interfacial_energy * (1.0 - localsettings.static_anisotropy_coff * cos(4.0 * (beta)));
                interfacial_mobility[i][j] = localsettings.interfacial_mobility * (1.0 + localsettings.kinetic_anisotropy_coff * cos(4.0 * beta));
            }
        }
}

void Phasefieldsolver::Solve_equations(const Settings& localsettings, const Drivingforcesolver& local_drivingforce, const Interruptsolver& local_interrupt)
{
    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    for (int i = 1; i != localsettings._Nx -1; ++i)
        for (int j = 1; j != localsettings._Ny - 1; ++j)
        {    
            //Note, in the following step, phi_old should be 0~1 due to the square root
            if (phi_old[i][j] < 0.0) phi_old[i][j] = 0.0;
            if (phi_old[i][j] > 1.0) phi_old[i][j] = 1.0;
            double Vm = 1.0e-5;

            double term1 = interfacial_energy[i][j] * (phi_laplacian[i][j]
                + M_PI * M_PI /localsettings.interface_thickness / localsettings.interface_thickness * (phi_old[i][j] - 1.0 /2.0));

            double term2 = M_PI / localsettings.interface_thickness * sqrt(phi_old[i][j] * (1.0 - phi_old[i][j]))
                * local_drivingforce.drivingforce_field[i][j] / Vm;
                        
            dphi_dt[i][j] = interfacial_mobility[i][j] * (term1 + term2);

            phi[i][j] = dphi_dt[i][j] * localsettings.time_step + phi_old[i][j]+ local_interrupt.term_interrupt[i][j];//add term_interrupt
        }
}

void Phasefieldsolver::Update_fields(const Settings& localsettings)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        { 
            if(phi[i][j] < 0.0) phi[i][j] = 0.0; 
            if(phi[i][j] > 1.0) phi[i][j] = 1.0;           
            phi_old[i][j] = phi[i][j];
        }
    this->Update_interface_flag(localsettings); 
}

void Phasefieldsolver::Solve(const Settings& localsettings, const Orientationfieldsolver& local_orientation,
           const Drivingforcesolver& local_drivingforce, const Interruptsolver& local_interrupt)
{
    Solve_derivatives(localsettings, local_orientation);
    Solve_equations(localsettings, local_drivingforce, local_interrupt);
    //Update_fields(localsettings);
}