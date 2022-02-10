#include <iostream>
#include <cmath>
#include "../include/Vector2D.h"
#include "../include/Settings.h"
#include "../include/Phasefieldsolver.h"
#include "../include/Orientationfieldsolver.h"
#include "../include/Solutefieldsolver.h"

using namespace std;

Solutefieldsolver::Solutefieldsolver(const Settings& localsettings)
{
    this->Create_fields(localsettings);
    this->Solutefield_initialize(localsettings);
}

void Solutefieldsolver::Create_fields(const Settings& localsettings)
{
    c_overall.Resize(localsettings._Nx, localsettings._Ny);
    c_overall_old.Resize(localsettings._Nx, localsettings._Ny);
    c_phase_1.Resize(localsettings._Nx, localsettings._Ny);
    c_phase_2.Resize(localsettings._Nx, localsettings._Ny);
    dc_dt.Resize(localsettings._Nx, localsettings._Ny);
}


void Solutefieldsolver::solve(const Settings& localsettings, const Phasefieldsolver& local_phasefield)
{
    for (int i = 1; i != localsettings._Nx - 1; ++i)
        for (int j = 1; j != localsettings._Ny - 1; ++j)
        {
            double phi_ipjp = (local_phasefield.phi[i + 1][j + 1] + local_phasefield.phi[i][j + 1]
                + local_phasefield.phi[i + 1][j] + local_phasefield.phi[i][j]) / 4.0;//dispersed
            double phi_ipjm = (local_phasefield.phi[i + 1][j] + local_phasefield.phi[i][j]
                + local_phasefield.phi[i + 1][j - 1] + local_phasefield.phi[i][j - 1]) / 4.0;
            double phi_imjp = (local_phasefield.phi[i][j + 1] + local_phasefield.phi[i - 1][j + 1]
                + local_phasefield.phi[i - 1][j] + local_phasefield.phi[i][j]) / 4.0;
            double phi_imjm = (local_phasefield.phi[i][j] + local_phasefield.phi[i - 1][j]
                + local_phasefield.phi[i - 1][j - 1] + local_phasefield.phi[i][j - 1]) / 4.0;

            //Current at the right-hand side-------------------------------------------------------------------
            double phi_here = (local_phasefield.phi[i + 1][j] + local_phasefield.phi[i][j]) / 2.0;
            double Cl_here = (c_phase_1[i + 1][j] + c_phase_1[i][j]) / 2.0;
            double Cs_here = (c_phase_2[i + 1][j] + c_phase_2[i][j]) / 2.0;
            double dphi_dt_here = (local_phasefield.dphi_dt[i + 1][j] + local_phasefield.dphi_dt[i][j]) / 2.0;
            double dphi_dx_here = (local_phasefield.phi[i + 1][j] - local_phasefield.phi[i][j]) / localsettings.grid_size;
            double dphi_dy_here = (phi_ipjp - phi_ipjm) / localsettings.grid_size;
            double MAG_here = sqrt(dphi_dx_here * dphi_dx_here + dphi_dy_here * dphi_dy_here);
            double dcphase1_dx_here = (c_phase_1[i + 1][j] - c_phase_1[i][j]) / localsettings.grid_size;
            double dcphase2_dx_here = (c_phase_2[i + 1][j] - c_phase_2[i][j]) / localsettings.grid_size;
            //Dl*(1-phi)*grad_cl
            double term1_R = localsettings.diffusion_coff_liquid * (1 - phi_here) * dcphase1_dx_here;
            //Ds*phi*grad_cs
            double term2_R = localsettings.diffusion_coff_solid * phi_here * dcphase2_dx_here;
            //calculate anti-trapping current term
            double term3_R = 0.0;
            if (MAG_here > 1e-10 / localsettings.grid_size)
            {
                term3_R = localsettings.interface_thickness / M_PI * sqrt(phi_here * (1 - phi_here))
                    * (Cl_here - Cs_here) * dphi_dt_here * dphi_dx_here / MAG_here;
            }
            double JR = term1_R + term2_R + term3_R;

            //Current at the left-hand side-------------------------------------------------------------------
            phi_here = (local_phasefield.phi[i][j] + local_phasefield.phi[i - 1][j]) / 2.0;
            Cl_here = (c_phase_1[i - 1][j] + c_phase_1[i][j]) / 2.0;
            Cs_here = (c_phase_2[i - 1][j] + c_phase_2[i][j]) / 2.0;
            dphi_dt_here = (local_phasefield.dphi_dt[i - 1][j] + local_phasefield.dphi_dt[i][j]) / 2.0;
            dphi_dx_here = (local_phasefield.phi[i][j] - local_phasefield.phi[i - 1][j]) / localsettings.grid_size;
            dphi_dy_here = (phi_imjp - phi_imjm) / localsettings.grid_size;
            MAG_here = sqrt(dphi_dx_here * dphi_dx_here + dphi_dy_here * dphi_dy_here);
            dcphase1_dx_here = (c_phase_1[i][j] - c_phase_1[i - 1][j]) / localsettings.grid_size;
            dcphase2_dx_here = (c_phase_2[i][j] - c_phase_2[i - 1][j]) / localsettings.grid_size;
            //Dl*(1-phi)*grad_cl
            double term1_L = localsettings.diffusion_coff_liquid * (1 - phi_here) * dcphase1_dx_here;
            //Ds*phi*grad_cs
            double term2_L = localsettings.diffusion_coff_solid * phi_here * dcphase2_dx_here;
            //calculate anti-trapping current term
            double term3_L = 0.0;
            if (MAG_here > 1e-10 / localsettings.grid_size)
            {
                term3_L = localsettings.interface_thickness / M_PI * sqrt(phi_here * (1 - phi_here))
                    * (Cl_here - Cs_here) * dphi_dt_here * dphi_dx_here / MAG_here;
            }
            double JL = term1_L + term2_L + term3_L;

            //Current at the top side-------------------------------------------------------------------
            phi_here = (local_phasefield.phi[i][j] + local_phasefield.phi[i][j + 1]) / 2.0;
            Cl_here = (c_phase_1[i][j] + c_phase_1[i][j + 1]) / 2.0;
            Cs_here = (c_phase_2[i][j] + c_phase_2[i][j + 1]) / 2.0;
            dphi_dt_here = (local_phasefield.dphi_dt[i][j] + local_phasefield.dphi_dt[i][j + 1]) / 2.0;
            dphi_dy_here = (local_phasefield.phi[i][j + 1] - local_phasefield.phi[i][j]) / localsettings.grid_size;
            dphi_dx_here = (phi_ipjp - phi_imjp) / localsettings.grid_size;
            MAG_here = sqrt(dphi_dx_here * dphi_dx_here + dphi_dy_here * dphi_dy_here);
            double dcphase1_dy_here = (c_phase_1[i][j + 1] - c_phase_1[i][j]) / localsettings.grid_size;
            double dcphase2_dy_here = (c_phase_2[i][j + 1] - c_phase_2[i][j]) / localsettings.grid_size;
            //Dl*(1-phi)*grad_cl
            double term1_T = localsettings.diffusion_coff_liquid * (1 - phi_here) * dcphase1_dy_here;
            //Ds*phi*grad_cs
            double term2_T = localsettings.diffusion_coff_solid * phi_here * dcphase2_dy_here;
            //calculate anti-trapping current term
            double term3_T = 0.0;
            if (MAG_here > 1e-10 / localsettings.grid_size)
            {
                term3_T = localsettings.interface_thickness / M_PI * sqrt(phi_here * (1 - phi_here))
                    * (Cl_here - Cs_here) * dphi_dt_here * dphi_dy_here / MAG_here;
            }
            double JT = term1_T + term2_T + term3_T;

            //Current at the bottom side-------------------------------------------------------------------
            phi_here = (local_phasefield.phi[i][j] + local_phasefield.phi[i][j - 1]) / 2.0;
            Cl_here = (c_phase_1[i][j] + c_phase_1[i][j - 1]) / 2.0;
            Cs_here = (c_phase_2[i][j] + c_phase_2[i][j - 1]) / 2.0;
            dphi_dt_here = (local_phasefield.dphi_dt[i][j] + local_phasefield.dphi_dt[i][j - 1]) / 2.0;
            dphi_dy_here = (local_phasefield.phi[i][j] - local_phasefield.phi[i][j - 1]) / localsettings.grid_size;
            dphi_dx_here = (phi_ipjm - phi_imjm) / localsettings.grid_size;
            MAG_here = sqrt(dphi_dx_here * dphi_dx_here + dphi_dy_here * dphi_dy_here);
            dcphase1_dy_here = (c_phase_1[i][j] - c_phase_1[i][j - 1]) / localsettings.grid_size;
            dcphase2_dy_here = (c_phase_2[i][j] - c_phase_2[i][j - 1]) / localsettings.grid_size;
            //Dl*(1-phi)*grad_cl
            double term1_B = localsettings.diffusion_coff_liquid * (1 - phi_here) * dcphase1_dy_here;
            //Ds*phi*grad_cs
            double term2_B = localsettings.diffusion_coff_solid * phi_here * dcphase2_dy_here;
            //calculate anti-trapping current term
            double term3_B = 0.0;
            if (MAG_here > 1e-10 / localsettings.grid_size)
            {
                term3_B = localsettings.interface_thickness / M_PI * sqrt(phi_here * (1 - phi_here))
                    * (Cl_here - Cs_here) * dphi_dt_here * dphi_dy_here / MAG_here;
            }
            double JB = term1_B + term2_B + term3_B;
            dc_dt[i][j] = localsettings.time_step / localsettings.grid_size * (JR - JL + JT - JB);
            c_overall[i][j] = c_overall_old[i][j] + dc_dt[i][j];
            //For debug-------------------------------------------------------------------------
            if (c_overall[i][j] <= 0 || c_overall[i][j] >= 1)
            {
            }
            //---------------------------------------------------------------------------------------------            
        }
}

void Solutefieldsolver::Update_fields(const Settings& localsettings)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            c_overall_old[i][j] = c_overall[i][j];
        }
}

void Solutefieldsolver::Solutefield_initialize(const Settings& localsettings)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            this->c_overall[i][j] = localsettings.initial_concentration;
            this->c_overall_old[i][j] = localsettings.initial_concentration;
        }
}