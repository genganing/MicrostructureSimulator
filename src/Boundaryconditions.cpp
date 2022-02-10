#include <iostream>
#include <fstream>
#include <string>
#include "../include/Settings.h"
#include "../include/Boundaryconditions.h"
#include "../include/Phasefieldsolver.h"
#include "../include/Solutefieldsolver.h" 

using namespace std;

//no-flux: 0, symmetry: 1, periodic: 2
void Boundaryconditions::set(Phasefieldsolver& localphasefield, Solutefieldsolver& localsolutefield, const Settings& localsetting,
            int flag_right, int flag_left, int flag_top, int flag_bottom, int process_rank)
{
    if(process_rank == 0)
    {
        switch (flag_bottom)
        {
        case 0:
            {
                for (int i = 0; i != localsetting._Nx; ++i)
                {
                    localphasefield.phi[i][0] = localphasefield.phi[i][1];
                    localphasefield.phi_old[i][0] = localphasefield.phi_old[i][1];
                    localsolutefield.c_overall[i][0] = localsolutefield.c_overall[i][1];
                    localsolutefield.c_overall_old[i][0] = localsolutefield.c_overall_old[i][1];
                    localphasefield.dphi_dt[i][0] = localphasefield.dphi_dt[i][1];
                    localsolutefield.c_phase_1[i][0] = localsolutefield.c_phase_1[i][1];
                    localsolutefield.c_phase_2[i][0] = localsolutefield.c_phase_2[i][1];
                }
            }
            break;
        case 1:
            {
                for (int i = 0; i != localsetting._Nx; ++i)
                {
                    localphasefield.phi[i][0] = localphasefield.phi[i][2];
                    localphasefield.phi_old[i][0] = localphasefield.phi_old[i][2];
                    localphasefield.dphi_dt[i][0] = localphasefield.dphi_dt[i][2];
                    localsolutefield.c_overall[i][0] = localsolutefield.c_overall[i][2];
                    localsolutefield.c_overall_old[i][0] = localsolutefield.c_overall_old[i][2];
                    localsolutefield.c_phase_1[i][0] = localsolutefield.c_phase_1[i][2];
                    localsolutefield.c_phase_2[i][0] = localsolutefield.c_phase_2[i][2];
                }
            }
            break;
        case 2:
            {
                
            }
            break;
        case 3:
            {                
                //no treatment, because of the MPI data 
            }
            break;
        default:
            break;
        }        
    }

    if(process_rank == localsetting.numprocess - 1)
    {
        switch (flag_top)
        {
        case 0:
            {
                for (int i = 0; i != localsetting._Nx; ++i)
                {
                    localphasefield.phi[i][localsetting._Ny-1] = localphasefield.phi[i][localsetting._Ny-2];
                    localphasefield.phi_old[i][localsetting._Ny-1] = localphasefield.phi_old[i][localsetting._Ny-2];
                    localsolutefield.c_overall[i][localsetting._Ny-1] = localsolutefield.c_overall[i][localsetting._Ny-2];
                    localsolutefield.c_overall_old[i][localsetting._Ny-1] = localsolutefield.c_overall_old[i][localsetting._Ny-2];
                    localphasefield.dphi_dt[i][localsetting._Ny - 1] = localphasefield.dphi_dt[i][localsetting._Ny - 2];
                    localsolutefield.c_phase_1[i][localsetting._Ny - 1] = localsolutefield.c_phase_1[i][localsetting._Ny - 2];
                    localsolutefield.c_phase_2[i][localsetting._Ny - 1] = localsolutefield.c_phase_2[i][localsetting._Ny - 2];
                }
            }
            break;
        case 1:
            {
                for (int i = 0; i != localsetting._Nx; ++i)
                {
                    localphasefield.phi[i][localsetting._Ny-1] = localphasefield.phi[i][localsetting._Ny-3];
                    localphasefield.phi_old[i][localsetting._Ny - 1] = localphasefield.phi_old[i][localsetting._Ny - 3];
                    localphasefield.dphi_dt[i][localsetting._Ny - 1] = localphasefield.dphi_dt[i][localsetting._Ny - 3];
                    localsolutefield.c_overall[i][localsetting._Ny-1] = localsolutefield.c_overall[i][localsetting._Ny-3];
                    localsolutefield.c_overall_old[i][localsetting._Ny-1] = localsolutefield.c_overall_old[i][localsetting._Ny-3];
                    localsolutefield.c_phase_1[i][localsetting._Ny - 1] = localsolutefield.c_phase_1[i][localsetting._Ny - 3];
                    localsolutefield.c_phase_2[i][localsetting._Ny - 1] = localsolutefield.c_phase_2[i][localsetting._Ny - 3];
                }
            }
            break;
        case 2:
            {
                
            }
            break;
        case 3:
            {                
            //no treatment, because of the MPI data 
            }
            break;
        default:
            break;
        }
    }
    {
        switch (flag_left)
        {
        case 0:
            {
                for (int j = 0; j != localsetting._Ny; ++j)
                {
                    localphasefield.phi[0][j] = localphasefield.phi[1][j];
                    localphasefield.phi_old[0][j] = localphasefield.phi_old[1][j];
                    localphasefield.dphi_dt[0][j] = localphasefield.dphi_dt[1][j];
                    localsolutefield.c_overall[0][j] = localsolutefield.c_overall[1][j];
                    localsolutefield.c_overall_old[0][j] = localsolutefield.c_overall_old[1][j];
                    localphasefield.dphi_dt[0][j] = localphasefield.dphi_dt[1][j];
                    localsolutefield.c_phase_1[0][j] = localsolutefield.c_phase_1[1][j];
                    localsolutefield.c_phase_2[0][j] = localsolutefield.c_phase_2[1][j];
                }
            }
            break;
        case 1:
            {
                for (int j = 0; j != localsetting._Ny; ++j)
                {
                    localphasefield.phi[0][j] = localphasefield.phi[2][j];
                    localphasefield.phi_old[0][j] = localphasefield.phi_old[2][j];
                    localphasefield.dphi_dt[0][j] = localphasefield.dphi_dt[2][j];
                    localsolutefield.c_overall[0][j] = localsolutefield.c_overall[2][j];
                    localsolutefield.c_overall_old[0][j] = localsolutefield.c_overall_old[2][j];
                    localsolutefield.c_phase_1[0][j] = localsolutefield.c_phase_1[2][j];
                    localsolutefield.c_phase_2[0][j] = localsolutefield.c_phase_2[2][j];
                }
            }
            break;
        case 2:
            {
                
            }
            break;
        case 3:
            {
                
            }
        break;
        case 4:
            {
                for (int j = 0; j != localsetting._Ny; ++j)
                {
                    localphasefield.phi[0][j] = localphasefield.phi[localsetting._Nx-2][j];
                    localphasefield.phi_old[0][j] = localphasefield.phi_old[localsetting._Nx-2][j];
                    localsolutefield.c_overall[0][j] = localsolutefield.c_overall[localsetting._Nx-2][j];
                    localsolutefield.c_overall_old[0][j] = localsolutefield.c_overall_old[localsetting._Nx-2][j];
                }
            }
            break;
        default:
            break;
        }
        
        switch (flag_right)
        {
        case 0:
            {
                for (int j = 0; j != localsetting._Ny; ++j)
                {
                    localphasefield.phi[localsetting._Nx-1][j] = localphasefield.phi[localsetting._Nx-2][j];
                    localphasefield.phi_old[localsetting._Nx-1][j] = localphasefield.phi_old[localsetting._Nx-2][j];
                    localsolutefield.c_overall[localsetting._Nx-1][j] = localsolutefield.c_overall[localsetting._Nx-2][j];
                    localsolutefield.c_overall_old[localsetting._Nx-1][j] = localsolutefield.c_overall_old[localsetting._Nx-2][j];
                    localphasefield.dphi_dt[localsetting._Nx - 1][j] = localphasefield.dphi_dt[localsetting._Nx - 2][j];
                    localsolutefield.c_phase_1[localsetting._Nx - 1][j] = localsolutefield.c_phase_1[localsetting._Nx - 2][j];
                    localsolutefield.c_phase_2[localsetting._Nx - 1][j] = localsolutefield.c_phase_2[localsetting._Nx - 2][j];
                }
            }
            break;
        case 1:
            {
                for (int j = 0; j != localsetting._Ny; ++j)
                {
                    localphasefield.phi[localsetting._Nx - 1][j] = localphasefield.phi[localsetting._Nx - 3][j];
                    localphasefield.phi_old[localsetting._Nx - 1][j] = localphasefield.phi_old[localsetting._Nx - 3][j];
                    localphasefield.dphi_dt[localsetting._Nx - 1][j] = localphasefield.dphi_dt[localsetting._Nx - 3][j];
                    localsolutefield.c_overall[localsetting._Nx - 1][j] = localsolutefield.c_overall[localsetting._Nx - 3][j];
                    localsolutefield.c_overall_old[localsetting._Nx - 1][j] = localsolutefield.c_overall_old[localsetting._Nx - 3][j];
                    localsolutefield.c_phase_1[localsetting._Nx - 1][j] = localsolutefield.c_phase_1[localsetting._Nx - 3][j];
                    localsolutefield.c_phase_2[localsetting._Nx - 1][j] = localsolutefield.c_phase_2[localsetting._Nx - 3][j];
                }
            }
            break;
        case 2:
            {
                
            }
            break;
        case 3:
            {
                
            }
        break;
        case 4:
            {
                for (int j = 0; j != localsetting._Ny; ++j)
                {
                    localphasefield.phi[localsetting._Nx-1][j] = localphasefield.phi[1][j];
                    localphasefield.phi_old[localsetting._Nx-1][j] = localphasefield.phi_old[1][j];
                    localsolutefield.c_overall[localsetting._Nx-1][j] = localsolutefield.c_overall[1][j];
                    localsolutefield.c_overall_old[localsetting._Nx-1][j] = localsolutefield.c_overall_old[1][j];
                }
            }
            break;
        default:
            break;
        }
    }
}