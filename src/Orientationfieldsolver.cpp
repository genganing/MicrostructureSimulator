#include <iostream>
#include "../include/Settings.h"
#include "../include/Orientationfieldsolver.h"
#include "mpich/mpi.h"

using namespace std;

Orientationfieldsolver::Orientationfieldsolver(const Settings& localsettings)
{
    this->Create_fields(localsettings);
}

Orientationfieldsolver::Orientationfieldsolver(const Settings& localsettings, double _orientation)
{
    this->Create_fields(localsettings);
    this->Initialize(localsettings,_orientation);
}

void Orientationfieldsolver::Create_fields(const Settings& localsettings)
{
    this->orientation_field.Resize(localsettings._Nx,localsettings._Ny);
}

void Orientationfieldsolver::Initialize(const Settings& localsettings, double _orientation)
{
    for (int i = 0; i != localsettings._Nx; ++i)
        for (int j = 0; j != localsettings._Ny; ++j)
        {
            this->orientation_field[i][j]=_orientation;
        }
}

void Orientationfieldsolver::solve(const Settings& localsettings, double _orientation)
{

}

void Orientationfieldsolver::communicate(const Settings& local_setting, int process_rank)
{
    int _Nx = local_setting._Nx;
    int _Ny = local_setting._Ny;
    int num_process = local_setting.numprocess;
    MPI_Status status;

    double *array_send_up_orientation = new double[_Nx];
    double *array_send_down_orientation = new double[_Nx];
    double *array_recieve_up_orientation = new double[_Nx];
    double *array_recieve_down_orientation = new double[_Nx];

    for(int i = 0; i!=_Nx;++i)
    {
        array_send_up_orientation[i] = orientation_field[i][_Ny - 2];
        array_send_down_orientation[i] = orientation_field[i][1];
    }

    if(num_process == 1)
    {
        for (int i = 0; i != _Nx;++i)
        {
            orientation_field[i][0] = orientation_field[i][_Ny - 2];
            orientation_field[i][_Ny - 1] = orientation_field[i][1];
        }
    }    
    else
    {     
        if(process_rank == 0)
        {
            MPI_Send(array_send_up_orientation, _Nx, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down_orientation, _Nx, MPI_DOUBLE, num_process - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(array_recieve_up_orientation, _Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_orientation, _Nx, MPI_DOUBLE, num_process - 1, 1, MPI_COMM_WORLD, &status);
        }

        if (process_rank == num_process - 1)
        {
            MPI_Send(array_send_up_orientation,_Nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down_orientation,_Nx, MPI_DOUBLE, process_rank-1, 0, MPI_COMM_WORLD);
            MPI_Recv(array_recieve_up_orientation,_Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_orientation,_Nx, MPI_DOUBLE, process_rank-1, 1, MPI_COMM_WORLD, &status);
        }
        
        if (process_rank != 0 && process_rank != num_process - 1)
        {
            MPI_Send(array_send_up_orientation, _Nx, MPI_DOUBLE, process_rank + 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down_orientation, _Nx, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(array_recieve_up_orientation, _Nx, MPI_DOUBLE, process_rank + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_orientation, _Nx, MPI_DOUBLE, process_rank - 1, 1, MPI_COMM_WORLD, &status);
        }

        for (int i = 0; i != _Nx;++i)
        {
            orientation_field[i][_Ny - 1] = array_recieve_up_orientation[i];
            orientation_field[i][0] = array_recieve_down_orientation[i];
        }         
    }  
    delete[] array_send_up_orientation;
    delete[] array_send_down_orientation;
    delete[] array_recieve_up_orientation;
    delete[] array_recieve_down_orientation; 
}