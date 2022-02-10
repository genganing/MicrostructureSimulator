#include "../include/MPIcommunication.h"
#include "../include/Settings.h"
#include "mpich/mpi.h"
#include <iostream>

using namespace std;

void MPIcommunication::Communication(Vector2D<double>& _arr, const Settings& local_setting, int process_rank)
{
    int _Nx = local_setting._Nx;
    int _Ny = local_setting._Ny;
    int num_process = local_setting.numprocess;
    MPI_Status status;

    double* array_send_up = new double[_Nx];
    double* array_send_down = new double[_Nx];
    double* array_recieve_up = new double[_Nx];
    double* array_recieve_down = new double[_Nx];

    for (int i = 1; i != _Nx-1;++i)
    {
        array_send_up[i] = _arr[i][_Ny - 2];
        array_send_down[i] = _arr[i][1];
    }

    if (num_process == 1)
    {
        for (int i = 1; i != _Nx-1;++i)
        {
            _arr[i][0] = _arr[i][_Ny - 2];
            _arr[i][_Ny - 1] = _arr[i][1];
        }
    }
    else
    {
        if (process_rank == 0)
        {
            MPI_Send(array_send_up, _Nx, MPI_DOUBLE, process_rank + 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down, _Nx, MPI_DOUBLE, num_process - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(array_recieve_up, _Nx, MPI_DOUBLE, process_rank + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down, _Nx, MPI_DOUBLE, num_process - 1, 1, MPI_COMM_WORLD, &status);
        }

        if (process_rank == num_process - 1)
        {
            MPI_Send(array_send_up,_Nx, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down, _Nx, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(array_recieve_up,_Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down, _Nx, MPI_DOUBLE, process_rank - 1, 1, MPI_COMM_WORLD, &status);
        }        

        if (process_rank != 0 && process_rank != num_process - 1)
        {
            MPI_Send(array_send_up, _Nx, MPI_DOUBLE, process_rank + 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down, _Nx, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(array_recieve_up, _Nx, MPI_DOUBLE, process_rank + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down, _Nx, MPI_DOUBLE, process_rank - 1, 1, MPI_COMM_WORLD, &status);            
        }

        for (int i = 1; i != _Nx-1; ++i)
        {
            _arr[i][_Ny - 1] = array_recieve_up[i];
            _arr[i][0] = array_recieve_down[i];
        }         
    }
    delete[] array_send_up;
    delete[] array_send_down;
    delete[] array_recieve_up;
    delete[] array_recieve_down;
}

void MPIcommunication::Communication(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, Orientationfieldsolver& local_orientation,
    const Settings& local_setting, int process_rank)
{
    int _Nx = local_setting._Nx;
    int _Ny = local_setting._Ny;
    int num_process = local_setting.numprocess;
    MPI_Status status;

    double* array_send_up_phi = new double[_Nx];
    double* array_send_up_conc = new double[_Nx];
    double* array_send_up_orientation = new double[_Nx];

    double* array_send_down_phi = new double[_Nx];
    double* array_send_down_conc = new double[_Nx];
    double* array_send_down_orientation = new double[_Nx];

    double* array_recieve_up_phi = new double[_Nx];
    double* array_recieve_up_conc = new double[_Nx];
    double* array_recieve_up_orientation = new double[_Nx];

    double* array_recieve_down_phi = new double[_Nx];
    double* array_recieve_down_conc = new double[_Nx];
    double* array_recieve_down_orientation = new double[_Nx];

    for (int i = 1; i != _Nx-1;++i)
    {
        array_send_up_phi[i] = local_phasefield.phi[i][_Ny - 2];
        array_send_up_conc[i] = local_solutefield.c_overall[i][_Ny - 2];
        array_send_up_orientation[i] = local_orientation.orientation_field[i][_Ny - 2];

        array_send_down_phi[i] = local_phasefield.phi[i][1];
        array_send_down_conc[i] = local_solutefield.c_overall[i][1];
        array_send_down_orientation[i] = local_orientation.orientation_field[i][1];
    }

    if (num_process == 1)
    {
        for (int i = 1; i != _Nx-1;++i)
        {
            local_phasefield.phi[i][0] = local_phasefield.phi[i][_Ny - 2];
            local_solutefield.c_overall[i][0] = local_solutefield.c_overall[i][_Ny - 2];
            local_orientation.orientation_field[i][0] = local_orientation.orientation_field[i][_Ny - 2];

            local_phasefield.phi[i][_Ny - 1] = local_phasefield.phi[i][1];
            local_solutefield.c_overall[i][_Ny - 1] = local_solutefield.c_overall[i][1];
            local_orientation.orientation_field[i][_Ny - 1] = local_orientation.orientation_field[i][1];
        }
    }
    else
    {
        if (process_rank == 0)
        {
            MPI_Send(array_send_up_phi, _Nx, MPI_DOUBLE, process_rank + 1, 3, MPI_COMM_WORLD);
            MPI_Send(array_send_up_conc, _Nx, MPI_DOUBLE, process_rank + 1, 4, MPI_COMM_WORLD);
            MPI_Send(array_send_up_orientation, _Nx, MPI_DOUBLE, process_rank + 1, 5, MPI_COMM_WORLD);

            MPI_Recv(array_recieve_up_phi, _Nx, MPI_DOUBLE, process_rank + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_up_conc, _Nx, MPI_DOUBLE, process_rank + 1, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_up_orientation, _Nx, MPI_DOUBLE, process_rank + 1, 2, MPI_COMM_WORLD, &status);

            MPI_Send(array_send_down_phi,_Nx, MPI_DOUBLE, num_process - 1, 0, MPI_COMM_WORLD);
            MPI_Send(array_send_down_conc,_Nx, MPI_DOUBLE, num_process - 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down_orientation,_Nx, MPI_DOUBLE, num_process - 1, 2, MPI_COMM_WORLD);

            MPI_Recv(array_recieve_down_phi,_Nx,MPI_DOUBLE,num_process - 1,3,MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_conc,_Nx,MPI_DOUBLE,num_process - 1,4,MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_orientation,_Nx,MPI_DOUBLE,num_process - 1,5,MPI_COMM_WORLD, &status);
        }
        else if (process_rank == num_process - 1)
        {
            MPI_Send(array_send_up_phi,_Nx, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            MPI_Send(array_send_up_conc,_Nx, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
            MPI_Send(array_send_up_orientation,_Nx, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
            
            MPI_Recv(array_recieve_up_phi,_Nx,MPI_DOUBLE,0,0,MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_up_conc,_Nx,MPI_DOUBLE,0,1,MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_up_orientation,_Nx,MPI_DOUBLE,0,2,MPI_COMM_WORLD, &status);

            MPI_Send(array_send_down_phi, _Nx, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD);
            MPI_Send(array_send_down_conc, _Nx, MPI_DOUBLE, process_rank - 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down_orientation, _Nx, MPI_DOUBLE, process_rank - 1, 2, MPI_COMM_WORLD);

            MPI_Recv(array_recieve_down_phi, _Nx, MPI_DOUBLE, process_rank - 1, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_conc, _Nx, MPI_DOUBLE, process_rank - 1, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_orientation, _Nx, MPI_DOUBLE, process_rank - 1, 5, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Send(array_send_up_phi, _Nx, MPI_DOUBLE, process_rank + 1, 3, MPI_COMM_WORLD);
            MPI_Send(array_send_up_conc, _Nx, MPI_DOUBLE, process_rank + 1, 4, MPI_COMM_WORLD);
            MPI_Send(array_send_up_orientation, _Nx, MPI_DOUBLE, process_rank + 1, 5, MPI_COMM_WORLD);

            MPI_Recv(array_recieve_up_phi, _Nx, MPI_DOUBLE, process_rank + 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_up_conc, _Nx, MPI_DOUBLE, process_rank + 1, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_up_orientation, _Nx, MPI_DOUBLE, process_rank + 1, 2, MPI_COMM_WORLD, &status);

            MPI_Send(array_send_down_phi, _Nx, MPI_DOUBLE, process_rank - 1, 0, MPI_COMM_WORLD);
            MPI_Send(array_send_down_conc, _Nx, MPI_DOUBLE, process_rank - 1, 1, MPI_COMM_WORLD);
            MPI_Send(array_send_down_orientation, _Nx, MPI_DOUBLE, process_rank - 1, 2, MPI_COMM_WORLD);

            MPI_Recv(array_recieve_down_phi, _Nx, MPI_DOUBLE, process_rank - 1, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_conc, _Nx, MPI_DOUBLE, process_rank - 1, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(array_recieve_down_orientation, _Nx, MPI_DOUBLE, process_rank - 1, 5, MPI_COMM_WORLD, &status);           
        }
        for (int i = 1; i != _Nx-1;++i)
        {
            local_phasefield.phi[i][_Ny - 1] = array_recieve_up_phi[i];
            local_solutefield.c_overall[i][_Ny - 1] = array_recieve_up_conc[i];
            local_orientation.orientation_field[i][_Ny - 1] = array_recieve_up_orientation[i];

            local_phasefield.phi[i][0] = array_recieve_down_phi[i];
            local_solutefield.c_overall[i][0] = array_recieve_down_conc[i];
            local_orientation.orientation_field[i][0] = array_recieve_down_orientation[i];
        }         
    }
    delete[] array_send_up_phi;
    delete[] array_send_up_conc;
    delete[] array_send_up_orientation;
    
    delete[] array_send_down_phi;
    delete[] array_send_down_conc;
    delete[] array_send_down_orientation;
    
    delete[] array_recieve_up_phi;
    delete[] array_recieve_up_conc;
    delete[] array_recieve_up_orientation;
    
    delete[] array_recieve_down_phi;
    delete[] array_recieve_down_conc;
    delete[] array_recieve_down_orientation;
}