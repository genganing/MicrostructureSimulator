#ifndef MPICOMMUNICATION_H
#define MPICOMMUNICATION_H

#include "Phasefieldsolver.h"
#include "Orientationfieldsolver.h"
#include "Phasefieldsolver.h"
#include "Solutefieldsolver.h"
#include "mpich/mpi.h"
#include "Vector2D.h"

class MPIcommunication
{
    public:
    MPIcommunication(){};

    void Communication(Phasefieldsolver& local_phasefield, Solutefieldsolver& local_solutefield, Orientationfieldsolver& local_orientation, 
                       const Settings& local_setting, int process_rank);
    void Communication(Vector2D<double>& _arr, const Settings& localsettings, int process_rank);
    ~MPIcommunication() {};
};
#endif