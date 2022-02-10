#include <iostream>
#include "../../libtorch/include/torch/script.h"
#include <memory>
#include "../include/Settings.h"
#include <string.h>
#include "mpich/mpi.h"
#include "../include/MPIcommunication.h"
#include "../include/Phasefieldsolver.h"
#include "../include/Solutefieldsolver.h"
#include "../include/Orientationfieldsolver.h"
#include "../include/Fieldwriter.h"
#include "../include/Temperaturefieldsolver.h"
#include "../include/Fieldinitializer.h"
#include "../include/Boundaryconditions.h"

using namespace std;

int main(int argc, char **argv)
{
    int process_rank;
    int numprocs;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    Settings setting("../inputfile/inputfile",numprocs, process_rank);
    Phasefieldsolver PFsolver(setting);
    Solutefieldsolver SFsolver(setting);
    Interruptsolver IRsolver(setting);
    Orientationfieldsolver OFsolver(setting,0.0);
    Temperaturefieldsolver TFsolver(setting,917.0);//initial temperature
    Drivingforcesolver DFsolver(setting);
    Boundaryconditions Boundarysetting;
    MPIcommunication MPIcom;
    
    Fieldinitializer F_initializer;
    F_initializer.Nucleation_fixed_loc(PFsolver, SFsolver, OFsolver, setting, process_rank, 400, 400, 2, 0.0);//Coordinates, radius and preferred orientation of initial nucleation point

    MPIcom.Communication(SFsolver.c_phase_1, setting, process_rank);
    MPIcom.Communication(SFsolver.c_phase_2, setting, process_rank);
    MPIcom.Communication(PFsolver.phi, setting, process_rank);
    MPIcom.Communication(PFsolver.phi_old, setting, process_rank);
    MPIcom.Communication(PFsolver.dphi_dt, setting, process_rank);
    MPIcom.Communication(PFsolver.orientation, setting, process_rank);
    MPIcom.Communication(OFsolver.orientation_field, setting, process_rank);
    MPIcom.Communication(SFsolver.c_overall, setting, process_rank);
    MPIcom.Communication(SFsolver.c_overall_old, setting, process_rank);
    Boundarysetting.set(PFsolver, SFsolver, setting, 1, 1, 1, 1, process_rank);

    PFsolver.Update_fields(setting);
    SFsolver.Update_fields(setting);
    PFsolver.Update_interface_flag(setting); 

    for (int calc_step = 0; calc_step <= 30000; ++calc_step)
    {
        Boundarysetting.set(PFsolver, SFsolver, setting, 1, 1, 1, 1, process_rank);
        DFsolver.Calculate_c_phase(setting,PFsolver,SFsolver,TFsolver);
        MPIcom.Communication(SFsolver.c_phase_1, setting, process_rank);
        MPIcom.Communication(SFsolver.c_phase_2, setting, process_rank);
        Boundarysetting.set(PFsolver, SFsolver, setting, 1, 1, 1, 1, process_rank);
        
        DFsolver.Calculate_driving_force(setting,PFsolver,SFsolver,TFsolver);
        IRsolver.Term_interrupt(setting, PFsolver, 0, 1.0e-3);//Disturbance

        PFsolver.Solve(setting,OFsolver,DFsolver,IRsolver);
        PFsolver.Update_fields(setting);
        
        MPIcom.Communication(PFsolver.phi, setting, process_rank);
        MPIcom.Communication(PFsolver.phi_old, setting, process_rank);
        MPIcom.Communication(PFsolver.dphi_dt, setting, process_rank);
        MPIcom.Communication(PFsolver.orientation, setting, process_rank);
        Boundarysetting.set(PFsolver, SFsolver, setting, 1, 1, 1, 1, process_rank);

        SFsolver.solve(setting,PFsolver);
        SFsolver.Update_fields(setting);
        MPIcom.Communication(OFsolver.orientation_field, setting, process_rank);
        MPIcom.Communication(SFsolver.c_overall, setting, process_rank);
        MPIcom.Communication(SFsolver.c_overall_old, setting, process_rank);
        Boundarysetting.set(PFsolver, SFsolver, setting, 1, 1, 1, 1, process_rank);
        
        PFsolver.Update_fields(setting);
        SFsolver.Update_fields(setting);

        if (calc_step % 2000 == 0)
        {
            string output_filename = "../outputfile/" + to_string(calc_step);
            cout << output_filename << endl; 
            Fieldwriter Fwriter(output_filename);
            for (int i = 0; i!= numprocs; ++i)
            {
                if(i == process_rank)
                {
                    Fwriter.Write(PFsolver,SFsolver,OFsolver,DFsolver, TFsolver, setting,process_rank);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        if(process_rank == 0) cout << calc_step << endl; 
    }
    MPI_Finalize ();
    return 0;
}