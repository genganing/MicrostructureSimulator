#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <string.h>

using namespace std;

class Settings
{
    public:
    //Total system size. width: Nx, higth: Ny
    int Nx = 0;
    int Ny = 0;

    double interface_thickness;
    double grid_size;
    double time_step;

    double initial_concentration;
    double interfacial_energy;
    double interfacial_mobility;
    double kinetic_anisotropy_coff;
    double static_anisotropy_coff;

    double diffusion_coff_solid;
    double diffusion_coff_liquid;

    string Pytorch_module_path;

    //MPI settings
    int numprocess = 1;
    int halo_layer = 1;


    int _Nx;
    int _Ny;
    Settings(){};
    Settings(string input_filename, int num_process, int process_rank);

    void Mpisetting(int num_process, int process_rank);
    void initializer(string input_filename, int process_rank);
    ~Settings() {};
};
#endif