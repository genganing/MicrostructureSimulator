#include <iostream>
#include <fstream>
#include <string>
#include "../include/Settings.h"

using namespace std;

Settings::Settings(string input_filename, int num_process, int process_rank)
{
    this->initializer(input_filename, process_rank);
    this->numprocess = num_process;
    this->_Nx = this->Nx + 2*this->halo_layer;
    this->_Ny = this->Ny/numprocess + 2*this->halo_layer;
}

void Settings::Mpisetting(int num_process, int process_rank) 
{

}

void Settings::initializer(string input_filename, int process_rank)
{
    ifstream input_file(input_filename,ios::in);
    if(!input_file)
    {
        cerr << "Can't open the input file." << endl;
        exit(-1);
    }

    input_file >> this->Nx;
    input_file >> this->Ny;
    input_file >> this->interface_thickness;
    input_file >> this->grid_size;
    input_file >> this->time_step;
    input_file >> this->initial_concentration;
    input_file >> this->interfacial_energy;
    input_file >> this->interfacial_mobility;
    input_file >> this->kinetic_anisotropy_coff;
    input_file >> this->static_anisotropy_coff;
    input_file >> this->diffusion_coff_solid;
    input_file >> this->diffusion_coff_liquid;
    input_file >> this->halo_layer;
    input_file >> this->Pytorch_module_path;

    if (process_rank == 0)
    {
        cout << "Nx: " << this->Nx << endl;
        cout << "Ny: " << this->Ny << endl;
        cout << "interface thickness: " << this->interface_thickness << endl;
        cout << "grid size: " << this->grid_size << endl;
        cout << "time step: " << this->time_step << endl;
        cout << "initial concentration: " << this->initial_concentration << endl;
        cout << "interfacial energy: " << this->interfacial_energy << endl;
        cout << "interfacial mobility: " << this->interfacial_mobility << endl;
        cout << "kinetic anisotropy mobility: " << this->kinetic_anisotropy_coff << endl;
        cout << "static anisotropy mobility: " << this->static_anisotropy_coff << endl;
        cout << "solid diffusion coefficient: " << this->diffusion_coff_solid << endl;
        cout << "liquid diffusion coefficient: "<< this->diffusion_coff_liquid << endl;
        cout << "halo layer thickness (dx): "<< this->halo_layer << endl;
        cout << "Pytorch module path: "<< this->Pytorch_module_path << endl;
    }
}