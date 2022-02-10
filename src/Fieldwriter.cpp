#include <iostream>
#include <fstream>
#include <string.h>
#include "../include/Phasefieldsolver.h"
#include "../include/Solutefieldsolver.h"
#include "../include/Orientationfieldsolver.h"
#include "../include/Settings.h"
#include "../include/Fieldwriter.h"
#include "../include/Interrupt.h"
#include "../include/Drivingforcesolver.h"
#include "../include/Temperaturefieldsolver.h"

using namespace std;

Fieldwriter::Fieldwriter(string _output_filename)
{
    this->output_filename = _output_filename + ".dat";
}

void Fieldwriter::Write(const Phasefieldsolver& localphasefield, const Solutefieldsolver& localsolutefield, 
                        const Orientationfieldsolver& local_orientation, const Drivingforcesolver& localdrivingforce, const Temperaturefieldsolver& localTemp,
                        const Settings& localsetting,int process_rank)
{
    if(process_rank == 0)
    {
        output_file.open(this->output_filename, ios::out);
        if(!output_file)
        {
            cout<<this->output_filename<<" file open failed!"<<endl;
        }
        output_file << "VARIABLES = " << "\"i\" " << "\"j\" " << "\"phi\" " << "\"C\" " << "\"interface_flag\" " << endl;
        output_file << "ZONE" << " I=" << localsetting.Nx << " J=" << localsetting.Ny << " F=POINT" << endl;
    }
    else
    {
        output_file.open(this->output_filename, ios::app);
    }

    int global_index_x = 0;
    int global_index_y = 0;

    for(int j = 1; j!=localsetting._Ny - 1;++j)
        for(int i = 1; i!=localsetting._Nx - 1;++i)
        {
            global_index_x = i;
            global_index_y = process_rank * (localsetting._Ny - 2) + j;
            
            output_file << global_index_x << " " << global_index_y << " " << localphasefield.phi[i][j] << " " << localsolutefield.c_overall[i][j] << " " << localphasefield.interface_flag[i][j] << endl;
        }
    output_file.close();
}