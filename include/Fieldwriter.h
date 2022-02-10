#ifndef FIELDWRITER_H
#define FIELDWRITER_H

#include "Phasefieldsolver.h"
#include "Solutefieldsolver.h"
#include "Orientationfieldsolver.h"
#include "Drivingforcesolver.h"
#include "Interrupt.h"
#include "Temperaturefieldsolver.h"
#include "Settings.h"
#include <string.h>
#include <fstream>

using namespace std;


class Fieldwriter
{
//no-flux: 0, symmetry: 1, periodic: 2
public:
    Fieldwriter(){};
    Fieldwriter(string _output_filename);

    void Write(const Phasefieldsolver& localphasefield, const Solutefieldsolver& localsolutefield, 
               const Orientationfieldsolver& local_orientation, const Drivingforcesolver& localdrivingforce, const Temperaturefieldsolver& localTemp,
               const Settings& localsetting,int process_rank);

    string output_filename;
    ofstream output_file;
    ~Fieldwriter() {};
};
#endif