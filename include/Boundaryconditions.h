#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "Phasefieldsolver.h"
#include "Solutefieldsolver.h"
#include "Settings.h"

using namespace std;
class Boundaryconditions
{
//no-flux: 0, symmetry: 1, periodic: 2
//right, left, top, bottom
public:
    Boundaryconditions(){};
    void set(Phasefieldsolver& localphasefield, Solutefieldsolver& localsolutefield, const Settings& localsetting,
            int flag_right, int flag_left, int flag_top, int flag_bottom, int process_rank);
    ~Boundaryconditions() {};
};
#endif