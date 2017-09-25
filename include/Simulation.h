#ifndef SIMULATION_H
#define SIMULATION_H

#include "Fields.h"
#include "Kinematics.h"
#include "Electrodynamics.h"

class Simulation
{

public:
	
	Simulation(x, y, z, int threadCount=1);

	~Simulation();

	DensityField3D getSimField();

	void insertFerriocs(IFerroic);

	void iterate();

	void iterateFor(int iterations);

private:

	int xDim;

};

#endif