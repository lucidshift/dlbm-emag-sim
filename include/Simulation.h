#ifndef SIMULATION_H
#define SIMULATION_H

#include "Fields.h"
#include "IFerroic.h"
#include "Electrodynamics.h"

class Simulation
{

public:
	
	Simulation(int x, int y, int z, int threadCount=1);

	~Simulation();

	Fields::DensityField3D getSimField();

	//void insertFerriocs(physics::IFerroic ferrioc);

	void iterate();

	void iterateFor(int iterations);

private:

	int xDim;

};

#endif