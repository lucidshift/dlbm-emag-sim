#ifndef ELECTRODYNAMICS_H
#define ELECTRODYNAMICS_H

#include <vector>
//TODO: Add gpu vector types

#include "Fields.h"

namespace math
{

class Electrodynamics
{

public:

	static double getFieldEnergy(Fields::DensityVectorField3D&  bfield);
	static bool calcOrbitalCurrent(Fields::DensityVectorField3D&  field); //TODO: update types
	static double calcForce(Fields::DensityVectorField3D&  field); //TODO: update types
	static bool calcHField(Fields::DensityVectorField3D&  field); //TODO: update types
	static bool calcBField(Fields::DensityVectorField3D&  field); //TODO: update types
	static bool calcMu(Fields::DensityVectorField3D&  field); //TODO: update types

private:

	Electrodynamics() {};

};

}

#endif