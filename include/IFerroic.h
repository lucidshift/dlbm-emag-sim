#ifndef IFERROIC_H
#define IFERROIC_H

#include <vector>
//TODO: Add gpu vector types

#include "Fields.h"

namespace physics
{

class IFerroic
{

public:

	struct Properties {
		bool mobile;
	};

	virtual Properties getProperties() = 0;

	virtual bool getField(Fields::DensityVectorField3D& bfield) = 0;

private:

};

}

#endif