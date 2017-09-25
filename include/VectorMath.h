#ifndef VECTORMATH_H
#define VECTORMATH_H

#include <vector>
//TODO: Add gpu vector types

#include "Fields.h"

namespace math
{

class VectorMath
{

public:

	static Fields::DensityVectorField3D curl(Fields::DensityVectorField3D& field1, Fields::DensityVectorField3D& field2);

	static bool scalarMultiply(scalar double, Fields::DensityVectorField3D& field1);

private:

	VectorMath() {};

};

}

#endif