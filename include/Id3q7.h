#ifndef ID3Q7_H
#define ID3Q7_H

#include "Fields.h"

using physics::Fields;

class Id3q7
{

public:

	Id3q7(int xDim, int yDim, int zDim);

	~Id3q7();

	virtual bool iterate() = 0;

	virtual bool loadSource(Fields::DensityField3D inputSource) = 0;	//TODO: Make return and input types vectors

	virtual Fields::DensityField3D getArray() = 0;	//TODO: Make return and input types vectors

	virtual Fields::DensityField2D getSlice(int zDimSlice) = 0;	//TODO: Make return and input types vectors

private:

};

#endif