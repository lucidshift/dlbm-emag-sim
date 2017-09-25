#ifndef D3Q7_H
#define D3Q7_H

#include <stdint.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "Fields.h"
#include "Id3q7.h"

using physics::Fields;

static const int T = 1;
static const int velCount = 7;
static const std::vector<double> velDist = {0.25,0.125,0.125,0.125,0.125,0.125,0.125};	

class d3q7 //: public Id3q7  //TODO: Fix this inheritance
{

public:

	d3q7(int xDim, int yDim, int zDim);

	~d3q7();

	bool iterate();

	bool loadSource(Fields::DensityField3D inputSource);	//TODO: Make return and input types vectors

	Fields::DensityField3D getArray();	//TODO: Make return and input types vectors

	Fields::DensityField2D getSlice(int zDimSlice);	//TODO: Make return and input types vectors

private:

	int xDim;
	int yDim;
	int zDim;

	Fields::DensityField3D _rho;
	Fields::DensityVectorField3D _rhoVector;
	Fields::DensityField2D _rhoDisplay;
	Fields::DensityField3D _source;

	Fields::DensityVectorField2D _leftBuffer; 		//Velocity 1
	Fields::DensityVectorField2D _rightBuffer;		//Velocity 2
	Fields::DensityVectorField2D _topBuffer; 		//Velocity 3
	Fields::DensityVectorField2D _bottomBuffer;	    //Velocity 4
	Fields::DensityVectorField2D _frontBuffer; 	    //Velocity 5
	Fields::DensityVectorField2D _backBuffer;		//Velocity 6

	void collision();

	void stream();

	void density();
};

#endif