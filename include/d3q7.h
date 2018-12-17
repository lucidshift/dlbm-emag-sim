#ifndef D3Q7_H
#define D3Q7_H

#include <stdint.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

static const int T = 1;
static const int velCount = 7;
static const std::vector<double> velDist = {0.25,0.125,0.125,0.125,0.125,0.125,0.125};

class d3q7
{

public:

	typedef double * DensityField2D;
	typedef double * DensityVectorField2D;
	typedef double * DensityField3D;
	typedef double * DensityVectorField3D;

	d3q7(int xDim, int yDim, int zDim);

	~d3q7();

	bool iterate();

	bool loadSource(DensityField3D inputSource);	//TODO: Make return and input types vectors

	DensityField3D getArray();	//TODO: Make return and input types vectors

	DensityField2D getSlice(int zDimSlice);	//TODO: Make return and input types vectors

private:

	int xDim;
	int yDim;
	int zDim;

	DensityField3D _rho;
	DensityVectorField3D _rhoVector;
	DensityField2D _rhoDisplay;
	DensityField3D _source;

	DensityVectorField2D _leftBuffer; 		//Velocity 1
	DensityVectorField2D _rightBuffer;		//Velocity 2
	DensityVectorField2D _topBuffer; 		//Velocity 3
	DensityVectorField2D _bottomBuffer;	    //Velocity 4
	DensityVectorField2D _frontBuffer; 	    //Velocity 5
	DensityVectorField2D _backBuffer;		//Velocity 6

	void collision();

	void stream();

	void density();
};

#endif