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

	typedef std::vector< std::vector<double> > DensityField2D;
	typedef std::vector< std::vector< std::vector<double> > > DensityVectorField2D;
	typedef std::vector< std::vector< std::vector<double> > > DensityField3D;
	typedef std::vector< std::vector< std::vector< std::vector<double> > > > DensityVectorField3D;

	d3q7(int xDim, int yDim, int zDim);

	~d3q7();

	bool iterate();

	bool loadSource(DensityField3D& inputSource);

	DensityField3D * getArray();

	DensityField2D * getSlice(int zDimSlice);

private:

	int xDim;
	int yDim;
	int zDim;

	DensityField3D rho;
	DensityVectorField3D rhoVector;
	DensityField2D rhoDisplay;
	DensityField3D source;

	DensityVectorField2D leftBuffer; 		//Velocity 1
	DensityVectorField2D rightBuffer;		//Velocity 2
	DensityVectorField2D topBuffer; 		//Velocity 3
	DensityVectorField2D bottomBuffer;	    //Velocity 4
	DensityVectorField2D frontBuffer; 	    //Velocity 5
	DensityVectorField2D backBuffer;		//Velocity 6

	double * test;

	void collision();

	void stream();

	void density();

	double * createArray2D(double * a, int xSize, int ySize);

};

#endif