#include <math.h>
#include <cstring>

#include "d3q7.h"

d3q7::d3q7(int xDim, int yDim, int zDim) :
	xDim(xDim),
	yDim(yDim),
	zDim(zDim)
{

/*	rho.reserve(sizeof(double)*xDim*yDim*zDim);
	rhoVector.reserve(sizeof(double)*velCount*xDim*yDim*zDim);
	rhoDisplay.reserve(sizeof(double)*yDim*zDim);
	source.reserve(sizeof(double)*xDim*yDim*zDim);

	leftBuffer.reserve(sizeof(double)*yDim*zDim); 		
	rightBuffer.reserve(sizeof(double)*yDim*zDim);		
	topBuffer.reserve(sizeof(double)*xDim*zDim); 		
	bottomBuffer.reserve(sizeof(double)*xDim*zDim);	    
	frontBuffer.reserve(sizeof(double)*xDim*yDim); 	    
	backBuffer.reserve(sizeof(double)*xDim*yDim);	*/	

/*	double ** test;
	test = createArray2D(10,10);*/

	double * test;
	test = new double [10]; 

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){

			test[x] = (double) x*y;
		}
	}

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){

			printf("%5.5f\n", test[x]);
		}
	}	

}

d3q7::~d3q7(){};

bool d3q7::iterate()
{
	collision();
	stream();
	density();

	return true;
}

bool d3q7::loadSource(DensityField3D& inputSource)
{

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
			for(int z=0; z<zDim; z++){

				source[x][y][z] = inputSource[x][y][z];
			}
		}
	}

	return true;
}

d3q7::DensityField3D * d3q7::getArray()
{
	return &rho;
}

d3q7::DensityField2D * d3q7::getSlice(int zDimSlice)
{
	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){

			rhoDisplay[x][y] = rho[x][y][zDimSlice];
		}
	}

	return &rhoDisplay;
}

//Private members

void d3q7::collision()
{

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
			for(int z=0; z<zDim; z++){
				for(int i=0; i<velCount; i++){
					//Finew = Fi + (1/T * (Fi0 - Fi)) + source
					rhoVector[i][x][y][z] = rhoVector[i][x][y][z] + ((1./T) * ((velDist[i] * rho[x][y][z]) - rhoVector[i][x][y][z])) + source[x][y][z] * velDist[i];
				}
			}
		}
	}
}


void d3q7::stream()
{

	//Velocity 1
	for(int y=yDim-1; y>=0; y--){
		for(int z=0; z<zDim; z++){
			leftBuffer[1][y][z] = rhoVector[1][xDim-1][y][z];
		}
	}
	memmove(&rhoVector[1][1][0][0], &rhoVector[1][0][0][0], (sizeof(double)*(xDim-1)*(yDim)*(zDim)));	

	//Velocity 2
	for(int x=xDim-1; x>=0; x--){
		for(int z=0; z<zDim; z++){
			bottomBuffer[2][x][z] = rhoVector[2][x][yDim-1][z];
		}	
	}
	memmove(&rhoVector[2][0][1][0], &rhoVector[2][0][0][0], (sizeof(double)*(xDim)*(yDim-1)*(zDim)));

	//Velocity 3
	for(int y=0; y<yDim; y++){
		for(int z=0; z<zDim; z++){
			rightBuffer[3][y][z] = rhoVector[3][0][y][z];
		}
	}
	memmove(&rhoVector[3][0][0][0], &rhoVector[3][1][0][0], (sizeof(double)*(xDim-1)*(yDim)*(zDim)));

	//Velocity 4
	for(int x=0; x<xDim; x++){
		for(int z=0; z<zDim; z++){
			topBuffer[4][x][z] = rhoVector[4][x][0][z];
		}
	}
	memmove(&rhoVector[4][0][0][0], &rhoVector[4][0][1][0], (sizeof(double)*(xDim)*(yDim-1)*(zDim)));

	//Velocity 5
	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
			frontBuffer[5][x][y] = rhoVector[5][x][y][zDim-1];
		}
	}
	memmove(&rhoVector[5][0][0][1], &rhoVector[5][0][0][0], (sizeof(double)*(xDim)*(yDim)*(zDim-1)));

	//Velocity 6
	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
			backBuffer[6][x][y] = rhoVector[6][x][y][0];
		}
	}
	memmove(&rhoVector[6][0][0][0], &rhoVector[6][0][0][1], (sizeof(double)*(xDim)*(yDim)*(zDim-1)));	


	//Merge velocities from buffers
	for(int x=0; x<xDim; x++){
		for(int z=0; z<zDim; z++){
			rhoVector[2][x][0][z] = bottomBuffer[2][x][z];
			rhoVector[4][x][yDim-1][z] = topBuffer[4][x][z];
		}
	}

	for(int y=0; y<yDim; y++){
		for(int z=0; z<zDim; z++){
			rhoVector[1][0][y][z] = leftBuffer[1][y][z];
			rhoVector[3][xDim-1][y][z] = rightBuffer[3][y][z];
		}
	}

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
			rhoVector[5][x][y][0] = frontBuffer[5][x][y];
			rhoVector[6][x][y][zDim-1] = backBuffer[6][x][y];
		}
	}	
}


void d3q7::density()
{

    for (int x=0;x<xDim;x++){
      	for (int y=0;y<yDim;y++){
			for(int z=0; z<zDim; z++){

	      		rho[x][y][z] = 0;

				for(int i=0; i<velCount; i++){
					rho[x][y][z] += rhoVector[i][x][y][z];
				}
			}
		}
  	}
}

double ** d3q7::createArray2D(int xSize, int ySize)
{
	double ** field2D;
	field2D = static_cast<double**>(malloc(xSize * sizeof(double*)));

	for(int i=0; i<xSize; i++)
	{
		field2D[i] = static_cast<double*>(malloc(ySize * sizeof(double)));
	}

	return field2D;
}