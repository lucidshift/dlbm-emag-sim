#include <math.h>

#include "d3q7.h"

#define velCount 7

int xDim;
int yDim;
int zDim;

double T = 1;
double velDist[velCount] = {0.25,0.125,0.125,0.125,0.125,0.125,0.125};

double * rho;
double * rhoVector;
double * rhoDisplay;
double * source;

double * leftBuffer; 		//Velocity 1
double * rightBuffer;		//Velocity 2
double * topBuffer; 		//Velocity 3
double * bottomBuffer;	    //Velocity 4
double * frontBuffer; 	    //Velocity 5
double * backBuffer;		//Velocity 6


void d3q7(int xDim, int yDim, int zDim)
{

	xDim = xDim;
    yDim = yDim;
	zDim = zDim;

	//Allocate memory based on dimensions
	double rhoTemp[xDim][yDim][zDim];
	double rhoVectorTemp[velCount][xDim][yDim][zDim];
	double rhoDisplayTemp[xDim][yDim];

	double leftBufferTemp[velCount][yDim][zDim]; 	//Velocity 1
	double rightBufferTemp[velCount][yDim][zDim];	//Velocity 2
	double topBufferTemp[velCount][xDim][zDim]; 	//Velocity 3
	double bottomBufferTemp[velCount][xDim][zDim];	//Velocity 4
	double frontBufferTemp[velCount][xDim][yDim]; 	//Velocity 5
	double backBufferTemp[velCount][xDim][yDim];	//Velocity 6

	//Assign pointer locations
	rho = &rhoTemp;
	rhoVector = &rhoVectorTemp;
	rhoDisplay= &rhoDisplayTemp;

	leftBuffer = &leftBufferTemp; 		
	rightBuffer = &rightBufferTemp;		
	topBuffer = &topBufferTemp; 		
	bottomBuffer = &bottomBufferTemp;	    
	frontBuffer = &frontBufferTemp; 	    
	backBuffer = &backBufferTemp;		

}

bool iterate()
{
	collision();
	stream();
	density();
}

bool loadSource(double* source)
{
	source = &source;
}

double* getArray()
{
	return &rho;
}

double* getSlice(int zDimSlice)
{
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){

			rhoDisplay[x][y] = rhoVector[x][y][zDimSlice];
		}
	}

	return &rhoDisplay;
}

//"Private" members

void collision()
{

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				for(int i=0; i<velCount; i++){
					//Finew = Fi + (1/T * (Fi0 - Fi)) + source
					rhoVector[i][x][y][z] = rhoVector[i][x][y][z] + ((1./T) * ((velDist[i] * rho[x][y][z]) - rhoVector[i][x][y][z])) + source[x][y][z] * velDist[i];
				}
			}
		}
	}
}


void stream()
{

	//Velocity 1
	for(int y=YDIM-1; y>=0; y--){
		for(int z=0; z<ZDIM; z++){
			leftBuffer[1][y][z] = rhoVector[1][XDIM-1][y][z];
		}
	}
	memmove(&rhoVector[1][1][0][0], &rhoVector[1][0][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));	

	//Velocity 2
	for(int x=XDIM-1; x>=0; x--){
		for(int z=0; z<ZDIM; z++){
			bottomBuffer[2][x][z] = rhoVector[2][x][YDIM-1][z];
		}	
	}
	memmove(&rhoVector[2][0][1][0], &rhoVector[2][0][0][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 3
	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			rightBuffer[3][y][z] = rhoVector[3][0][y][z];
		}
	}
	memmove(&rhoVector[3][0][0][0], &rhoVector[3][1][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));

	//Velocity 4
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			topBuffer[4][x][z] = rhoVector[4][x][0][z];
		}
	}
	memmove(&rhoVector[4][0][0][0], &rhoVector[4][0][1][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 5
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			frontBuffer[5][x][y] = rhoVector[5][x][y][ZDIM-1];
		}
	}
	memmove(&rhoVector[5][0][0][1], &rhoVector[5][0][0][0], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));

	//Velocity 6
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			backBuffer[6][x][y] = rhoVector[6][x][y][0];
		}
	}
	memmove(&rhoVector[6][0][0][0], &rhoVector[6][0][0][1], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));	


	//Merge velocities from buffers
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			rhoVector[2][x][0][z] = bottomBuffer[2][x][z];
			rhoVector[4][x][YDIM-1][z] = topBuffer[4][x][z];
		}
	}

	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			rhoVector[1][0][y][z] = leftBuffer[1][y][z];
			rhoVector[3][XDIM-1][y][z] = rightBuffer[3][y][z];
		}
	}

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			rhoVector[5][x][y][0] = frontBuffer[5][x][y];
			rhoVector[6][x][y][ZDIM-1] = backBuffer[6][x][y];
		}
	}	
}


void density()
{

    for (int x=0;x<XDIM;x++){
      	for (int y=0;y<YDIM;y++){
			for(int z=0; z<ZDIM; z++){

	      		rho[x][y][z] = 0;

				for(int i=0; i<velCount; i++){
					rho[x][y][z] += rhoVector[i][x][y][z];
				}
			}
		}
  	}
}

