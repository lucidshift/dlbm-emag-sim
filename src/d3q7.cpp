#include <math.h>
#include <cstring>

#include "d3q7.h"

d3q7::d3q7(int xDim, int yDim, int zDim) :
	xDim(xDim),
	yDim(yDim),
	zDim(zDim)
{
	_rho = new double [xDim * yDim * zDim];
	_rhoVector = new double [velCount * xDim * yDim * zDim];
	_rhoDisplay = new double [yDim * zDim];
	_source = new double [xDim * yDim * zDim];

	_leftBuffer = new double [velCount * yDim * zDim]; 		
	_rightBuffer = new double [velCount * yDim * zDim];		
	_topBuffer = new double [velCount * xDim * zDim]; 		
	_bottomBuffer = new double [velCount * xDim * zDim];	    
	_frontBuffer = new double [velCount * xDim * yDim]; 	    
	_backBuffer = new double [velCount * xDim * yDim];		
}

d3q7::~d3q7(){};

bool d3q7::iterate()
{
	collision();
	stream();
	density();

	return true;
}

bool d3q7::loadSource(DensityField3D inputSource)
{
	double (*source)[yDim][zDim] = (double(*)[yDim][zDim]) _source;
	double (*inputSourceTemp)[yDim][zDim] = (double(*)[yDim][zDim]) inputSource;

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
			for(int z=0; z<zDim; z++){

				source[x][y][z] = inputSourceTemp[x][y][z];
			}
		}
	}

	return true;
}

d3q7::DensityField3D d3q7::getArray()
{
	return _rho;
}

d3q7::DensityField2D d3q7::getSlice(int zDimSlice)
{
	double (*rho)[yDim][zDim] = (double(*)[yDim][zDim]) _rho;
	double (*rhoDisplay)[zDim] = (double(*)[zDim]) _rhoDisplay;

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){

			rhoDisplay[x][y] = rho[x][y][zDimSlice];
		}
	}

	return _rhoDisplay;
}

//Private members

void d3q7::collision()
{
	double (*rho)[yDim][zDim] = (double(*)[yDim][zDim]) _rho;
	double (*rhoVector)[xDim][yDim][zDim] = (double(*)[xDim][yDim][zDim]) _rhoVector;
	double (*source)[yDim][zDim] = (double(*)[yDim][zDim]) _source;

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
	double (*rhoVector)[xDim][yDim][zDim] = (double(*)[xDim][yDim][zDim]) _rhoVector;

	double (*leftBuffer)[yDim][zDim] = (double(*)[yDim][zDim]) _leftBuffer;
	double (*rightBuffer)[yDim][zDim] = (double(*)[yDim][zDim]) _rightBuffer;
	double (*topBuffer)[xDim][zDim] = (double(*)[xDim][zDim]) _topBuffer;
	double (*bottomBuffer)[xDim][zDim] = (double(*)[xDim][zDim]) _bottomBuffer;
	double (*frontBuffer)[xDim][yDim] = (double(*)[xDim][yDim]) _frontBuffer;
	double (*backBuffer)[xDim][yDim] = (double(*)[xDim][yDim]) _backBuffer;

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
	double (*rho)[yDim][zDim] = (double(*)[yDim][zDim]) _rho;
	double (*rhoVector)[xDim][yDim][zDim] = (double(*)[xDim][yDim][zDim]) _rhoVector;

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
