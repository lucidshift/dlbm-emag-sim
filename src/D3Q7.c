#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define XDIM 50
#define YDIM 50
#define ZDIM 50
#define velCount 7
#define u0 1

int xDim = XDIM;
int yDim = YDIM;
//int velCount = 5;
int iterationCount = 0;
int cpuCount=0;

double T = 1;
double velDist[velCount] = {0.25,0.125,0.125,0.125,0.125,0.125,0.125};

double A[velCount][XDIM][YDIM][ZDIM]; 	//D2Q5, X,Y  
double Source[XDIM][YDIM][ZDIM];	//Includes u0
double rho[XDIM][YDIM][ZDIM];
double B[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
double rhoDisplay[XDIM][YDIM];
double BDisplay[XDIM][YDIM][2];	//x,y (location), and 2D vector 

double leftBuffer[velCount][YDIM][ZDIM]; 	//Velocity 1
double rightBuffer[velCount][YDIM][ZDIM];	//Velocity 2
double topBuffer[velCount][XDIM][ZDIM]; 	//Velocity 3
double bottomBuffer[velCount][XDIM][ZDIM];	//Velocity 4
double frontBuffer[velCount][XDIM][YDIM]; 	//Velocity 5
double backBuffer[velCount][XDIM][YDIM];	//Velocity 6

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double slice[XDIM];
double sliceTh[XDIM];
double rhoTh[XDIM][YDIM];
int sliceLoc = 33;
double k = 1;
int rhoThreq=0;
int sliceReq=0;
int zSlice = 25;

struct current {double x; double y; double z;} J[XDIM][YDIM][ZDIM]; //Current vector at a given location
double JMagnitude = 1 * u0;


//Function Declaration
void Collision();
void Initialize();
void Stream();
void Curl();
void Iterate();
void GetSlice();
void GenerateSourcePath();


void Initialize(){

	/*
	for(int z=0; z<ZDIM; z++){
		Source[17][17][z] = 1 * u0;
		Source[34][17][z] = -1 * u0;
		Source[17][34][z] = -1 * u0;	
		Source[34][34][z] = 1 * u0;	
	}
	*/

	GenerateSourcePath();
	for(int x=0; x<(XDIM); x++){
		for(int y=0; y<(YDIM); y++){
			for(int z=0; z<(ZDIM); z++){
												//Don't forget you flipped these (x <-> z)
				Source[x][y][z] = JMagnitude * J[x][y][z].x;


			}
		}
	}	


	for(int x=0; x<(XDIM); x++){
		for(int y=0; y<(YDIM); y++){
			for(int z=0; z<(ZDIM); z++){
				for(int i=0; i<velCount; i++){

					A[i][x][y][z] = 0; //velDist[i] + Source[i][x][y];

					leftBuffer[i][y][z] = 0; 	//Velocity 1
					rightBuffer[i][y][z] = 0;	//Velocity 2
					topBuffer[i][x][z] = 0; 	//Velocity 3
					bottomBuffer[i][x][z] = 0;	//Velocity 4
					frontBuffer[i][x][y] = 0; 	//Velocity 5
					backBuffer[i][x][y] = 0;	//Velocity 6					
				}
			}
		}
	}		

    for (int x=0;x<(XDIM);x++){
      	for (int y=0;y<(YDIM);y++){
			for(int z=0; z<(ZDIM); z++){
		      		rho[x][y][z] = 0;
			}
		}
  	}

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

}


void GenerateSourcePath(){

	register int cx = 25;
	register int cy = 25;
	register int radius = 15;
	int xLow = cx - radius;
	int xHigh = cx + radius;
	//int yLow = cy - radius;
	//int yHigh = cy + radius; 

	//Create the "wire" for current to flow through
	for(int x=xLow; x<xHigh; x++){
		int y = sqrt((radius*radius) - (x-cx)*(x-cx)) + cy;
		for(int z=0; z<ZDIM; z++){

			//2D, -90 degree rotation in the x-y plane
			J[x][y][z].x = -y;
			J[x][y][z].y = x;
			J[x][y][z].z = 0;

			//Bottom half of the circle
			y = -y;
			J[x][y][z].x = -y;
			J[x][y][z].y = x;
			J[x][y][z].z = 0;			
		}
	}
	

}


void Collision(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				for(int i=0; i<velCount; i++){
					//Finew = Fi + (1/T * (Fi0 - Fi)) + Source
					A[i][x][y][z] = A[i][x][y][z] + ((1./T) * ((velDist[i] * rho[x][y][z]) - A[i][x][y][z])) + Source[x][y][z] * velDist[i];
				}
			}
		}
	}
}


void Stream(){

	//Velocity 1
	for(int y=YDIM-1; y>=0; y--){
		for(int z=0; z<ZDIM; z++){
			leftBuffer[1][y][z] = A[1][XDIM-1][y][z];
		}
	}
	memmove(&A[1][1][0][0], &A[1][0][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));	

	//Velocity 2
	for(int x=XDIM-1; x>=0; x--){
		for(int z=0; z<ZDIM; z++){
			bottomBuffer[2][x][z] = A[2][x][YDIM-1][z];
		}	
	}
	memmove(&A[2][0][1][0], &A[2][0][0][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 3
	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			rightBuffer[3][y][z] = A[3][0][y][z];
		}
	}
	memmove(&A[3][0][0][0], &A[3][1][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));

	//Velocity 4
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			topBuffer[4][x][z] = A[4][x][0][z];
		}
	}
	memmove(&A[4][0][0][0], &A[4][0][1][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 5
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			frontBuffer[5][x][y] = A[5][x][y][0];
		}
	}
	memmove(&A[5][0][0][1], &A[5][0][0][0], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));

	//Velocity 6
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			backBuffer[5][x][y] = A[5][x][y][ZDIM-1];
		}
	}
	memmove(&A[6][0][0][0], &A[6][0][0][1], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));	


	//Merge velocities from buffers
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			A[2][x][0][z] = bottomBuffer[2][x][z];
			A[4][x][YDIM-1][z] = topBuffer[4][x][z];
		}
	}

	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			A[1][0][y][z] = leftBuffer[1][y][z];
			A[3][XDIM-1][y][z] = rightBuffer[3][y][z];
		}
	}

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			A[5][x][y][0] = frontBuffer[5][x][y];
			A[6][x][y][ZDIM-1] = backBuffer[6][x][y];
		}
	}	
}


void Density(){

    for (int x=0;x<XDIM;x++){
      	for (int y=0;y<YDIM;y++){
			for(int z=0; z<ZDIM; z++){
	      		rho[x][y][z] = 0;
				for(int i=0; i<velCount; i++){
					rho[x][y][z] += A[i][x][y][z];
				}
			}
		}
  	}
}

void Curl(){

	int dx;
	int dy;
	
    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){
			for(int z=0; z<(ZDIM-1); z++){
      			dy = rho[x][y][z] - rho[x][y-1][z];	//Az
      			dx = rho[x][y][z] - rho[x-1][y][z];	//Az

      			B[x][y][zSlice][0] = dy;
      			B[x][y][zSlice][1] = -dx;
      		}
		}
  	}	 
}

void DeltaTime(){

	struct timespec spec;
	clock_gettime(CLOCK_REALTIME, &spec);

	t.past2 = t.past1;
	t.past1 = t.now;
	t.now = spec.tv_nsec / 1000000.0;
	t.dt = t.now - t.past1;
	printf ("deltaTime = %2.5f ms\n", t.dt);

}

void GetSlice(){

	double distance;
	int xOff;
	int yOff;
 	int zOff;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
		rhoTh[x][y] = 0;
		}
	}

	for(int x=0; x<XDIM; x++){
		//for(int y=0; y<YDIM; y++){
		int y = sliceLoc;	
		int z = zSlice;	
			for(int xx=0; xx<XDIM; xx++){
				for(int yy=0; yy<YDIM; yy++){	
					for(int zz=0; zz<ZDIM; zz++){
						for(int i=-3; i<4; i++){
							for(int j=-3; j<4; j++){
								for(int k=-3; k<4; k++){

									xOff = i * XDIM;
									yOff = j * YDIM;
									zOff = k * ZDIM;

									if(x!=xx || y!=yy){
										distance = sqrt(((x-xx+xOff)*(x-xx+xOff)) + ((y-yy+yOff)*(y-yy+yOff)) + ((z-zz+zOff)*(z-zz+zOff)));
									}
																
									rhoTh[x][y] += (Source[xx][yy][zz] / distance) / (3.1415926 * (T - 0.5));  // (J/d) / (4pi * 0.25 * M) 
								}
							}
						}
					}
				}
			}
		//}
		printf("X=%d\n",x);		
	}
	/*
	for(int x=0; x<XDIM; x++){
		printf("sliceTh[%d] = %5.5f\n",x,sliceTh[x]);
	}	
	*/
}



void GetGraph(){

	if(rhoThreq && sliceReq){ //ANDed because otherwise it constantly calculates
		GetSlice();
		rhoThreq = 0;
	}

	if(sliceReq){
		sliceReq=0;
		for(int x=0; x<XDIM; x++){
			slice[x] = rho[x][sliceLoc][zSlice];			
			sliceTh[x] = k * rhoTh[x][sliceLoc];
		}
	}

}

void Visualize(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){

			rhoDisplay[x][y] = rho[x][y][zSlice];

			//Projection into 2D
			BDisplay[x][y][0] = B[x][y][zSlice][0];
			BDisplay[x][y][1] = B[x][y][zSlice][1];
			
		}
	}
}


void Iterate(){

	Density();	
	Collision();
	Stream();
	Curl();
	Visualize();
	
	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	//DeltaTime();
	
}

//-----Main Routines-----
int main(){
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=100;
  	int drawCount=0;


	DefineGraphNxN_R("XY Density",&rhoDisplay[0][0],&xDim,&yDim,NULL);  
	DefineGraphNxN_R("Theoretical Density",&rhoTh[0][0],&xDim,&yDim,NULL);  		
	DefineGraphNxN_RxR("B Field(X,Y)",&BDisplay[0][0][0],&xDim,&yDim,NULL);
	DefineGraphN_R("Rho(x)",&slice[0],&xDim,NULL);
	DefineGraphN_R("RhoTh(x)",&sliceTh[0],&xDim,&sliceReq);


	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineDouble("Tau", &T);
    DefineInt("Slice", &sliceLoc);
    DefineBool("GetSlice",&rhoThreq);
    DefineDouble("k", &k);
    DefineInt("Z plane", &zSlice);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");    
  	DefineBool("Quit",&done);
  	EndMenu();

	Initialize();

  	while (!done){
  		if(drawCount>repeat){
      		Events(1);
      		GetGraph();
      		DrawGraphs();		
      		drawCount = 0;
  		}
	    if(cont){
	    	Iterate();
	    }
	    if(reset){
	    	Initialize();
	    	reset=0;
	    }

	    drawCount++;
  }
}