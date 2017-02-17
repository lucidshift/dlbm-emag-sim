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
#define velCount 5

int xDim = XDIM;
int yDim = YDIM;
//int velCount = 5;
int iterationCount = 0;
int cpuCount=0;

double T = 1;
double velDist[5] = {0.5,0.125,0.125,0.125,0.125};

double A[5][XDIM][YDIM]; //D2Q5, X,Y  
double Source[XDIM][YDIM];
double rho[XDIM][YDIM];
double rhoPast[XDIM][YDIM];
double B[XDIM][YDIM][2];
double errorTotal;
double error[XDIM][YDIM];

double leftBuffer[velCount][YDIM], rightBuffer[velCount][YDIM];
double topBuffer[velCount][XDIM], bottomBuffer[velCount][XDIM];

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double slice[XDIM];
double sliceTh[XDIM];
double rhoTh[XDIM][YDIM];
int sliceLoc = 33;
double k = 0.6332573978; // 100/(4pi)^2 when J = 1
int rhoThreq=0;
int sliceReq=0;

//Function Declaration
void Collision();
void Initialize();
void Stream4();
void Curl();
void Iterate();
void GetSlice();


void Initialize(){


	Source[17][17] = 1;
	Source[34][17] = -1;
	Source[17][34] = -1;	
	Source[34][34] = 1;		

	for(int x=0; x<(XDIM-1); x++){
		for(int y=0; y<(YDIM-1); y++){
			for(int i=0; i<5; i++){

				A[i][x][y] = 0; //velDist[i] + (velDist[i] * Source[x][y]);

			}
		}
	}		

    for (int x=0;x<(XDIM-1);x++){
      	for (int y=0;y<(YDIM-1);y++){
      		rho[x][y] = 0;
      		//B[0][x][y] = 0;
      		//B[1][x][y] = 0;
		}
  	}

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

}


void Collision(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int i=0; i<5; i++){
				//Finew = Fi + (1/T * (Fi0 - Fi)) + Source
				A[i][x][y] = A[i][x][y] + ((1./T) * ((velDist[i] * rho[x][y]) - A[i][x][y])) + Source[x][y] * velDist[i];


			}
		}
	}
}


void Stream4(){

	//Velocity 1
	for(int y=YDIM-1; y>=0; y--){
		leftBuffer[1][y] = A[1][XDIM-1][y];
	}
	memmove(&A[1][1][0], &A[1][0][0], (sizeof(double)*(XDIM-1)*(YDIM)));	

	//Velocity 2
	for(int x=XDIM-1; x>=0; x--){
		bottomBuffer[2][x] = A[2][x][YDIM-1];
	}
	memmove(&A[2][0][1], &A[2][0][0], (sizeof(double)*(XDIM)*(YDIM-1)));

	//Velocity 3
	for(int y=0; y<YDIM; y++){
		rightBuffer[3][y] = A[3][0][y];
	}
	memmove(&A[3][0][0], &A[3][1][0], (sizeof(double)*(XDIM-1)*(YDIM)));

	//Velocity 4
	for(int x=0; x<XDIM; x++){
		topBuffer[4][x] = A[4][x][0];
	}
	memmove(&A[4][0][0], &A[4][0][1], (sizeof(double)*(XDIM)*(YDIM-1)));

	//Merge velocities from buffers
	for(int x=0; x<XDIM; x++){
		A[2][x][0] = bottomBuffer[2][x];
		A[4][x][YDIM-1] = topBuffer[4][x];
	}

	for(int y=0; y<YDIM; y++){
		A[1][0][y] = leftBuffer[1][y];
		A[3][XDIM-1][y] = rightBuffer[3][y];
	}
}


void Density(){

    for (int x=0;x<XDIM;x++){
      	for (int y=0;y<YDIM;y++){

      		rho[x][y] = 0;
			for(int i=0; i<5; i++){
				rho[x][y] += A[i][x][y];
			}
		}
  	}
}


void Curl(){

	int dx;
	int dy;
	
    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){

      			dy = rho[x][y] - rho[x][y-1];	//Az
      			dx = rho[x][y] - rho[x-1][y];	//Az

      			B[x][y][0] = dy;
      			B[x][y][1] = -dx;
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

	double r[4];
	double qTest = 1;
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
										distance = sqrt(((x-xx+xOff)*(x-xx+xOff)) + ((y-yy+yOff)*(y-yy+yOff)) + ((0-zz+zOff)*(0-zz+zOff)));
									}
																
									rhoTh[x][y] += (Source[xx][yy] / distance) / (4 * 3.1415926 * 0.25 * (T - 0.5));  // (J/d) / (4pi * 0.25 * M) 
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


void Iterate(){

	Density();	
	Collision();
	Stream4();
	Curl();
	ComputeError();
	
	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	//DeltaTime();
	
}

void GetGraph(){

	if(rhoThreq && sliceReq){ //ANDed because otherwise it constantly calculates
		GetSlice();
		rhoThreq = 0;
	}

	if(sliceReq){
		sliceReq=0;
		for(int x=0; x<XDIM; x++){
			slice[x] = rho[x][sliceLoc];			
			sliceTh[x] = k * rhoTh[x][sliceLoc];
		}
	}

}

void ComputeError(){

	errorTotal=0;
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
		error[x][y] = (rho[x][y] - rhoPast[x][y])*(rho[x][y] - rhoPast[x][y]);
		errorTotal += error[x][y];
		rhoPast[x][y] = rho[x][y];
		}
	}	

	printf("Error Total = %d\n", errorTotal);

}

//-----Main Routines-----
int main(){
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=100;
  	int drawCount=0;


	DefineGraphNxN_R("XY Density",&rho[0][0],&xDim,&yDim,NULL);  
	DefineGraphNxN_R("Theoretical Density",&rhoTh[0][0],&xDim,&yDim,NULL);  		
	DefineGraphNxN_RxR("B Field(X,Y)",&B[0][0][0],&xDim,&yDim,NULL);
	DefineGraphN_R("Rho(x)",&slice[0],&xDim,NULL);
	DefineGraphN_R("RhoTh(x)",&sliceTh[0],&xDim,&sliceReq);
	DefineGraphN_R("Difference(x)",&error[0][sliceLoc],NULL);


	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineDouble("Tau", &T);
    DefineInt("Slice", &sliceLoc);
    DefineBool("GetSlice",&rhoThreq);
    DefineDouble("k", &k);
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