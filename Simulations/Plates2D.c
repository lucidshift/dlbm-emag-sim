#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define XDIM 50
#define YDIM 50
#define ZDIM 1
#define velCount 5

int xDim = XDIM;
int yDim = YDIM;
//int velCount = 5;
int iterationCount = 0;
int cpuCount=0;


double T = 1;
double velDist[5] = {0.5,0.125,0.125,0.125,0.125};

double A[3][5][XDIM][YDIM]; //D2Q5, 3 sims, 5 velo, X,Y  
double Source[3][5][XDIM][YDIM]; // 3 sims, 5 velo, x,y,
double rho[XDIM][YDIM][3];
double B[XDIM][YDIM][3]; //x, y, 3 simulations

//Arrays for visualizing
double rhoX[XDIM][YDIM];
double rhoY[XDIM][YDIM];
double rhoZ[XDIM][YDIM];
double BXY[XDIM][YDIM][2];
//double BY[XDIM][YDIM][2];

double BZ[XDIM][YDIM];	//testing this

double leftBuffer[velCount][YDIM][3], rightBuffer[velCount][YDIM][3];
double topBuffer[velCount][XDIM][3], bottomBuffer[velCount][XDIM][3];

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double slice[XDIM];
double sliceTh[XDIM];
double rhoTh[XDIM][YDIM];
int sliceLoc = 33;
double k = 1.5;
int rhoThreq=0;
int sliceReq=0;


//Function Declaration
void Collision();
void Initialize();
void Stream4(int idx);
void Curl();
void Iterate();
void Density(int idx);
void Visualize();
void GetSlice();


void Initialize(){


	for(int i=1; i<5; i++){
		for(int j=17; j<35; j++){

			//Vertical Plates
			Source[0][i][17][j] = 10 * velDist[i];
			Source[0][i][34][j] = -10 * velDist[i];

			//Horizontal Plates
			Source[1][i][j][17] = -10 * velDist[i];
			Source[1][i][j][34] = 10 * velDist[i];			
		}
	}


	for(int x=0; x<(XDIM-1); x++){
		for(int y=0; y<(YDIM-1); y++){
			for(int i=0; i<5; i++){
				for(int j=0; j<3; j++){

					A[j][i][x][y] = velDist[i];// + Source[i][x][y];
				}
			}
		}
	}		

    for (int x=0;x<(XDIM-1);x++){
      	for (int y=0;y<(YDIM-1);y++){
      		for(int i=0; i<3; i++){

      			rho[x][y][i] = 0;
      			B[x][y][i] = 0;

      		}
		}
  	}

    for (int i=0;i<5;i++){
      	for (int y=0;y<(YDIM-1);y++){
      		for(int j=0; j<3; j++){
	      		leftBuffer[i][y][j] = 0;
	      		rightBuffer[i][y][j] = 0;
	      		topBuffer[i][y][j] = 0;
	      		bottomBuffer[i][y][j] = 0;
      		}
		}
  	}  	

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);
}


void Collision(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){	
			for(int i=0; i<5; i++){
				for(int j=0; j<2; j++){		//2 simulations for now, 3 later

					//Finew = Fi + (1/T * (Fi0 - Fi)) + Source
					A[j][i][x][y] = A[j][i][x][y] + ((1./T) * ((velDist[i] * rho[x][y][j]) - A[j][i][x][y])) + Source[j][i][x][y];
				}
			}
		}
	}
}


void Stream4(int idx){

	//Velocity 1
	for(int y=YDIM-1; y>=0; y--){
		leftBuffer[1][y][idx] = A[idx][1][XDIM-1][y];
	}
	memmove(&A[idx][1][1][0], &A[idx][1][0][0], (sizeof(double)*(XDIM-1)*(YDIM)));	

	//Velocity 2
	for(int x=XDIM-1; x>=0; x--){
		bottomBuffer[2][x][idx] = A[idx][2][x][YDIM-1];
	}
	memmove(&A[idx][2][0][1], &A[idx][2][0][0], (sizeof(double)*(XDIM)*(YDIM-1)));

	//Velocity 3
	for(int y=0; y<YDIM; y++){
		rightBuffer[3][y][idx] = A[idx][3][0][y];
	}
	memmove(&A[idx][3][0][0], &A[idx][3][1][0], (sizeof(double)*(XDIM-1)*(YDIM)));

	//Velocity 4
	for(int x=0; x<XDIM; x++){
		topBuffer[4][x][idx] = A[idx][4][x][0];
	}
	memmove(&A[idx][4][0][0], &A[idx][4][0][1], (sizeof(double)*(XDIM)*(YDIM-1)));

	//Merge velocities from buffers
	for(int x=0; x<XDIM; x++){
		A[idx][2][x][0] = bottomBuffer[2][x][idx];
		A[idx][4][x][YDIM-1] = topBuffer[4][x][idx];
	}

	for(int y=0; y<YDIM; y++){
		A[idx][1][0][y] = leftBuffer[1][y][idx];
		A[idx][3][XDIM-1][y] = rightBuffer[3][y][idx];
	}
}

void Density(int idx){

    for (int x=0;x<XDIM;x++){
      	for (int y=0;y<YDIM;y++){

      		rho[x][y][idx] = 0;
			for(int i=0; i<5; i++){
				rho[x][y][idx] += A[idx][i][x][y];
			}
		}
  	}
}

void Curl(){

	int dx;
	int dy;
	//int dz;
	
    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){

  			dy = rho[x][y][0] - rho[x][y-1][0];	
  			dx = rho[x][y][0] - rho[x-1][y][0];	
  			//Z direction
  			B[x][y][2] = dy-dx;			//This is equals!


  			dy = rho[x][y][1] - rho[x][y-1][1];	
  			dx = rho[x][y][1] - rho[x-1][y][1];	
  			//Z direction
  			B[x][y][2] += dy-dx;		//This is plus equals!      			

		}
  	}	 
}


void Visualize(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){

			rhoX[x][y] = rho[x][y][0];
			rhoY[x][y] = rho[x][y][1];

			BZ[x][y] = B[x][y][2];
			//BXY[x][y][2] = B[x][y][2];

		}
	}
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
																
									rhoTh[x][y] += (Source[0][1][xx][yy] / distance);
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


void DeltaTime(){

	struct timespec spec;
	clock_gettime(CLOCK_REALTIME, &spec);

	t.past2 = t.past1;
	t.past1 = t.now;
	t.now = spec.tv_nsec / 1000000.0;
	t.dt = t.now - t.past1;
	printf ("deltaTime = %2.5f ms\n", t.dt);

}

void GetGraph(){

	if(rhoThreq && sliceReq){ //ANDed because otherwise it constantly calculates
		GetSlice();
		rhoThreq = 0;
	}

	if(sliceReq){
		sliceReq=0;
		for(int x=0; x<XDIM; x++){
			slice[x] = rho[x][sliceLoc][0];			
			sliceTh[x] = k * rhoTh[x][sliceLoc];
		}
	}

}

void Iterate(){

	Density(0);
	Density(1);	
	Collision();
	Stream4(0);
	Stream4(1);
	Curl(0);
	Curl(1);
	
	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	DeltaTime();
	
}

//-----Main Routines-----
int main(){
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=100;
  	int drawCount=0;

	DefineGraphNxN_R("X Density",&rhoX[0][0],&xDim,&yDim,NULL);  	
	DefineGraphNxN_R("Y Density",&rhoY[0][0],&xDim,&yDim,NULL);  
	DefineGraphNxN_R("BZ Field(Curl of X,Y)",&BZ[0][0],&xDim,&yDim,NULL);
	DefineGraphN_R("Rho(x)",&slice[0],&xDim,NULL);
	DefineGraphN_R("RhoTh(x)",&sliceTh[0],&xDim,&sliceReq);

	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineDouble("Tau", &T);
    DefineInt("Slice", &sliceLoc);
    DefineInt("GetSlice",&rhoThreq);
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
      		Visualize();
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