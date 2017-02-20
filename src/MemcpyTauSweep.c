#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define XDIM 50
#define YDIM 50
#define velCount 5

int xDim = XDIM;
int yDim = YDIM;
//int velCount = 5;
int iterationCount = 0;
int cpuCount=0;

double T = 1;
double velDist[5] = {0.5,0.125,0.125,0.125,0.125};

double A[5][XDIM][YDIM]; //D2Q5, X,Y  
double Source[5][XDIM][YDIM];
double rho[XDIM][YDIM];
double B[XDIM][YDIM];
double dB[XDIM][YDIM];
double errorB;
int errorCount;
double slice[XDIM];
double sliceTh[XDIM];
double slope[XDIM];
double sourceVal = 10;
int sourcePos1 = 17;
int sourcePos2 = 34;
int sliceThrequest;



double leftBuffer[velCount][YDIM], rightBuffer[velCount][YDIM];
double topBuffer[velCount][XDIM], bottomBuffer[velCount][XDIM];

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;


//Main loop variables
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=1;
  	int drawCount=0;


//Function Declaration
void Collision();
void Initialize();
void Stream4();
void Stream2();
void Curl();
void Iterate();
void GetSlope();
void GetSlice();
void PrintSlice();


void Initialize(){

	cont = 1;

	for(int i=1; i<5; i++){
		for(int y=0; y<YDIM; y++){
			for(int x=0; x<xDim; x++){
				Source[i][x][y] = 0;
				Source[i][x][y] = 0;
			}
		}
	}

	
	for(int i=1; i<5; i++){
		for(int j=0; j<YDIM; j++){
			Source[i][sourcePos1][j] = sourceVal * velDist[i];
			Source[i][sourcePos2][j] = -sourceVal * velDist[i];
		}
	}
	
	/*
	for(int j=0; j<YDIM; j++){	
		Source[2][17][j] = 10 * velDist[2];
		Source[4][34][j] = -10 * velDist[4];
	}
	*/


	for(int x=0; x<(XDIM-1); x++){
		for(int y=0; y<(YDIM-1); y++){
			for(int i=0; i<5; i++){

				A[i][x][y] = velDist[i];// + Source[i][x][y];

			}
		}
	}		

    for (int x=0;x<(XDIM-1);x++){
      	for (int y=0;y<(YDIM-1);y++){
      		rho[x][y] = 0;
      		B[x][y] = 0;
		}
  	}

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

}


void Collision(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
		//printf("rho[%d][%d]= %10.10f\n", x,y,rho[x][y]);	
			for(int i=0; i<5; i++){

				//if(x==25 && y==25){
				//	printf("A[%d][%d][%d]= %10.10f\n", i,x,y,A[i][x][y]);
				//}

				//Finew = Fi + (1/T * (Fi0 - Fi)) + Source
				A[i][x][y] = A[i][x][y] + ((1./T) * ((velDist[i] * rho[x][y]) - A[i][x][y])) + Source[i][x][y];

				//if(x==25 && y==25){
				//	printf("VeloDist= %10.10f\n", velDist[i]);
				//	printf("Anew[%d][%d][%d]= %10.10f\n", i,x,y,A[i][x][y]);
				//}

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
	errorB = 0;
	
    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){

      		dx = rho[x][y] - rho[x-1][y];
      		dy = rho[x][y] - rho[x][y-1];

      		dB[x][y] = B[x][y];
      		B[x][y] = dx - dy;
      		dB[x][y] = dB[x][y] - B[x][y];
      		//printf("dB = %10.10f\n", dB[x][y]);
      		errorB += abs(dB[x][y]);
		}
  	}	 
  	printf("B error = %10.10f\n", errorB);
}


void deltaB(){

	if((abs(errorB) < .0001) && iterationCount > 10){
		errorCount++;
		if(errorCount > 100){
			cont = 0;
			printf("B Field has stabilized.\n");
			GetSlope();
			PrintSlice();
		}
	}
	else{
		errorCount = 0;
	}

}


void GetSlice(){

	int y = YDIM / 2;
	for(int x=0; x<XDIM-1; x++){

		slice[x] = rho[x][y];

	}
}


void PrintSlice(){

	printf("Slice[");
	for(int x=0; x<XDIM-1; x++){
		printf(" %5.5f,",slice[x]);
	}	
	printf("]\n");
}


void GetSlope(){

	for(int x=0; x<XDIM-1; x++){
		slope[x] = slice[x+1] - slice[x]; // dy / dx, dx = 1
	}

	printf("Slope[");
	for(int x=0; x<XDIM-1; x++){
		printf(" %5.5f,",slope[x]);
	}	
	printf("]\n");

}

void TauSweep(double start, double end, double step){



}

void TheoryLine(){
	double mu;
	int dx1 = sourcePos2-sourcePos1;
	int dx2 = sourcePos1 + XDIM - sourcePos2;

	mu = sourceVal / (/*Removed (2*)*/(T - .5)*( 1./dx1 + 1./dx2 ));

	for(int i=sourcePos1; i<=sourcePos2; i++){
		sliceTh[i] = mu - (i-sourcePos1) * 2*(mu/dx1);
	}

	for(int i=sourcePos2; i<XDIM; i++){
		sliceTh[i] = -mu + (i-sourcePos2) * 2*(mu/dx2);
	}

	for(int i=0; i<=sourcePos1; i++){
		sliceTh[sourcePos1 - i] = mu - i * 2*(mu/dx2);
	}		


}


void Iterate(){

	Density();	
	Collision();
	Stream4();
	//deltaB();
	Curl();
	GetSlice();
	
	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	
}

//-----Main Routines-----
int main(){


	DefineGraphNxN_R("XY Density",&rho[0][0],&xDim,&yDim,NULL);  	
	DefineGraphNxN_R("B Field(Z)",&B[0][0],&xDim,&yDim,NULL);
	DefineGraphN_R("Rho vs X",&slice[0],&xDim,NULL);
	DefineGraphN_R("TheoryLine",&sliceTh[0],&xDim,&sliceThrequest);

	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    StartMenu("Init",0);
    DefineInt("Source 1",&sourcePos1);
    DefineInt("Source 2",&sourcePos2);
    DefineDouble("Source Value", &sourceVal);
    DefineFunction("Reset", &Initialize);
    EndMenu();
    DefineInt("Repeat", &repeat);
    DefineDouble("Tau", &T);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");
  	DefineBool("Quit",&done);
  	EndMenu();

	Initialize();

  	while (!done){
  		if(drawCount>repeat){
      		Events(1);
      	    if(sliceThrequest){
	  			TheoryLine();
	  			sliceThrequest = 0;
	    	}
      		DrawGraphs();
      		drawCount = 0;
  		}
	    if(cont){
	    	Iterate();
		}

	    drawCount++;
  }
}
