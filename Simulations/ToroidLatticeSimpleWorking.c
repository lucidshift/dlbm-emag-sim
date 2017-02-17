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
double velDist[velCount] = {0.5,0.125,0.125,0.125,0.125};

double A[5][XDIM][YDIM]; //D2Q5, X,Y  
double Source[5][XDIM][YDIM];
double rho[XDIM][YDIM];
double B[XDIM][YDIM];

double leftBuffer[velCount][YDIM], rightBuffer[velCount][YDIM];
double topBuffer[velCount][XDIM], bottomBuffer[velCount][XDIM];

struct offset {int x; int y;};
struct offset off[velCount];

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

//Function Declaration
void Collision();
void Initialize();
void Stream();
void Stream2();
void Curl();
void Iterate();
double GetTime();
void DeltaTime();


void Initialize(){

	for(int i=1; i<5; i++){
		for(int j=0; j<yDim; j++){
			Source[i][17][j] = 10 * velDist[i];
			Source[i][34][j] = -10 * velDist[i];
		}
	}

	for(int x=0; x<(xDim-1); x++){
		for(int y=0; y<(yDim-1); y++){
			for(int i=0; i<5; i++){

				A[i][x][y] = velDist[i];// + Source[i][x][y];

			}
		}
	}		

    for (int x=0;x<(xDim-1);x++){
      	for (int y=0;y<(yDim-1);y++){
      		rho[x][y] = 0;
      		B[x][y] = 0;
		}
  	}
	
	for(int i=0; i<5; i++){

		off[i].x = 0;
		off[i].y = 0;
	}

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

}


void Collision(){

	for(int x=0; x<xDim; x++){
		for(int y=0; y<yDim; y++){
		//printf("rho[%d][%d]= %10.10f\n", x,y,rho[x][y]);	
			for(int i=0; i<5; i++){

				//if(x==25 && y==25){
				//	printf("A[%d][%d][%d]= %10.10f\n", i,x,y,A[i][x][y]);
				//}

				//Finew = Fi + (1/T * (Fi0 - Fi)) + Source
				A[i][(x + off[i].x) % (xDim-1)][(y + off[i].y) % (yDim-1)] = 
				A[i][(x + off[i].x) % (xDim-1)][(y + off[i].y) % (yDim-1)] + ((1./T) * ((velDist[i] * rho[x][y]) - 
				A[i][(x + off[i].x) % (xDim-1)][(y + off[i].y) % (yDim-1)])) + Source[i][x][y];

				//printf("offX = %d\n", off[i].x);
				//printf("offY = %d\n", off[i].y);

				//if(x==25 && y==25){
				//	printf("VeloDist= %10.10f\n", velDist[i]);
				//	printf("Anew[%d][%d][%d]= %10.10f\n", i,x,y,A[i][x][y]);
				//}

			}
		}
	}
}


void Stream(){

	for(int x=1; x<xDim-2; x++){
		for(int y=1; y<yDim-2; y++){

			A[1][(xDim-1)-x][y] = A[1][(xDim-1)-x-1][y];
			A[2][x][(yDim-1)-y] = A[2][x][(yDim-1)-y-1];
			A[3][x][y] = A[3][x+1][y];
			A[4][x][y] = A[4][x][y+1];

		}
	}

	//x=0
	for(int y=1; y<yDim-2; y++){

		A[1][0][y] = A[1][xDim-1][y];		//Wrapping Condition
		A[2][0][(yDim-1)-y] = A[2][0][(yDim-1)-y-1];
		A[3][0][y] = A[3][0+1][y];
		A[4][0][y] = A[4][0][y+1];

	}


	//x=xDim-1
	for(int y=1; y<yDim-2; y++){

		A[1][xDim-1][y] = A[1][xDim-1-1][y];
		A[2][xDim-1][(yDim-1)-y] = A[2][xDim-1][(yDim-1)-y-1];
		A[3][xDim-1][y] = A[3][0][y];			//Wrapping Condition
		A[4][xDim-1][y] = A[4][xDim-1][y+1];

	}

	//y=0
	for(int x=1; x<xDim-2; x++){

		A[1][(xDim-1)-x][0] = A[1][(xDim-1)-x-1][0];
		A[2][x][0] = A[2][x][yDim-1];		//Wrapping Condition
		A[3][x][0] = A[3][x+1][0];
		A[4][x][0] = A[4][x][0+1];

	}

	//y=yDim-1
	for(int x=1; x<xDim-2; x++){

		A[1][(xDim-1)-x][yDim-1] = A[1][(xDim-1)-x-1][yDim-1];
		A[2][x][yDim-1] = A[2][x][yDim-1-1];
		A[3][x][yDim-1] = A[3][x+1][yDim-1];
		A[4][x][yDim-1] = A[4][x][0]; 		//Wrapping Condition
	
	}

	//Corners
	A[1][0][0] = A[1][xDim-1][0];
	A[2][0][0] = A[2][0][yDim-1];	

	A[1][0][yDim-1] = A[1][xDim-1][yDim-1];	
	A[4][0][yDim-1] = A[4][0][0];

	A[3][xDim-1][yDim-1] = A[3][0][yDim-1];
	A[4][xDim-1][yDim-1] = A[4][xDim-1][0];	

	A[2][xDim-1][0] = A[2][xDim-1][yDim-1];	
	A[3][xDim-1][0] = A[3][0][0];



}

void Stream2(){

	//Velocity 1
	for(int x=xDim-1; x>0; x--){
		for(int y=yDim-1; y>0; y--){

			if(x==xDim-1){
				leftBuffer[1][y] = A[1][x][y];
			}
			else{
				A[1][x][y] = A[1][x-1][y];
			}
		}
	}

	//Velocity 2
	for(int x=xDim-1; x>0; x--){
		for(int y=yDim-1; y>0; y--){

			if(y==yDim-1){
				bottomBuffer[2][x] = A[2][x][y];
			}
			else{
				A[2][x][y] = A[2][x][y-1];
			}
		}
	}

	//Velocity 3
	for(int x=0; x<xDim-1; x++){
		for(int y=0; y<yDim-1; y++){

			if(x==0){
				rightBuffer[3][y] = A[3][x][y];
			}
			else{
				A[3][x][y] = A[3][x+1][y];
			}
		}
	}	

	//Velocity 4
	for(int x=0; x<xDim-1; x++){
		for(int y=0; y<yDim-1; y++){

			if(y==0){
				topBuffer[4][x] = A[4][x][y];
			}
			else{
				A[4][x][y] = A[4][x][y+1];
			}
		}
	}

	//Merge velocities from buffers
	for(int x=0; x<xDim; x++){

		A[2][x][0] = bottomBuffer[2][x];
		A[4][x][yDim-1] = topBuffer[4][x];
	}

	for(int y=0; y<yDim; y++){

		A[1][0][y] = leftBuffer[1][y];
		A[3][xDim-1][y] = rightBuffer[3][y];
	}
}

void Stream3(){

	//off[0].x = 0;

	if(off[1].x-- < 0){
		off[1].x = xDim - 1;
	}


	if(off[2].y-- < 0){
		off[2].y = yDim - 1;
	}
	

	if(off[3].x++ > xDim-1){
		off[3].x = 0;
	}


	if(off[4].y++ > yDim-1){
		off[4].y = 0;
	}


}

void Density(){

    for (int x=0;x<xDim;x++){
      	for (int y=0;y<yDim;y++){

      		rho[x][y] = 0;
			for(int i=0; i<5; i++){
				rho[x][y] += A[i][(x + off[i].x) % (xDim-1)][(y + off[i].y) % (yDim-1)];
			}
		}
  	}
}

void Curl(){

	 

}

void ShearFlowSetup(){

	for(int x=0; x<(xDim-1)/2; x++){
		for(int y=0; y<(yDim-1)/2; y++){
			for(int i=0; i<5; i++){

				A[3][x][y] = 5;// + Source[i][x][y];

			}
		}
	}	

	for(int x=(xDim-1)/2; x<(xDim-1); x++){
		for(int y=(yDim-1)/2; y<(yDim-1); y++){
			for(int i=0; i<5; i++){

				A[1][x][y] = 5;// + Source[i][x][y];

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


void Iterate(){

	Collision();
	Stream3();
	Curl();
	Density();

	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	DeltaTime();

}

//-----Main Routines-----
int main(){
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=1;
  	int drawCount=0;

	DefineGraphNxN_R("XY Density",&rho[0][0],&xDim,&yDim,NULL);  	

	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineGraph(contour2d_,"Densities");
  	DefineBool("Quit",&done);
  	EndMenu();

	Initialize();

  	while (!done){
  		if(drawCount>repeat){
      		Events(1);
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