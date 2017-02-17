#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define pi 3.141592653
#define MPERM 4 * pi * 0.0000001
#define XDIM 50
#define YDIM 50
#define ZDIM 100
#define velCount 7
#define u0 1
#define u .6		//Technically this is relative permeability when u0 is anything other than 4pi * 10^-7 


int xDim = XDIM;
int yDim = YDIM;
int zDim = ZDIM;
//int velCount = 5;
int iterationCount = 0;
int cpuCount=0;

double T = 1;
double velDist[velCount] = {0.25,0.125,0.125,0.125,0.125,0.125,0.125};

double A[3][velCount][XDIM][YDIM][ZDIM]; 	//Simulation,D3Q7, X,Y,Z  
double Source[3][XDIM][YDIM][ZDIM];	//Includes u0
double rho[3][XDIM][YDIM][ZDIM];
double B[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
double rhoDisplay[XDIM][YDIM];
double rhoYZDisplay[ZDIM][YDIM];
double BDisplay[XDIM][YDIM][2];	//x,y (location), and 2D vector 
double BYZDisplay[ZDIM][XDIM][2];	//z,y (location), and 2D vector 
double BZDisplay[XDIM][YDIM];
double mu[XDIM][YDIM][ZDIM];


double leftBuffer[3][velCount][YDIM][ZDIM]; 	//Velocity 1
double rightBuffer[3][velCount][YDIM][ZDIM];	//Velocity 2
double topBuffer[3][velCount][XDIM][ZDIM]; 	//Velocity 3
double bottomBuffer[3][velCount][XDIM][ZDIM];	//Velocity 4
double frontBuffer[3][velCount][XDIM][YDIM]; 	//Velocity 5
double backBuffer[3][velCount][XDIM][YDIM];	//Velocity 6

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double slice[XDIM];
double sliceTh[XDIM];
double rhoTh[XDIM][YDIM];
int sliceLoc = 33;
double k = 1;
int Simulation = 0;
int rhoThreq=0;
int sliceReq=0;
int zSlice = 25;
int ySlice = 25;
int xSlice = 25;

struct current {double x; double y; double z;} J[XDIM][YDIM][ZDIM]; //Current vector at a given location
double JMagnitude = .111111111;
int sourceOn = 1;


//Function Declaration
void Collision(int s);
void Initialize();
void Stream(int s);
void Curl();
void Iterate();
void GetSlice(int s);
void GenerateSourcePath(int cx, int cy, int cz, int radius, int thickness); //Coil center location (x,y,z), radius, and wire thickness
void Density(int s);
void Visualize(int s);
void ConfigurePermeability();
void GetFieldEnergy();



void Initialize(){

	for(int z=17; z<(ZDIM-17); z=z+3){
		int r=10;
		//for(int r=11; r<25; r=r+3){

			printf("(%d,%d)\n",z,r);
			GenerateSourcePath(25, 25, z, r, 3);
		//}
	}
	

	//GenerateSourcePath(25, 25, 25, 5, 3);

	ConfigurePermeability();

	for(int x=0; x<(XDIM); x++){
		for(int y=0; y<(YDIM); y++){
			for(int z=0; z<(ZDIM); z++){
												
				Source[0][x][y][z] = JMagnitude * J[x][y][z].x;
				Source[1][x][y][z] = JMagnitude * J[x][y][z].y;
				Source[2][x][y][z] = JMagnitude * J[x][y][z].z;
			}
		}
	}	


	for(int s=0; s<3; s++){
		for(int x=0; x<(XDIM); x++){
			for(int y=0; y<(YDIM); y++){
				for(int z=0; z<(ZDIM); z++){
					for(int i=0; i<velCount; i++){

						A[s][i][x][y][z] = 0; //velDist[i] + Source[i][x][y];

						leftBuffer[s][i][y][z] = 0; 	//Velocity 1
						rightBuffer[s][i][y][z] = 0;	//Velocity 2
						topBuffer[s][i][x][z] = 0; 	//Velocity 3
						bottomBuffer[s][i][x][z] = 0;	//Velocity 4
						frontBuffer[s][i][x][y] = 0; 	//Velocity 5
						backBuffer[s][i][x][y] = 0;	//Velocity 6					
					}
				}
			}
		}	
	}	

	for(int s=0; s<3; s++){
	    for (int x=0;x<(XDIM);x++){
	      	for (int y=0;y<(YDIM);y++){
				for(int z=0; z<(ZDIM); z++){
			      		rho[s][x][y][z] = 0;
				}
			}
	  	}
  	}

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

}


void GenerateSourcePath(int cx, int cy, int cz, int radius, int thickness){

	int dx;
	int dy;		
	int xPos;
	int yPos;

	//Create the "wire" for current to flow through
	for(double theta=0; theta<(2*pi); theta+=.01){
		for(int r=(radius-(thickness/2)); r<(radius+((thickness+1)/2)); r++){
			for(int z=(cz-1); z<(cz+2); z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				//2D, 90 degree rotation in the x-y plane
				J[xPos][yPos][z].x = -dy;
				J[xPos][yPos][z].y = dx;
				J[xPos][yPos][z].z = 0;
			}				
		}
	}
}

void ConfigurePermeability(){

	int cx = 25;
	int cy = 25;
	int cz = 25;
	int radius = 10;

	int dx;
	int dy;		
	int xPos;
	int yPos;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				mu[x][y][z] = u0;
			}
		}
	}

	
	for(double theta=0; theta<(2*pi); theta+=.1){
		for(int r=0; r<(radius); r++){
			for(int z=5; z<95; z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				mu[xPos][yPos][z] = u;

			}				
		}
	}
	
	
	

}


void Collision(int s){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				for(int i=0; i<velCount; i++){
					//Finew = Fi + (1/T * (Fi0 - Fi)) + Source, where T = mu +.5
					A[s][i][x][y][z] = A[s][i][x][y][z] 
										+ (((mu[x][y][z] / T)) * ((velDist[i] * rho[s][x][y][z]) - A[s][i][x][y][z])) 
										+ Source[s][x][y][z] * velDist[i] * sourceOn;
				}
			}
		}
	}
}


void Stream(int s){

	//Velocity 1
	for(int y=YDIM-1; y>=0; y--){
		for(int z=0; z<ZDIM; z++){
			leftBuffer[s][1][y][z] = A[s][1][XDIM-1][y][z];
		}
	}
	memmove(&A[s][1][1][0][0], &A[s][1][0][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));	

	//Velocity 2
	for(int x=XDIM-1; x>=0; x--){
		for(int z=0; z<ZDIM; z++){
			bottomBuffer[s][2][x][z] = A[s][2][x][YDIM-1][z];
		}	
	}
	memmove(&A[s][2][0][1][0], &A[s][2][0][0][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 3
	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			rightBuffer[s][3][y][z] = A[s][3][0][y][z];
		}
	}
	memmove(&A[s][3][0][0][0], &A[s][3][1][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));

	//Velocity 4
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			topBuffer[s][4][x][z] = A[s][4][x][0][z];
		}
	}
	memmove(&A[s][4][0][0][0], &A[s][4][0][1][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 5
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			frontBuffer[s][5][x][y] = A[s][5][x][y][0];
		}
	}
	memmove(&A[s][5][0][0][1], &A[s][5][0][0][0], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));

	//Velocity 6
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			backBuffer[s][5][x][y] = A[s][5][x][y][ZDIM-1];
		}
	}
	memmove(&A[s][6][0][0][0], &A[s][6][0][0][1], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));	


	//Merge velocities from buffers
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			A[s][2][x][0][z] = bottomBuffer[s][2][x][z];
			A[s][4][x][YDIM-1][z] = topBuffer[s][4][x][z];
		}
	}

	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			A[s][1][0][y][z] = leftBuffer[s][1][y][z];
			A[s][3][XDIM-1][y][z] = rightBuffer[s][3][y][z];
		}
	}

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			A[s][5][x][y][0] = frontBuffer[s][5][x][y];
			A[s][6][x][y][ZDIM-1] = backBuffer[s][6][x][y];
		}
	}	
}


void Density(int s){

    for (int x=0;x<XDIM;x++){
      	for (int y=0;y<YDIM;y++){
			for(int z=0; z<ZDIM; z++){
	      		rho[s][x][y][z] = 0;
				for(int i=0; i<velCount; i++){
					rho[s][x][y][z] += A[s][i][x][y][z];
				}
			}
		}
  	}
}

void Curl(){

	int i;
	int j;
	int k;

    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){
			for(int z=1; z<(ZDIM-1); z++){

      			i = (rho[2][x][y][z] - rho[2][x][y-1][z]) - (rho[1][x][y][z] - rho[1][x][y][z-1]);
      			j = (rho[0][x][y][z] - rho[0][x][y][z-1]) - (rho[2][x][y][z] - rho[2][x-1][y][z]);
      			k = (rho[1][x][y][z] - rho[1][x-1][y][z]) - (rho[0][x][y][z] - rho[0][x][y-1][z]);

      			B[x][y][z][0] = i;
      			B[x][y][z][1] = j;
      			B[x][y][z][2] = k;
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

void GetSlice(int s){

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
																
									rhoTh[x][y] += mu[xx][yy][zz] * (Source[s][xx][yy][zz] / distance) / (3.1415926 * (T));  // (J/d) / (4pi * 0.25 * M) 
								}
							}
						}
					}
				}
			}
		//}
		printf("X=%d\n",x);		
	}
}



void GetGraph(int s){

	if(rhoThreq && sliceReq){ //ANDed because otherwise it constantly calculates
		GetSlice(s);
		rhoThreq = 0;
	}

	if(sliceReq){
		sliceReq=0;
		for(int x=0; x<XDIM; x++){
			slice[x] = rho[s][x][sliceLoc][zSlice];			
			sliceTh[x] = k * rhoTh[x][sliceLoc];
		}
	}

}

void Visualize(int s){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){

				rhoDisplay[x][y] = rho[s][x][y][zSlice];
				rhoYZDisplay[z][y] = rho[s][xSlice][y][z]; //
				
				//Projection into 2D
				BDisplay[x][y][0] = B[x][y][zSlice][0];
				BDisplay[x][y][1] = B[x][y][zSlice][1];
				BYZDisplay[z][y][0] = B[xSlice][y][z][2]; //
				BZDisplay[x][y] = B[x][y][zSlice][2];
				BYZDisplay[z][y][1] = B[xSlice][y][z][1]; //
			}
		}
	}
}


void GetFieldEnergy(){

	double eTotal = 0;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
			
							// 		sum(B^2 * dv * 1/(2u0))																	//Scale Factor to 1mm * factor for B field
				eTotal += (((B[x][y][z][0] * B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2])) * 0.000000001) / (2*u0);


			}
		}
	}

	printf("Total Field Energy = %10.10f\n", eTotal);
}


void Iterate(){

	Density(0);	
	Density(1);	
	Density(2);		
	Collision(0);
	Collision(1);
	Collision(2);		
	Stream(0);
	Stream(1);
	Stream(2);		
	Curl();
	GetFieldEnergy();
	
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
	DefineGraphNxN_R("B Field(Z)",&BZDisplay[0][0],&xDim,&yDim,NULL);
	DefineGraphNxN_R("YZ Density",&rhoYZDisplay[0][0],&zDim,&yDim,NULL);
	DefineGraphNxN_RxR("B Field(Y,Z)",&BYZDisplay[0][0][0],&zDim,&yDim,NULL);	
	DefineGraphN_R("Rho(x)",&slice[0],&xDim,NULL);
	DefineGraphN_R("RhoTh(x)",&sliceTh[0],&xDim,&sliceReq);


	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineDouble("Tau", &T);
    DefineInt("Simulation", &Simulation);
    DefineInt("Slice", &sliceLoc);
    DefineBool("GetSlice",&rhoThreq);
    DefineDouble("J", &JMagnitude);
    DefineBool("Source", &sourceOn);
    DefineInt("Z plane", &zSlice);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");    
  	DefineBool("Quit",&done);
  	EndMenu();

	Initialize();

  	while (!done){
  		if(drawCount>repeat){
      		Events(1);
      		Visualize(Simulation);
      		GetGraph(Simulation);
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