#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <pthread.h>
#include <stdio.h>

#define pi 3.141592653
#define MPERM 4 * pi * 0.0000001
#define XDIM 50
#define YDIM 50
#define ZDIM 100
#define velCount 7
#define u0 1
#define u 1.2		//Technically this is relative permeability when u0 is anything other than 4pi * 10^-7 
#define SIZEX 4000


int xDim = XDIM;
int yDim = YDIM;
int zDim = ZDIM;
//int velCount = 5;
int iterationCount = 0;
int cpuCount=0;

float T = 1;
float velDist[velCount] = {0.25,0.125,0.125,0.125,0.125,0.125,0.125};

float A[3][velCount][XDIM][YDIM][ZDIM]; 	//Simulation,D3Q7, X,Y,Z  
float Source[3][XDIM][YDIM][ZDIM];	//Includes u0
float rho[3][XDIM][YDIM][ZDIM];
float B[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
float H[XDIM][YDIM][ZDIM][3];
double rhoDisplay[XDIM][YDIM];
double rhoYZDisplay[ZDIM][YDIM];
double BDisplay[XDIM][YDIM][2];	//x,y (location), and 2D vector 
double BYZDisplay[ZDIM][XDIM][2];	//z,y (location), and 2D vector 
double BZDisplay[XDIM][YDIM];
float mu[XDIM][YDIM][ZDIM];


float leftBuffer[3][velCount][YDIM][ZDIM]; 	//Velocity 1
float rightBuffer[3][velCount][YDIM][ZDIM];	//Velocity 2
float topBuffer[3][velCount][XDIM][ZDIM]; 	//Velocity 3
float bottomBuffer[3][velCount][XDIM][ZDIM];	//Velocity 4
float frontBuffer[3][velCount][XDIM][YDIM]; 	//Velocity 5
float backBuffer[3][velCount][XDIM][YDIM];	//Velocity 6

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double slice[XDIM];
double sliceTh[XDIM];
double rhoTh[XDIM][YDIM];
int sliceLoc = 33;
float k = 1;
int Simulation = 0;
int rhoThreq=0;
int sliceReq=0;
int zSlice = 25;
int ySlice = 25;
int xSlice = 25;

struct current {float x; float y; float z;} J[XDIM][YDIM][ZDIM]; //Current vector at a given location
float JInduced[3][XDIM][YDIM][ZDIM];
float JMagnitude = .111111111;
int sourceOn = 1;

double EMFDisplay[XDIM][YDIM][2];
struct volume {int xLow; int xHigh; int yLow; int yHigh; int zLow; int zHigh; int simulation;};

int sizeX = SIZEX;
double HUrCurve[SIZEX];	//Collection of x,y points
double BHCurve[SIZEX]; 
double muSlice[XDIM][YDIM];


//Function Declaration
void *Collision(void *threadBox);
void CollisionThreaded();
void Initialize();
void Stream(int s);
void Curl();
void Iterate();
void GetSlice(int s);
void GenerateSourcePath(int cx, int cy, int cz, int radius, int thickness); //Coil center location (x,y,z), radius, and wire thickness
void *Density(void *threadBox);
void DensityThreaded();
void *Visualize(void *threadBox);
void VisualizeThreaded(int s);
void ConfigurePermeability(float a);
void GetFieldEnergy();
void CalcInducedEMF();
void CalcHField();
void CalcMu();



void Initialize(){

	for(int z=17; z<(ZDIM-17); z=z+3){
		int r=10;
		//for(int r=11; r<25; r=r+3){

			printf("(%d,%d)\n",z,r);
			GenerateSourcePath(25, 25, z, r, 3);
		//}
	}
	
	ConfigurePermeability(u);

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

	for(int x=0; x<(XDIM); x++){

		BHCurve[x] = 0;
		HUrCurve[x] = 0;		
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
	for(float theta=0; theta<(2*pi); theta+=.01){
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

void ConfigurePermeability(float a){

	int cx = 25;
	int cy = 25;
	//int cz = 25;
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

	
	for(float theta=0; theta<(2*pi); theta+=.1){
		for(float r=0; r<(radius); r+=.5){
			for(int z=5; z<95; z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				mu[xPos][yPos][z] = a;

			}				
		}
		//printf("Rad= %5.5f\n",theta);
	}
	
	
	

}


void *Collision(void *threadBox){

	struct volume v = *((struct volume *) threadBox);
    free(threadBox);
    int s = v.simulation;

	for(int x= v.xLow; x< v.xHigh; x++){
		for(int y= v.yLow; y< v.yHigh; y++){
			for(int z= v.zLow; z< v.zHigh; z++){
				for(int i=0; i<velCount; i++){
					//Finew = Fi + (1/uT * (Fi0 - Fi)) + (Source - Source Induced)
					A[s][i][x][y][z] = A[s][i][x][y][z] 
										+ ((1./(mu[x][y][z] * T)) * ((velDist[i] * rho[s][x][y][z]) - A[s][i][x][y][z])) 
										+ (Source[s][x][y][z] + JInduced[s][x][y][z]) * velDist[i] * sourceOn;
				}
			}
		}
	}

	pthread_exit(0);
}


void Stream(int s){

	//Velocity 1
	for(int y=YDIM-1; y>=0; y--){
		for(int z=0; z<ZDIM; z++){
			leftBuffer[s][1][y][z] = A[s][1][XDIM-1][y][z];
		}
	}
	memmove(&A[s][1][1][0][0], &A[s][1][0][0][0], (sizeof(float)*(XDIM-1)*(YDIM)*(ZDIM)));	

	//Velocity 2
	for(int x=XDIM-1; x>=0; x--){
		for(int z=0; z<ZDIM; z++){
			bottomBuffer[s][2][x][z] = A[s][2][x][YDIM-1][z];
		}	
	}
	memmove(&A[s][2][0][1][0], &A[s][2][0][0][0], (sizeof(float)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 3
	for(int y=0; y<YDIM; y++){
		for(int z=0; z<ZDIM; z++){
			rightBuffer[s][3][y][z] = A[s][3][0][y][z];
		}
	}
	memmove(&A[s][3][0][0][0], &A[s][3][1][0][0], (sizeof(float)*(XDIM-1)*(YDIM)*(ZDIM)));

	//Velocity 4
	for(int x=0; x<XDIM; x++){
		for(int z=0; z<ZDIM; z++){
			topBuffer[s][4][x][z] = A[s][4][x][0][z];
		}
	}
	memmove(&A[s][4][0][0][0], &A[s][4][0][1][0], (sizeof(float)*(XDIM)*(YDIM-1)*(ZDIM)));

	//Velocity 5
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			frontBuffer[s][5][x][y] = A[s][5][x][y][0];
		}
	}
	memmove(&A[s][5][0][0][1], &A[s][5][0][0][0], (sizeof(float)*(XDIM)*(YDIM)*(ZDIM-1)));

	//Velocity 6
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			backBuffer[s][5][x][y] = A[s][5][x][y][ZDIM-1];
		}
	}
	memmove(&A[s][6][0][0][0], &A[s][6][0][0][1], (sizeof(float)*(XDIM)*(YDIM)*(ZDIM-1)));	


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


void *Density(void *threadBox){

	struct volume v = *((struct volume *) threadBox);
    free(threadBox);
    int s = v.simulation;

	for(int x= v.xLow; x< v.xHigh; x++){
		for(int y= v.yLow; y< v.yHigh; y++){
			for(int z= v.zLow; z< v.zHigh; z++){
	      		rho[s][x][y][z] = 0;
				for(int i=0; i<velCount; i++){
					rho[s][x][y][z] += A[s][i][x][y][z];
				}
			}
		}
  	}

  	pthread_exit(0);
}

void Curl(){

	float B2, H2;
	float muTemp;

    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){
			for(int z=1; z<(ZDIM-1); z++){

      			B[x][y][z][0] = (rho[2][x][y][z] - rho[2][x][y-1][z]) - (rho[1][x][y][z] - rho[1][x][y][z-1]);
      			B[x][y][z][1] = (rho[0][x][y][z] - rho[0][x][y][z-1]) - (rho[2][x][y][z] - rho[2][x-1][y][z]);
      			B[x][y][z][2] = (rho[1][x][y][z] - rho[1][x-1][y][z]) - (rho[0][x][y][z] - rho[0][x][y-1][z]);

      			muTemp = mu[x][y][z];
				B2 = sqrt((B[x][y][z][0] * B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2]));

				if(x==27 && y==27 && z==27){
					printf("|B|=%10.10f\n", B2);
					printf("Mu=%10.10f\n",muTemp);

					if(((iterationCount*2)+1)<sizeX){
						BHCurve[iterationCount*2] = H2*100;	
						BHCurve[iterationCount*2 + 1] = B2;
						HUrCurve[iterationCount*2] = H2;
						HUrCurve[iterationCount*2 + 1] = muTemp;
						printf("%5.10f,%5.10f\n", H2, B2);
					}					
				}								
				/*
				if(B2<0.0375 && muTemp != u0){		//Don't change the permeability of free space, it doesn't change.
					mu[x][y][z] = 1.3;
				}
				if(B2<0.5625 && B2>0.0375 && muTemp != u0){
					mu[x][y][z] = 1.25;
				}	
				if(B2<7.0 && B2>0.5625 && muTemp != u0){
					mu[x][y][z] = 1.13;
				}							
				if(B2<7.5 && B2>7.0 && muTemp != u0){
					mu[x][y][z] = 1.12;
				}
				if(B2>=7.5 && muTemp != u0){
					mu[x][y][z] = 1.1;
				}	
				*/
				H[x][y][z][0] = B[x][y][z][0] / mu[x][y][z];
				H[x][y][z][1] = B[x][y][z][1] / mu[x][y][z];
				H[x][y][z][2] = B[x][y][z][2] / mu[x][y][z];

				H2 = sqrt((H[x][y][z][0] * H[x][y][z][0]) + (H[x][y][z][1]*H[x][y][z][1]) + (H[x][y][z][2]*H[x][y][z][2]));

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

	float distance;
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
																
									rhoTh[x][y] += mu[xx][yy][zz] * ((Source[s][xx][yy][zz] + JInduced[s][xx][yy][zz])/ distance) / (3.1415926 * (T));  // (J/d) / (4pi * 0.25 * M) 
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

void *Visualize(void *threadBox){

	struct volume v = *((struct volume *) threadBox);
    free(threadBox);
    int s = v.simulation;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){

				rhoDisplay[x][y] = rho[s][x][y][zSlice];
				rhoYZDisplay[z][y] = rho[s][xSlice][y][z]; //

				EMFDisplay[x][y][0] = JInduced[0][x][y][zSlice] * 10;
				EMFDisplay[x][y][1] = JInduced[1][x][y][zSlice] * 10;

				muSlice[x][y] = mu[x][y][zSlice];
				
				//Projection into 2D
				BDisplay[x][y][0] = B[x][y][zSlice][0];
				BDisplay[x][y][1] = B[x][y][zSlice][1];
				BYZDisplay[z][y][0] = B[xSlice][y][z][2]; //
				BZDisplay[x][y] = B[x][y][zSlice][2];
				BYZDisplay[z][y][1] = B[xSlice][y][z][1]; //
			}
		}
	}

	pthread_exit(0);
}


void GetFieldEnergy(){

	float eTotal = 0;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
			
							// 		sum(B^2 * dv * 1/(2u0))																	//Scale Factor to 1mm * factor for B field
				eTotal += (((B[x][y][z][0] * B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2])) * 0.000000001) / (2*u0);


			}
		}
	}

	//printf("Total Field Energy = %10.10f\n", eTotal);
}


void CalcInducedEMF(){

	float muFactor, dxMu, dyMu, dzMu;

	float dxRhoX, dxRhoY, dxRhoZ;
	float dyRhoX, dyRhoY, dyRhoZ;
	float dzRhoX, dzRhoY ,dzRhoZ;		

	for(int x=1; x<XDIM-1; x++){
		for(int y=1; y<YDIM-1; y++){
			for(int z=1; z<ZDIM-1; z++){

				muFactor = - 1./(mu[x][y][z] * mu[x][y][z]);
			
				dxMu = mu[x-1][y][z] - mu[x][y][z];
				dyMu = mu[x][y-1][z] - mu[x][y][z];
				dzMu = mu[x][y][z-1] - mu[x][y][z];

				dxRhoX = rho[0][x-1][y][z] - rho[0][x][y][z];
				dxRhoY = rho[1][x-1][y][z] - rho[1][x][y][z];
				dxRhoZ = rho[2][x-1][y][z] - rho[2][x][y][z];

				dyRhoX = rho[0][x][y-1][z] - rho[0][x][y][z];
				dyRhoY = rho[1][x][y-1][z] - rho[1][x][y][z];
				dyRhoZ = rho[2][x][y-1][z] - rho[2][x][y][z];

				dzRhoX = rho[0][x][y][z-1] - rho[0][x][y][z];
				dzRhoY = rho[1][x][y][z-1] - rho[1][x][y][z];
				dzRhoZ = rho[2][x][y][z-1] - rho[2][x][y][z];	

				JInduced[0][x][y][z] = muFactor * ((dxMu * dxRhoX) + (dyMu * dxRhoY) + (dzMu * dxRhoZ));
				JInduced[1][x][y][z] = muFactor * ((dxMu * dyRhoX) + (dyMu * dyRhoY) + (dzMu * dyRhoZ));
				JInduced[2][x][y][z] = muFactor * ((dxMu * dzRhoX) + (dyMu * dzRhoY) + (dzMu * dzRhoZ));			
			}
		}
	}
}


void CollisionThreaded(){

	int dz = (ZDIM / cpuCount);		//Each slice of work the threads will do.

	for(int s=0; s<3; s++){		//Hits all of the simulations.

		pthread_t pth[cpuCount];	// Thread Identifier

		for(int threads = 0; threads < cpuCount; threads++){

			struct volume * threadBox = malloc(sizeof(*threadBox));

			threadBox->simulation = s;
			threadBox->xLow = 0;
			threadBox->xHigh = XDIM;
			threadBox->yLow = 0;
			threadBox->yHigh = YDIM;		

			threadBox->zLow = threads * dz;
			threadBox->zHigh = (threads + 1) * dz;
			if(threads + 1 == cpuCount){
				threadBox->zHigh = ZDIM;
			}
			
			pthread_create(&pth[threads],NULL,Collision,(void *) threadBox);		//Make thread for a specific region.
		}

		for(int threads = 0; threads < cpuCount; threads++){
			if(pthread_join(pth[threads],NULL) > 0){
				printf("Threads aren't joining correctly.\n");
			}
		}
	}
}

void DensityThreaded(){

	int dz = (ZDIM / cpuCount);		//Each slice of work the threads will do.

	for(int s=0; s<3; s++){		//Hits all of the simulations.

		pthread_t pth[cpuCount];	// Thread Identifier

		for(int threads = 0; threads < cpuCount; threads++){

			struct volume * threadBox = malloc(sizeof(*threadBox));

			threadBox->simulation = s;
			threadBox->xLow = 0;
			threadBox->xHigh = XDIM;
			threadBox->yLow = 0;
			threadBox->yHigh = YDIM;		

			threadBox->zLow = threads * dz;
			threadBox->zHigh = (threads + 1) * dz;
			if(threads + 1 == cpuCount){
				threadBox->zHigh = ZDIM;
			}
			
			pthread_create(&pth[threads],NULL,Density,(void *) threadBox);		//Make thread for a specific region.
		}

		for(int threads = 0; threads < cpuCount; threads++){
			if(pthread_join(pth[threads],NULL) > 0){
				printf("Threads aren't joining correctly.\n");
			}
		}
	}
}


void VisualizeThreaded(int s){

	int dz = (ZDIM / cpuCount);		//Each slice of work the threads will do.
	pthread_t pth[cpuCount];	// Thread Identifier

	for(int threads = 0; threads < cpuCount; threads++){

		struct volume * threadBox = malloc(sizeof(*threadBox));

		threadBox->simulation = s;	//Only one simulation can be viewed
		threadBox->xLow = 0;
		threadBox->xHigh = XDIM;
		threadBox->yLow = 0;
		threadBox->yHigh = YDIM;		

		threadBox->zLow = threads * dz;
		threadBox->zHigh = (threads + 1) * dz;
		if(threads + 1 == cpuCount){
			threadBox->zHigh = ZDIM;
		}
		
		pthread_create(&pth[threads],NULL,Visualize,(void *) threadBox);		//Make thread for a specific region.
	}

	for(int threads = 0; threads < cpuCount; threads++){
		if(pthread_join(pth[threads],NULL) > 0){
			printf("Threads aren't joining correctly.\n");
		}
	}

}


void Iterate(){

	DensityThreaded();
	CollisionThreaded();

	Stream(0);
	Stream(1);
	Stream(2);		
	Curl();
	GetFieldEnergy();
	CalcInducedEMF();
	
	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	//DeltaTime();
	
}

//-----Main Routines-----
int main(){
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=1;
  	int drawCount=0;


	DefineGraphNxN_R("XY Density",&rhoDisplay[0][0],&xDim,&yDim,NULL);  
	DefineGraphNxN_R("Theoretical Density",&rhoTh[0][0],&xDim,&yDim,NULL);  	
	DefineGraphNxN_RxR("B Field(X,Y)",&BDisplay[0][0][0],&xDim,&yDim,NULL);
	DefineGraphNxN_R("B Field(Z)",&BZDisplay[0][0],&xDim,&yDim,NULL);
	DefineGraphNxN_R("YZ Density",&rhoYZDisplay[0][0],&zDim,&yDim,NULL);
	DefineGraphNxN_RxR("B Field(Y,Z)",&BYZDisplay[0][0][0],&zDim,&yDim,NULL);	
	DefineGraphN_R("Rho(x)",&slice[0],&xDim,NULL);
	DefineGraphN_R("RhoTh(x)",&sliceTh[0],&xDim,&sliceReq);
	DefineGraphNxN_RxR("EMFDisplay",&EMFDisplay[0][0][0],&xDim,&yDim,NULL);	
	DefineGraphNxN_R("MuSlice",&muSlice[0][0],&xDim,&yDim,NULL);
	DefineGraphN_RxR("BH Curve",&BHCurve[0],&sizeX,NULL);
	DefineGraphN_RxR("H-Ur Curve",&HUrCurve[0],&sizeX,NULL);


	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineFloat("Tau", &T);
    DefineInt("Simulation", &Simulation);
    DefineInt("Slice", &sliceLoc);
    DefineBool("GetSlice",&rhoThreq);
    DefineFloat("J", &JMagnitude);
    DefineBool("Source", &sourceOn);
    DefineInt("Z plane", &zSlice);
    DefineInt("Threads", &cpuCount);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");    
  	DefineBool("Quit",&done);
  	EndMenu();

	Initialize();

  	while (!done){
  		if(drawCount>repeat){
      		Events(1);
      		VisualizeThreaded(Simulation);
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