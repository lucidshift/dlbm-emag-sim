#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>


#define pi M_PI
#define MPERM (4 * pi * 0.0000001)
#define XDIM 80
#define YDIM 80
#define ZDIM 180
#define velCount 7
#define u0 MPERM
#define u  MPERM * 3119		//Technically this is relative permeability when u0 is anything other than 4pi * 10^-7 
#define SIZEX 8000


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
double H[XDIM][YDIM][ZDIM][3];
double rhoDisplay[XDIM][YDIM];
double rhoYZDisplay[ZDIM][YDIM];
double BDisplay[XDIM][YDIM][2];	//x,y (location), and 2D vector 
double BYZDisplay[ZDIM][XDIM][2];	//z,y (location), and 2D vector 
double BZDisplay[XDIM][YDIM];
double mu[XDIM][YDIM][ZDIM];

double leftBuffer[3][velCount][YDIM][ZDIM]; 	//Velocity 1
double rightBuffer[3][velCount][YDIM][ZDIM];	//Velocity 2
double topBuffer[3][velCount][XDIM][ZDIM]; 		//Velocity 3
double bottomBuffer[3][velCount][XDIM][ZDIM];	//Velocity 4
double frontBuffer[3][velCount][XDIM][YDIM]; 	//Velocity 5
double backBuffer[3][velCount][XDIM][YDIM];		//Velocity 6

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double slice[XDIM];
double sliceTh[XDIM];
double rhoTh[XDIM][YDIM];
int sliceLoc = 33;
double k = 1;
int Simulation = 1;
int rhoThreq=0;
int sliceReq=0;
int zSlice = 25;
int ySlice = 25;
int xSlice = 25;

struct current {double x; double y; double z;} J[XDIM][YDIM][ZDIM]; //Current vector at a given location
double JInduced[3][XDIM][YDIM][ZDIM];
double JMagnitude = .111111111;// * MPERM;
int sourceOn = 1;

double OrbitalDisplay[XDIM][YDIM][2];
struct volume {int xLow; int xHigh; int yLow; int yHigh; int zLow; int zHigh; int simulation;};

int sizeX = SIZEX;
double HUrCurve[SIZEX];	//Collection of x,y points
double BHCurve[SIZEX]; 
double muSlice[XDIM][YDIM];
double alpha = 1;
int changePermeability=0;
double uCenter = u;
double Energy[SIZEX];
double BMiddle;

//Variables for ferric materials
double Amag[3][velCount][XDIM][YDIM][ZDIM]; 	//Simulation,D3Q7, X,Y,Z 
double Bmag[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
double Hmag[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
double rhomag[3][XDIM][YDIM][ZDIM];
int ferricCore = 1;
double leftBuffermag[3][velCount][YDIM][ZDIM]; 	//Velocity 1
double rightBuffermag[3][velCount][YDIM][ZDIM];	//Velocity 2
double topBuffermag[3][velCount][XDIM][ZDIM]; 		//Velocity 3
double bottomBuffermag[3][velCount][XDIM][ZDIM];	//Velocity 4
double frontBuffermag[3][velCount][XDIM][YDIM]; 	//Velocity 5
double backBuffermag[3][velCount][XDIM][YDIM];		//Velocity 6
double Htotal[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
double HYZCoreDisplay[ZDIM][XDIM][2];	//z,y (location), and 2D vector
double Force[XDIM][YDIM][ZDIM][3]; //x,y,z (location), and 3D vector 
double Umag[XDIM][YDIM][ZDIM];
double ForceDisplay[ZDIM][XDIM][2];	//z,y (location), and 2D vector
double UmagTotal;
double ForceTotal[3];

double ForceGraph[ZDIM];
double EnergyDelta[20];
int EDeltaPointer = 0;
double EDeltaThreshold = .05; //Threshold factor for determining simulation convergence
double EDeltaMax;
int coilMidpoint;
double ForceVsPosition[ZDIM][3];
double ZForceDisplay[ZDIM];
int GlobalXYZ[3];


//Main routine variables
int done=0;
int cont=0;
int reset=0;
int repeat=20;
int drawCount=0;
int simExit = 0;
int netForceSim = 0;


//Function Declaration
void *Collision(void *threadBox);
void CollisionThreaded();
void Initialize();
void Stream(int s);
void Curl();
void Iterate();
void *AnalyticalSolution(void *threadBox);
void AnalyticalSolutionThreaded(int s);
void GenerateSourcePath(int cx, int cy, int cz, int radius, int thickness); //Coil center location (x,y,z), radius, and wire thickness
void *Density(void *threadBox);
void DensityThreaded();
void *Visualize(void *threadBox);
void VisualizeThreaded(int s);
void ConfigurePermeability(double a);
void GetFieldEnergy();
void CalcOrbitalCurrent();
void CalcHField();
void CalcMu();
void NetForceSimulation();
void ConfigureFreeSpace();
void ConfigureProjectile(int midpoint, int length, int radius, int permeability);
int SolutionConvergence();




void Initialize(){

	//Sets up u0
	ConfigureFreeSpace();

	//Coil Setup - This section needs work.
	for(int z=57; z<(123); z=z+3){		//Fixed number of coils!!
		int r=11;
		//for(int r=11; r<25; r=r+3){

			coilMidpoint = ((123 - 57) / 2) + 57;	//Change this when you change the coil dimensions!
			printf("(%d,%d)\n",z,r);
			GenerateSourcePath(40, 40, z, r, 3);
		//}
	}

	//Current setup
	for(int x=0; x<(XDIM); x++){
		for(int y=0; y<(YDIM); y++){
			for(int z=0; z<(ZDIM); z++){
												
				Source[0][x][y][z] = JMagnitude * J[x][y][z].x * MPERM;
				Source[1][x][y][z] = JMagnitude * J[x][y][z].y * MPERM;
				Source[2][x][y][z] = JMagnitude * J[x][y][z].z * MPERM;
			}
		}
	}	

	if(netForceSim == 0){	
		//Buffer Initializations
		for(int s=0; s<3; s++){
			for(int x=0; x<(XDIM); x++){
				for(int y=0; y<(YDIM); y++){
					for(int z=0; z<(ZDIM); z++){
						for(int i=0; i<velCount; i++){

							A[s][i][x][y][z] = 0; //velDist[i] + Source[i][x][y];

							leftBuffer[s][i][y][z] = 0; 	//Velocity 1
							rightBuffer[s][i][y][z] = 0;	//Velocity 2
							topBuffer[s][i][x][z] = 0; 		//Velocity 3
							bottomBuffer[s][i][x][z] = 0;	//Velocity 4
							frontBuffer[s][i][x][y] = 0; 	//Velocity 5
							backBuffer[s][i][x][y] = 0;		//Velocity 6					
						}
					}
				}
			}	
		}	

		//Density Initialization
		for(int s=0; s<3; s++){
		    for (int x=0;x<(XDIM);x++){
		      	for (int y=0;y<(YDIM);y++){
					for(int z=0; z<(ZDIM); z++){
				      		rho[s][x][y][z] = 0;
					}
				}
		  	}
		}
	}

		for(int s=0; s<3; s++){
			for(int x=0; x<(XDIM); x++){
				for(int y=0; y<(YDIM); y++){
					for(int z=0; z<(ZDIM); z++){

						rhomag[s][x][y][z] = 0;
						for(int i=0; i<velCount; i++){

							Amag[s][i][x][y][z] = 0; //velDist[i] + Source[i][x][y];

							leftBuffermag[s][i][y][z] = 0; 	//Velocity 1
							rightBuffermag[s][i][y][z] = 0;	//Velocity 2
							topBuffermag[s][i][x][z] = 0; 		//Velocity 3
							bottomBuffermag[s][i][x][z] = 0;	//Velocity 4
							frontBuffermag[s][i][x][y] = 0; 	//Velocity 5
							backBuffermag[s][i][x][y] = 0;		//Velocity 6					
						}
					}
				}
			}	
		}

  	//Graph Initialization
	for(int x=0; x<(XDIM); x++){

		BHCurve[x] = 0;
		HUrCurve[x] = 0;		
	}

	//Reset to prevent initialization loop on false solution convergence.
	for(int i=0; i<20; i++){
		EnergyDelta[i] = 0;
	}
	EDeltaPointer = 0;
	EDeltaMax = 0;

	//Initialize Force arrays
    for (int x=0;x<(XDIM);x++){
      	for (int y=0;y<(YDIM);y++){
			for(int z=0; z<(ZDIM); z++){

				Force[x][y][z][0] = 0;
				Force[x][y][z][1] = 0;
				Force[x][y][z][2] = 0;
			}
		}
  	}

	//System variable setup
	iterationCount = 0;
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

void ConfigurePermeability(double a){

	int cx = 40;
	int cy = 40;
	//int cz = 25;
	int radius = 9;

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
		for(double r=0; r<(radius); r+=.5){
			for(int z=25; z<ZDIM-25; z++){

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

void ConfigureFreeSpace(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				mu[x][y][z] = u0;
			}
		}
	}
}

void ConfigureProjectile(int midpoint, int length, int radius, int permeability){


	int cx = 40;
	int cy = 40;
	//int cz = 25;

	int dx;
	int dy;		
	int xPos;
	int yPos;

	if(midpoint-(length/2) < 0){
		printf("Projectile out of lower simulation bounds.\n");
	}

	if(midpoint+(length/2) >= ZDIM){
		printf("Projectile out of upper simulation bounds.\n");
	}	
	
	for(double theta=0; theta<(2*pi); theta+=.1){
		for(double r=0; r<(radius); r+=.5){
			for(int z=midpoint-(length/2); z<midpoint+(length/2); z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				mu[xPos][yPos][z] = permeability;

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

					//Finew = Fi + (Fi0 - Fi) + (Source)
					A[s][i][x][y][z] = A[s][i][x][y][z] 
										+ ((velDist[i] * rho[s][x][y][z]) - A[s][i][x][y][z]) 
										+ (Source[s][x][y][z] * velDist[i] * sourceOn);					
				}
			}
		}
	}


	if(ferricCore == 1){

		for(int x= v.xLow; x< v.xHigh; x++){
			for(int y= v.yLow; y< v.yHigh; y++){
				for(int z= v.zLow; z< v.zHigh; z++){
					for(int i=0; i<velCount; i++){

					//Finew = Fi + (Fi0 - Fi) + (Source)
					Amag[s][i][x][y][z] = Amag[s][i][x][y][z] 
									+ ((velDist[i] * rhomag[s][x][y][z]) - Amag[s][i][x][y][z]) 
									+ (JInduced[s][x][y][z] * velDist[i]);
						
					}
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
			frontBuffer[s][5][x][y] = A[s][5][x][y][ZDIM-1];
		}
	}
	memmove(&A[s][5][0][0][1], &A[s][5][0][0][0], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));

	//Velocity 6
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			backBuffer[s][6][x][y] = A[s][6][x][y][0];
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


	if(ferricCore == 1){

		//Velocity 1
		for(int y=YDIM-1; y>=0; y--){
			for(int z=0; z<ZDIM; z++){
				leftBuffermag[s][1][y][z] = Amag[s][1][XDIM-1][y][z];
			}
		}
		memmove(&Amag[s][1][1][0][0], &Amag[s][1][0][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));	

		//Velocity 2
		for(int x=XDIM-1; x>=0; x--){
			for(int z=0; z<ZDIM; z++){
				bottomBuffermag[s][2][x][z] = Amag[s][2][x][YDIM-1][z];
			}	
		}
		memmove(&Amag[s][2][0][1][0], &Amag[s][2][0][0][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

		//Velocity 3
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				rightBuffermag[s][3][y][z] = Amag[s][3][0][y][z];
			}
		}
		memmove(&Amag[s][3][0][0][0], &Amag[s][3][1][0][0], (sizeof(double)*(XDIM-1)*(YDIM)*(ZDIM)));

		//Velocity 4
		for(int x=0; x<XDIM; x++){
			for(int z=0; z<ZDIM; z++){
				topBuffermag[s][4][x][z] = Amag[s][4][x][0][z];
			}
		}
		memmove(&Amag[s][4][0][0][0], &Amag[s][4][0][1][0], (sizeof(double)*(XDIM)*(YDIM-1)*(ZDIM)));

		//Velocity 5
		for(int x=0; x<XDIM; x++){
			for(int y=0; y<YDIM; y++){
				frontBuffermag[s][5][x][y] = Amag[s][5][x][y][ZDIM-1];
			}
		}
		memmove(&Amag[s][5][0][0][1], &Amag[s][5][0][0][0], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));

		//Velocity 6
		for(int x=0; x<XDIM; x++){
			for(int y=0; y<YDIM; y++){
				backBuffermag[s][6][x][y] = Amag[s][6][x][y][0];
			}
		}
		memmove(&Amag[s][6][0][0][0], &Amag[s][6][0][0][1], (sizeof(double)*(XDIM)*(YDIM)*(ZDIM-1)));	


		//Merge velocities from buffers
		for(int x=0; x<XDIM; x++){
			for(int z=0; z<ZDIM; z++){
				Amag[s][2][x][0][z] = bottomBuffermag[s][2][x][z];
				Amag[s][4][x][YDIM-1][z] = topBuffermag[s][4][x][z];
			}
		}

		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				Amag[s][1][0][y][z] = leftBuffermag[s][1][y][z];
				Amag[s][3][XDIM-1][y][z] = rightBuffermag[s][3][y][z];
			}
		}

		for(int x=0; x<XDIM; x++){
			for(int y=0; y<YDIM; y++){
				Amag[s][5][x][y][0] = frontBuffermag[s][5][x][y];
				Amag[s][6][x][y][ZDIM-1] = backBuffermag[s][6][x][y];
			}
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

    if(ferricCore == 1){

		for(int x= v.xLow; x< v.xHigh; x++){
			for(int y= v.yLow; y< v.yHigh; y++){
				for(int z= v.zLow; z< v.zHigh; z++){
		      		rhomag[s][x][y][z] = 0;
					for(int i=0; i<velCount; i++){
						rhomag[s][x][y][z] += Amag[s][i][x][y][z];
					}
				}
			}
	  	}  

    }  	

  	pthread_exit(0);
}

void Curl(){

	double B2, H2, B2mag, H2mag;
	double muTemp;

    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){
			for(int z=1; z<(ZDIM-1); z++){

				muTemp = u;
				
				//B field is calculated and redimensionalized here.
      			B[x][y][z][0] = (4 * pi * ((rho[2][x][y+1][z] - rho[2][x][y-1][z]) - (rho[1][x][y][z+1] - rho[1][x][y][z-1]))) / 2;
      			B[x][y][z][1] = (4 * pi * ((rho[0][x][y][z+1] - rho[0][x][y][z-1]) - (rho[2][x+1][y][z] - rho[2][x-1][y][z]))) / 2;
      			B[x][y][z][2] = (4 * pi * ((rho[1][x+1][y][z] - rho[1][x-1][y][z]) - (rho[0][x][y+1][z] - rho[0][x][y-1][z]))) / 2;

				B2 = sqrt(((B[x][y][z][0]*B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2])));	

				H[x][y][z][0] = B[x][y][z][0] / MPERM;
				H[x][y][z][1] = B[x][y][z][1] / MPERM;
				H[x][y][z][2] = B[x][y][z][2] / MPERM;

				H2 = sqrt((H[x][y][z][0]*H[x][y][z][0]) + (H[x][y][z][1]*H[x][y][z][1]) + (H[x][y][z][2]*H[x][y][z][2]));


				//This is the H field, because the u0 factor was moved to the other side, multiplied into the JInduced current. Divided by two due to how the difference is calculated.
      			Hmag[x][y][z][0] = (4 * pi * ((rhomag[2][x][y+1][z] - rhomag[2][x][y-1][z]) - (rhomag[1][x][y][z+1] - rhomag[1][x][y][z-1]))) / 2;
      			Hmag[x][y][z][1] = (4 * pi * ((rhomag[0][x][y][z+1] - rhomag[0][x][y][z-1]) - (rhomag[2][x+1][y][z] - rhomag[2][x-1][y][z]))) / 2;
      			Hmag[x][y][z][2] = (4 * pi * ((rhomag[1][x+1][y][z] - rhomag[1][x-1][y][z]) - (rhomag[0][x][y+1][z] - rhomag[0][x][y-1][z]))) / 2;

      			H2mag = sqrt(((Hmag[x][y][z][0]*Hmag[x][y][z][0]) + (Hmag[x][y][z][1]*Hmag[x][y][z][1]) + (Hmag[x][y][z][2]*Hmag[x][y][z][2])));

      			Bmag[x][y][z][0] = Hmag[x][y][z][0] / mu[x][y][z];
      			Bmag[x][y][z][1] = Hmag[x][y][z][1] / mu[x][y][z];
      			Bmag[x][y][z][2] = Hmag[x][y][z][2] / mu[x][y][z];

				B2mag = sqrt(((Bmag[x][y][z][0]*Bmag[x][y][z][0]) + (Bmag[x][y][z][1]*Bmag[x][y][z][1]) + (Bmag[x][y][z][2]*Bmag[x][y][z][2])));


				/*
				//Permeability Curve
				if(B2<0.0 && muTemp != u0){		//Don't change the permeability of free space, it doesn't change.
					mu[x][y][z] = u;
				}
				if(B2>0.0 && B2<1.5 && muTemp != u0){
					mu[x][y][z] = (3484.714 * (B2*B2*B2*B2) - 6261.368 * (B2*B2*B2) - 9598.177 * (B2*B2) + 15953.587 * B2 + 3119.105);
				}	
				if(B2>=1.5 && muTemp != u0){
					mu[x][y][z] = 608;
				}	
				*/								

				//H = B/u0 - M
				Htotal[x][y][z][0] = H[x][y][z][0] - Hmag[x][y][z][0];
				Htotal[x][y][z][1] = H[x][y][z][1] - Hmag[x][y][z][1];
				Htotal[x][y][z][2] = H[x][y][z][2] - Hmag[x][y][z][2];


				if(x==27 && y==27 && z==50){
					//printf("|B|=%10.10f\n", B2);
					//printf("Mu=%10.10f\n",muTemp);

					BMiddle = B2;

					if(((iterationCount*2)+1)<sizeX){
						BHCurve[iterationCount*2] = H2;	
						BHCurve[iterationCount*2 + 1] = B2;
						HUrCurve[iterationCount*2] = H2;
						HUrCurve[iterationCount*2 + 1] = muTemp;
						printf("H2=%5.25f,B2=%5.25f\n", H2, B2);
						printf("H2mag=%5.25f,B2mag=%5.25f\n", H2mag, B2mag);
					}					
				}	

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

void *AnalyticalSolution(void *threadBox){

	struct volume v = *((struct volume *) threadBox);
    free(threadBox);
    int s = v.simulation;	

	double distance = 0;
	int xOff;
	int yOff;
 	int zOff;

	for(int x=v.xLow; x<v.xHigh; x++){

		if(v.xHigh == XDIM){
			printf("Iterations Remaining = %d\n",XDIM - x);
		}	

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

									if(x!=xx || y!=yy){																						//Scale Factor
										distance = sqrt(((x-xx+xOff)*(x-xx+xOff)) + ((y-yy+yOff)*(y-yy+yOff)) + ((z-zz+zOff)*(z-zz+zOff))) * 0.001;
									}
											
									//Does not currently have alpha factored in it, seems to break when that's the case.					
									rhoTh[x][y] += ((MPERM) / (4*pi)) * ((Source[s][xx][yy][zz] / MPERM) / distance);    		//Source has MPERM in it.
								}
							}
						}
					}
				}
			}
		//}	
	}

	pthread_exit(0);
}



void GetGraph(int s){

	if(rhoThreq && sliceReq){ //ANDed because otherwise it constantly calculates
		AnalyticalSolutionThreaded(s);
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


int SolutionConvergence(){

	double energyDeltaTemp = Energy[iterationCount] - Energy[iterationCount - 1];

	if(iterationCount < 19)
	{
		return 0; //Too Soon :)
	} 
	else if(netForceSim == 1)
	{
		return 1;
	}

	if(iterationCount < SIZEX){
		EnergyDelta[iterationCount] = energyDeltaTemp;
	}

	if(energyDeltaTemp > EDeltaMax){
		EDeltaMax = energyDeltaTemp;
	}

	if(energyDeltaTemp < (EDeltaMax * EDeltaThreshold)){
		printf("Solution Convergence.\n");
		return 1;
	}	
	
	return 0;
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

				OrbitalDisplay[x][y][0] = JInduced[0][x][y][zSlice];
				OrbitalDisplay[x][y][1] = JInduced[1][x][y][zSlice];

				muSlice[x][y] = mu[x][y][zSlice];
				
				//Projection into 2D
				BDisplay[x][y][0] = B[x][y][zSlice][0];
				BDisplay[x][y][1] = B[x][y][zSlice][1];
				BYZDisplay[z][y][0] = B[xSlice][y][z][2]; //
				BZDisplay[x][y] = B[x][y][zSlice][2];
				BYZDisplay[z][y][1] = B[xSlice][y][z][1]; //

				HYZCoreDisplay[z][y][0] = Hmag[xSlice][y][z][2];
				HYZCoreDisplay[z][y][1] = Hmag[xSlice][y][z][1];

				ForceDisplay[z][y][0] = Force[xSlice][y][z][2];
				ForceDisplay[z][y][1] = Force[xSlice][y][z][1];

			}
		}
	}

	pthread_exit(0);
}


void GetFieldEnergy(){

	double eTotal = 0;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
			
							// 		sum(B^2 * dv * 1/(2u0))															//dv is 1mm, hence the factor of 10^-9 for the B field
				eTotal += (((B[x][y][z][0] * B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2])) * 0.000000001) / (2 * MPERM);


			}
		}
	}
	Energy[iterationCount] = eTotal;
	printf("Total Field Energy = %10.15f\n", eTotal);
}


void CalcOrbitalCurrent(){

	//double muFactor; 
	double dxMu, dyMu, dzMu;
	//double badMu = - 1./(u0 * u0 * MPERM * MPERM);

	double dxRhoX, dxRhoY, dxRhoZ;
	double dyRhoX, dyRhoY, dyRhoZ;
	double dzRhoX, dzRhoY, dzRhoZ;		

	for(int x=1; x<XDIM-1; x++){
		for(int y=1; y<YDIM-1; y++){
			for(int z=1; z<ZDIM-1; z++){

				//muFactor = (mu[x][y][z] * mu[x][y][z] * MPERM * MPERM) / (1000 * 0.148809);

				/*
				if(iterationCount > 50 && x == 26){
					printf("muFactor=%10.10f\n",muFactor);
				}				
				*/
			
				dxMu = (mu[x-1][y][z] - mu[x+1][y][z]) / 2;
				dyMu = (mu[x][y-1][z] - mu[x][y+1][z]) / 2;
				dzMu = (mu[x][y][z-1] - mu[x][y][z+1]) / 2;

				/*
				if(dxMu > 0 && muFactor == badMu){
					muFactor = -1;//- 1./(mu[x-1][y][z] * mu[x-1][y][z]);
				}

				if(dyMu > 0 && muFactor == badMu){
					muFactor = -1;//- 1./(mu[x][y-1][z] * mu[x][y-1][z]);
				}

				if(dzMu > 0 && muFactor == badMu){
					muFactor = -1;//- 1./(mu[x][y][z-1] * mu[x][y][z-1]);
				}
				*/

				dxRhoX = (rho[0][x-1][y][z] - rho[0][x+1][y][z]) / 2;
				dxRhoY = (rho[1][x-1][y][z] - rho[1][x+1][y][z]) / 2;
				dxRhoZ = (rho[2][x-1][y][z] - rho[2][x+1][y][z]) / 2;

				dyRhoX = (rho[0][x][y-1][z] - rho[0][x][y+1][z]) / 2;
				dyRhoY = (rho[1][x][y-1][z] - rho[1][x][y+1][z]) / 2;
				dyRhoZ = (rho[2][x][y-1][z] - rho[2][x][y+1][z]) / 2;

				dzRhoX = (rho[0][x][y][z-1] - rho[0][x][y][z+1]) / 2;
				dzRhoY = (rho[1][x][y][z-1] - rho[1][x][y][z+1]) / 2;
				dzRhoZ = (rho[2][x][y][z-1] - rho[2][x][y][z+1]) / 2;	

				//Modified to test uJ hypothesis
				JInduced[0][x][y][z] = -(mu[x][y][z]/u0) * ((dxMu * dxRhoX) + (dyMu * dxRhoY) + (dzMu * dxRhoZ));// * muFactor;
				JInduced[1][x][y][z] = -(mu[x][y][z]/u0) * ((dxMu * dyRhoX) + (dyMu * dyRhoY) + (dzMu * dyRhoZ));// * muFactor;
				JInduced[2][x][y][z] = -(mu[x][y][z]/u0) * ((dxMu * dzRhoX) + (dyMu * dzRhoY) + (dzMu * dzRhoZ));// * muFactor;	

				/*
				if(iterationCount > 50 && x == 26){
					printf("(%10.10f, %10.10f, %10.10f)\n", JInduced[0][x][y][z], JInduced[1][x][y][z], JInduced[2][x][y][z]);
				}
				*/
						
			}
		}
	}
}


void CalcForce(){

	double uTemp = 0;
	double forceTemp[3] = {0,0,0};

	//Potential Energy Calculation
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){

				if(mu[x][y][z] != u0){
					Umag[x][y][z] = ((Hmag[x][y][z][0]*B[x][y][z][0]) + (Hmag[x][y][z][1]*B[x][y][z][1]) + (Hmag[x][y][z][2]*B[x][y][z][2]));
					uTemp += Umag[x][y][z];
				}
			}
		}
	}

	//Force Calculation
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){

				if(mu[x][y][z] != u0){
					Force[x][y][z][0] = (Umag[x+1][y][z] - Umag[x-1][y][z]) / 2;
					Force[x][y][z][1] = (Umag[x][y+1][z] - Umag[x][y-1][z]) / 2;
					Force[x][y][z][2] = (Umag[x][y][z+1] - Umag[x][y][z-1]) / 2;

					forceTemp[0] += Force[x][y][z][0];
					forceTemp[1] += Force[x][y][z][1];
					forceTemp[2] += Force[x][y][z][2];
				}
			}
		}
	}	

	UmagTotal = uTemp;
	ForceTotal[0] = forceTemp[0];
	ForceTotal[1] = forceTemp[1];
	ForceTotal[2] = forceTemp[2];

	//For multiposition simulations
	ForceVsPosition[GlobalXYZ[2]][0] = forceTemp[0];
	ForceVsPosition[GlobalXYZ[2]][1] = forceTemp[1];
	ForceVsPosition[GlobalXYZ[2]][2] = forceTemp[2];
	ZForceDisplay[GlobalXYZ[2]] = forceTemp[2];

	printf("UmagTotal=%5.15f\n", UmagTotal);
	printf("ForceTotal=(%5.15f, %5.15f, %5.15f)\n", forceTemp[0], forceTemp[1], forceTemp[2]);
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


void AnalyticalSolutionThreaded(int s){

	int dx = (XDIM / cpuCount);		//Each slice of work the threads will do.
	pthread_t pth[cpuCount];	// Thread Identifier

	for(int x=0; x<XDIM; x++){			//This only needs to be done once.
		for(int y=0; y<YDIM; y++){
		rhoTh[x][y] = 0;
		}
	}	

	for(int threads = 0; threads < cpuCount; threads++){

		struct volume * threadBox = malloc(sizeof(*threadBox));

		threadBox->simulation = s;	//Only one simulation can be viewed
		threadBox->zLow = 0;
		threadBox->zHigh = ZDIM;
		threadBox->yLow = 0;
		threadBox->yHigh = YDIM;		

		threadBox->xLow = threads * dx;
		threadBox->xHigh = (threads + 1) * dx;
		if(threads + 1 == cpuCount){
			threadBox->xHigh = XDIM;
		}
		
		pthread_create(&pth[threads],NULL,AnalyticalSolution,(void *) threadBox);		//Make thread for a specific region.
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
	CalcForce();
	GetFieldEnergy();
	CalcOrbitalCurrent();
	simExit = SolutionConvergence();

	iterationCount++;
	printf("Iteration: %d\n", iterationCount);
	//DeltaTime();
	
}

//-----Main Routines-----
int main(){
	DefineGraphNxN_R("XY Density",&rhoDisplay[0][0],&xDim,&yDim,NULL);  
	DefineGraphNxN_RxR("XY B-Field",&BDisplay[0][0][0],&xDim,&yDim,NULL);	
	DefineGraphNxN_R("Z B-Field",&BZDisplay[0][0],&xDim,&yDim,NULL);	
	DefineGraphNxN_RxR("Orbital Current",&OrbitalDisplay[0][0][0],&xDim,&yDim,NULL);	
	DefineGraphNxN_R("MuSlice",&muSlice[0][0],&xDim,&yDim,NULL);
	DefineGraphNxN_R("YZ Density",&rhoYZDisplay[0][0],&zDim,&yDim,NULL);
	DefineGraphNxN_RxR("YZ B-Field",&BYZDisplay[0][0][0],&zDim,&yDim,NULL);
	DefineGraphNxN_RxR("YZ H-Core",&HYZCoreDisplay[0][0][0],&zDim,&yDim,NULL);	
	DefineGraphNxN_RxR("YZ Force",&ForceDisplay[0][0][0],&zDim,&yDim,NULL);

	DefineGraphN_R("Rho(x)",&slice[0],&xDim,NULL);
	DefineGraphN_R("RhoTh(x)",&sliceTh[0],&xDim,&sliceReq);
	DefineGraphN_RxR("BH Curve",&BHCurve[0],&iterationCount,NULL);
	DefineGraphN_RxR("H-Ur Curve",&HUrCurve[0],&iterationCount,NULL);
	DefineGraphN_R("Energy",&Energy[0],&iterationCount,NULL);
	DefineGraphN_R("EnergyDelta",&EnergyDelta[0],&iterationCount,NULL);
	DefineGraphN_R("Force v Position",&ZForceDisplay[0],&zDim,NULL);


	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineDouble("Alpha", &alpha);
    DefineInt("Simulation", &Simulation);
    DefineInt("Slice", &sliceLoc);
    DefineBool("Analytical",&rhoThreq);
    DefineDouble("J", &JMagnitude);
    DefineDouble("K", &k);
    DefineBool("Source", &sourceOn);
    DefineInt("Z plane", &zSlice);
    DefineInt("Threads", &cpuCount);
    DefineBool("U", &changePermeability);
    DefineDouble("U", &uCenter);
    DefineDouble("BMiddle", &BMiddle);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");    
    DefineBool("Ferromagnetics",&ferricCore);
    DefineBool("Force Step", &simExit);
    DefineInt("Pojectile Pos", &GlobalXYZ[2]);
  	DefineBool("Quit",&done);
  	EndMenu();




	NetForceSimulation();
}


void NetForceSimulation(){

	//Setup
	int startPos = 12;	//Midpoint of the projectile
	int endPos = 0;		//Set later
	int length = 20;
	int radius = 4;


	GlobalXYZ[0] = 40;
	GlobalXYZ[1] = 40;
	GlobalXYZ[2] = startPos;

	Initialize();
	endPos = 168; //coilMidpoint + 5;

	for(int i=startPos; i<=endPos; i++){

		Initialize();
		ConfigureProjectile(i, length, radius, u);
		GlobalXYZ[2] = i;
		simExit = 0;
		drawCount = repeat + 1;
	
	  	while (simExit == 0){
	  		if(drawCount>repeat){
	      		Events(1);
	      		VisualizeThreaded(Simulation);
	      		GetGraph(Simulation);
	      		DrawGraphs();
	      		if(changePermeability){
	      			ConfigurePermeability(uCenter); //TODO: Remove/update this feature, it's broken if used.
	      			changePermeability = 0;
	      		}		

	      		drawCount = 0;
	  		}

		    if((cont == 1) && (simExit == 0)){
		    	Iterate();		    	
		    }

		    drawCount++;
	  	}

	  	netForceSim = 1; //Set after first iteration.
	  	printf("i=%d\n",i);
	}







	//Very end, holds program from exiting.
	while(!done){
		//Find pause command?
		Events(1);
	    VisualizeThreaded(Simulation);
	    GetGraph(Simulation);
	    DrawGraphs();
		DrawGraphs();
	}
}






