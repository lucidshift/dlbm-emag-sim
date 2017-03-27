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
void Initialize();
void Iterate();
void *AnalyticalSolution(void *threadBox);
void AnalyticalSolutionThreaded(int s);
void GenerateSourcePath(int cx, int cy, int cz, int radius, int thickness); //Coil center location (x,y,z), radius, and wire thickness
void ConfigurePermeability(double a);
void GetFieldEnergy();
void CalcOrbitalCurrent();
void CalcHField();
void CalcMu();
void NetForceSimulation();
void ConfigureFreeSpace();
void ConfigureProjectile(int midpoint, int length, int radius, int permeability);
int SolutionConvergence();


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






