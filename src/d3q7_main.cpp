#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <vector.h>
#include <thread>

extern "C" {
	#include <mygraph.h>
}

#include "d3q7.h"
class d3q7;

#define itCountMax 30

int iterationCount = 0;
int cpuCount=0;
int cubeMax = 500;

int xDim = 50;
int yDim = 50;
int zDim = 50;

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

double rhoDisplayArray[50][50];
void * rhoDisplayPointer;
int sliceLoc = 33;

/*d3q7::DensityField3D _rho;
d3q7::DensityField3D _inputSource;
d3q7::DensityField2D _rhoDisplay ;*/
double timeAvg [itCountMax];
double timeAvgTotal;

void deltaTime(){

	//TODO: Fix the wrapping error.
	struct timespec spec;
	clock_gettime(CLOCK_REALTIME, &spec);

	t.past2 = t.past1;
	t.past1 = t.now;
	t.now = spec.tv_nsec / 1000000.0;
	t.dt = t.now - t.past1;
	//printf ("deltaTime = %2.5f ms\n", t.dt);

}

//-----Main Routines-----
int main(){
	printf("Main startup.\n");
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=100;
  	int drawCount=0;
  	int threadCount=4;

	DefineGraphNxN_R("XY Density",&rhoDisplayArray[0][0],&xDim,&yDim,NULL);

	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineInt("Slice", &sliceLoc);
    DefineInt("Cube Dim Max", &cubeMax);
    DefineInt("Thread Count", &threadCount);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");
  	DefineBool("Quit",&done);
  	EndMenu();

  	//Simulation initialization
	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores (and hyperthreads) availible = %d\n", cpuCount);

	//Hold at startup for configuration.
	while(cont==0)
	{
		usleep(100000);
		Events(1);
	}

	while (!done){
		for(int i=50; i<cubeMax; i++){
			timeAvgTotal = 0;
			xDim = i;
			yDim = i;
			zDim = i;

			std::vector<d3q7*> simulation;

			d3q7::DensityField3D _rho = new double [xDim * yDim * zDim];
			d3q7::DensityField3D _inputSource = new double [xDim * yDim * zDim];
			double (*rho)[yDim][zDim] = (double(*)[yDim][zDim]) _rho;
			double (*inputSource)[yDim][zDim] = (double(*)[yDim][zDim]) _inputSource;

			//Loop through many sizes
			for(int z=0; z<zDim; z++){
				inputSource[17][17][z] = 1;
				inputSource[34][17][z] = -1;
				inputSource[17][34][z] = -1;
				inputSource[34][34][z] = 1;
			}

			for(int i=0; i<threadCount; i++)
			{
				d3q7 * simTemp = new d3q7(xDim, yDim, zDim);
				simulation.push_back(simTemp);
				//d3q7::DensityField2D _rhoDisplay = simulation[i].getSlice(sliceLoc);
			}
			//printf("Finished allocating memory.\n");

			for(int i=0; i<threadCount; i++)
			{
				simulation.at(i)->loadSource(_inputSource);
				//d3q7::DensityField2D _rhoDisplay = simulation[i].getSlice(sliceLoc);
			}
			//printf("Simulation initialization complete.\n");

			std::vector<std::thread> threads;
			deltaTime(); //Set the count for the iteration.
			for(int it=0; it<itCountMax; it++){

		  		if(drawCount>repeat){
		      		Events(1);

		      		//rhoDisplay* = simulation.getSlice();
		      		DrawGraphs();
		      		drawCount = 0;
		  		}

			    if(cont){

			    	for(int index=0; index<threadCount; index++)
			    	{
			    		threads.push_back( std::thread(&d3q7::iterate,  simulation.at(index)) );
/*			    		pthread_create(&pth[index],NULL,
			    			(void*) &(simulation[index].iterate),
			    			void);*/
			    	}

					for(int index = 0; index<threadCount; index++){
						threads[index].join();
					}

					// _rhoDisplay = simulation.getSlice(sliceLoc);
					// double (*rhoDisplay)[zDim] = (double(*)[zDim]) _rhoDisplay;

/*					for(int x=0; x<xDim; x++){
						for(int y=0; y<yDim; y++){

							rhoDisplayArray[x][y] = rhoDisplay[x][y];
						}
					}*/
					//rhoDisplayArray* = rhoDisplay*;	 //Could get dicey here.
			    	if(t.dt > 0)
			    	{
			    		timeAvg[it] = t.dt;
			    		timeAvgTotal += t.dt;

			    	}
			    	else
			    	{
			    		//Copy last iteration value so the average stays consistent
			    		timeAvg[it] = timeAvg[it-1];
			    		timeAvgTotal += timeAvg[it-1];
			    	}
				    deltaTime();	//TODO: Fix time wrapping, maybe rewrite whole thing.

				    threads.clear();
			    }
			    drawCount++;

				//Pause iteration.
				while(cont==0)
				{
					usleep(100000);
					Events(1);
				}
		  	}

			timeAvgTotal = timeAvgTotal / itCountMax;
			printf("Cube Size = %u, Total Volume = %u, Time delay = %7.7f ms, Elements per second = %7.0f\n", i, i*i*i*threadCount, timeAvgTotal, (i*i*i*threadCount)/timeAvgTotal );

	/*	  	delete [] _rho;
		  	delete [] _rhoDisplay;
		  	delete [] _inputSource;	 */
	  	}
	}
}