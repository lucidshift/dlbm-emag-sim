#include <math.h>
#include <time.h>
#include <mygraph.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "d3q7.h"

int iterationCount = 0;
int cpuCount=0;

int xDim = 50;
int yDim = 50;
int zDim = 50;

struct time {double now; double past1; double past2; double dt;};
struct timespec sDelay;
clock_t timeValue;
struct time t;

d3q7::DensityField2D rhoDisplay;
void * rhoDisplayPointer;
int sliceLoc = 33;

void deltaTime(){

	struct timespec spec;
	clock_gettime(CLOCK_REALTIME, &spec);

	t.past2 = t.past1;
	t.past1 = t.now;
	t.now = spec.tv_nsec / 1000000.0;
	t.dt = t.now - t.past1;
	printf ("deltaTime = %2.5f ms\n", t.dt);

}

//-----Main Routines-----
int main(){
  	int done=0;
  	int cont=0;
  	int reset=0;
  	int repeat=100;
  	int drawCount=0;

	DefineGraphNxN_R("XY Density",&rhoDisplay[0][0],&xDim,&yDim,NULL);  		

	StartMenu("Start",1);
    DefineBool("Run Simulation",&cont);
    DefineBool("Reset", &reset);
    DefineInt("Repeat", &repeat);
    DefineInt("Slice", &sliceLoc);
    DefineGraph(contour2d_,"Densities");
    DefineGraph(curve2d_,"Graphs");    
  	DefineBool("Quit",&done);
  	EndMenu();

  	//Simulation initialization
	d3q7::DensityField3D sourceTemp;
	sourceTemp.reserve(sizeof(double)*xDim*yDim*zDim);
	for(int z=0; z<zDim; z++){
		sourceTemp[17][17][z] = 1;
		sourceTemp[34][17][z] = -1;
		sourceTemp[17][34][z] = -1;	
		sourceTemp[34][34][z] = 1;	
	}

	d3q7 simulation = d3q7(xDim, yDim, zDim);
	simulation.loadSource(sourceTemp);

	//rhoDisplayTemp[xDim][yDim];

	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

  	while (!done){
  		if(drawCount>repeat){
      		Events(1);

      		//rhoDisplay* = simulation.getSlice();
      		DrawGraphs();		
      		drawCount = 0;
  		}
	    if(cont){
	    	simulation.iterate();
	    	deltaTime();
	    }

	    drawCount++;
  }
}