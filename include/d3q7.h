#ifndef D3Q7_H
#define D3Q7_H

#include <stdint.h>


void d3q7(int xDim, int yDim, int zDim);

bool iterate();

bool loadSource(double* source);

double* getArray();

double* getSlice(int zDimSlice);

//"Private" members

void collision();

void stream();

void density();


#endif