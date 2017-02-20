#ifndef _UFIELD_H
#define _UFIELD_H

#include "basicdef.h"

extern void ufield(int anzy,int xmin, int xmax, int ymin, int ymax,
		   double *u,double scale, int xscale,
		   int xofs, int yofs, int xsize, int ysize,
		   int veccol1, int veccol2, int vecstep, double vecangle, double vecwidth);

#endif




