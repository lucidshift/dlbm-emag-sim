#ifndef _DRAW3DGRAPH_H
#define _DRAW3DGRAPH_H

extern void draw3dgraph(int xsize,int ysize,grstr graph,double **gp, int *getgraph,
		 int nox, int noy,
		 int *newdim,int newdata, int comments, long no_iterations,
		 double persp_factor, double scalefakt, 
		 double xangle,double yangle,double zangle,
		 double xmax,double xmin,double ymax,double ymin,
		 double zmin,double zmax,double zfakt,
		 Meshp **mep,int *meinit);

extern void draw3dcurve(int xsize,int ysize,grstr graph,double **gp, int *getgraph,
	        int anzN_RxRxR, int comments,
		double persp_factor, double scalefakt, 
		double xangle,double yangle,
		double xmax,double xmin,int logx,double ymax,double ymin,int logy,
		double zmin,double zmax,int logz,double xfakt, double yfakt, double zfakt,
			char *Xlabel, char *Ylabel, char *Zlabel, grdat *mygrdat);
#endif
