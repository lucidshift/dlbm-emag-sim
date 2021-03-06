#ifndef _BASICDEF_H
#define _BASICDEF_H

# include	<X11/Xlib.h>
# include	<X11/Xutil.h>

#include "mydraw.h"
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399
#endif

typedef long (*indfn)(int,int);
typedef long (*indfxfn)(int ,int, int);


typedef struct nz1data {
  int comments;
  int Autoscale;
  int scaleonce;
  int rough;
  int NEWDIMENSIONS;
  int newdata;
  int nox;
  int noy;
  int *getgraph;
  double persp_factor;
  double scalefakt;
  double xangle;
  double yangle;
  double zfact;
  double zmin;
  double zmax;
  int copy[2];
  Meshp *me;
  int meinit;
} nz1d;

typedef struct nz2data {
  int comments;
  int koord;
  double koordxmin;
  double koordxmax;
  int koordlogx;
  double koordymin;
  double koordymax;
  int koordlogy;
  char Xlabel[500],Ylabel[500];

  int density;
  int DensityLegend;
  int BarInset;
  double BarHeight;
  double BarWidth;
  double BarXOfs;
  double BarYOfs;
  int hex;
  double density_min;
  double density_max;
  double density1_min;
  double density1_max;
  double density2_min;
  double density2_max;
  int newdata;
  int rough;
  int copy[2];
  int copyy;
  int nox;
  int noy;
  int *draw;
  int nocuts;
  int nocuts1;
  int maxcuts;
  int maxcuts1;
  int adjustcuts;
  int adjustdensity;
  double *cut;
  double *cut1;
  int *draw_ufield;
  int *draw_tfield;
  int *rough_t;
  int NewXsize;
  int NewYsize;
  int DoCenter;
  double xofs;
  double yofs;
  int magnify;
  /* Sept 2010: new two-density plot data */
  int drawtwodensity,d1,d2;
} nz2d;

typedef struct project3d_data {
  int comments;
  int density;
  double density_min;
  double density_max;
  int newdata;
  int copy[2];
  int copyy;
  int nox;
  int noy;
  int *drawX,*drawY,*drawZ;
  int nocuts;
  int maxcuts;
  int adjustcuts;
  int adjustdensity;
  int *rough_t;
  int NewXsize;
  int NewYsize;
  int DoCenter;
  double xofs;
  double yofs;
  int magnify;
} p3dd;

typedef struct nzcontdata {  
  int comments;
  int rough;
  int newdata;
  int *getgraph;
  double persp_factor;
  double scalefakt;
  double xangle;
  double yangle;
  int copy[3];
} nzcontour3d;

typedef struct nzfrdat {  
  int newdata;
  int *draw;
} nzfreedata;


typedef struct time_ {
  int comments;
  int newdata;
  int restart;
} mytime;

typedef struct nzgraph2d_ {
  int comments;
  int newdata;
  int *draw;
  int *draw1;
  int *draw2; 
  int rough;
  int copy;
  int scaling;
  double xmin,xmax,ymin,ymax;
  int logx, logy;
  char Xlabel[500],Ylabel[500];
  double size;
} nzgraph2d;

typedef struct nzgraph3d_ {
  int comments;
  int newdata;
  int *draw;
  int rough;
  int copy;
  double persp_fakt;
  double scalefakt;
  double xfakt;
  double yfakt;
  double zfakt;
  double xangle;
  double yangle;
  int scaling;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  int logx, logy,logz;
  char Xlabel[500],Ylabel[500],Zlabel[500];
  double size;
} nzgraph3d;

typedef union _nzdata{
  nz1d nz1;
  nz2d nz2;
  nzfreedata nzfree;
  p3dd p3d;
  nzcontour3d nzc;
  mytime time;
  nzgraph2d graph2d;
  nzgraph3d graph3d;
} nzdata;

/* Part of a grand plan that can at the current state not work together
   with ZPL. A more primitive approach is needed. 

typedef void (*gr1func)(int,double *);
typedef void (*gr2func)(int,int,double *);
typedef void (*gr3func)(int,int,int,double *);

typedef union grfunc_{
  gr1func gr1;
  gr2func gr2;
  gr3func gr3;
} grfunc;
*/


#endif








