#ifndef _MYDRAW_2D_H
#define _MYDRAW_2D_H

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#define doublebufferX 1
extern int schwarz();
extern int weiss();
extern int rot();
extern int gruen();
extern int blau();
extern int gelb();
extern void set_schwarz(int d);
extern void set_weiss(int d);
extern void set_rot(int d);
extern void set_gruen(int d);
extern void set_blau(int d);
extern void set_gelb(int d);

extern void initX(Display *md, Drawable mydr, Drawable mydoubledr, GC myGC, 
	   int mode, int xsize, int ysize);
extern void setXwindow(Display *md, Drawable mydr, Drawable mydoubledr, GC myGC, 
		int mode, int _xsize, int _ysize);
extern void setX();
extern void setdisplay(Display *md);

extern int initPS(char *filename, double xsize, double ysize, 
	    int PSletter_, int PSlandscape_,int *PSxsize, int *PSysize);
extern void setPS();


extern int ColorRampStart();
extern int ColorRampEnd();
extern int NumberOfColors();
extern void SetColorRampStart(int i);
extern void SetColorRampEnd  (int i);

extern void SetColorFieldStart(int i);
extern int ColorFieldStart();
extern void SetColorFieldX(int i);
extern int ColorFieldX();
extern void SetColorFieldY(int i);
extern int ColorFieldY();

extern int  GetFreeColors();
extern void setcolormap(XColor *map, int n);


extern void (*mydrawline)(int color, int x1, int y1, int x2,int y2);
extern void (*mypolygon)(int color, XPoint *points, int n);
extern void (*mypolygon_line)(int color, int color_line, XPoint *points, int n);
extern void (*myline_polygon)(int color_line, XPoint *points, int n);
extern void (*mytriangle_2)(int color, int color_line, XPoint *points);
extern void (*mycircle)(int color, int x, int y, int r);
extern void (*myfilledcircle)(int color, int x, int y, int r);
extern void (*densityfield)(int nox, int noy,double zmin, double zmax,
			    int  colmin, int colmax, double *f,
			    int xofs, int yofs, int xsize, int ysize,int hex);
extern void (*twodensityfield)(int nox, int noy,
			       double z1min, double z1max,double z2min, double z2max,
			       int colomin,int colx, int coly, double *f1, double *f2,
			       int xofs, int yofs, int xsize, int ysize);

extern void (*mytext)(int color,int x, int y,const char *text, int orient);
extern void (*myselectfont)(int font, char *testtext, double length);
extern void (*mysetfontdirection)(double alpha);

extern void (*myclear)();
extern void (*myshow)();

#endif /* _MYDRAW_2D_H */
