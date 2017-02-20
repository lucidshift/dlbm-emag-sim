#include <stdlib.h>
#include <math.h>
#include "color.h"
#include "mydraw_2d.h"


Colormap     cmap;

void (*ColRamp)(double ,double *,double *,double *);
void (*ColField)(double ,double ,double *,double *,double *);

void MyColRamp(double frac,double *red,double *green,double *blue){
  *red   = 0.991* (0.80*frac+0.2*frac);
  *green = 0.991* (0.80*frac+0.1); 
  *blue  = 0.991* (0.80*frac+0.2-0.2*frac);
}

void GreenRedColRamp(double frac,double *red,double *green,double *blue){
  *red   =  (1.0-frac);
  *green =  (0.80*frac+0.1); 
  *blue  =  (0.80*sin(M_PI*frac)+0.1);
}

void MyColField(double x, double y,double *red,double *green,double *blue){
  *red   = x;
  *green = y;
  *blue  = y*x;
}

void SetColorRamp(void (*NewColRamp)(double ,double *,double *,double *))
{
  unsigned int
    myncolors;
  XColor *mymap;
  int i,j,c;
  double red,green,blue;
  if (NewColRamp==NULL) ColRamp=MyColRamp;
  else ColRamp=NewColRamp;
  ColField=MyColField;

  myncolors=GetFreeColors();

  if (myncolors<6){
    printf("Not enough colors (%i) to allocate standard colors!",
	   myncolors);
  }
  SetColorRampStart(6);
  if (myncolors<256){
    printf("Not enough colors (%i) to allocate color ramp!/n",myncolors);
    SetColorRampEnd(myncolors-1);
  }
  else SetColorRampEnd(255);
  c=sqrt(myncolors-ColorRampEnd());
  SetColorFieldX(c);
  SetColorFieldY(c);
  if (myncolors<ColorRampEnd()+ColorFieldX()*ColorFieldY()){
    printf("Not enough colors (%i) to allocate color field!",
	   myncolors);
  }
  SetColorFieldStart(ColorRampEnd()+1);

  mymap = (XColor *) malloc(myncolors*sizeof(XColor));

  set_weiss(0);
  mymap [0].red   = 65000; 
  mymap [0].green = 65000;
  mymap [0].blue  = 65000;
  set_schwarz(1);
  mymap [1].red   = 0; 
  mymap [1].green = 0;
  mymap [1].blue  = 0;
  set_rot(2);
  mymap [2].red  = 65000; 
  mymap [2].green= 05000;
  mymap [2].blue = 0;
  set_gruen(3);
  mymap [3].red  = 0; 
  mymap [3].green= 65000;
  mymap [3].blue = 0;
  set_blau(4);
  mymap [4].red  = 0; 
  mymap [4].green= 20000;
  mymap [4].blue = 65000;
  set_gelb(5);
  mymap [5].red  = 65000; 
  mymap [5].green= 65000;
  mymap [5].blue = 30000;
  for (i=ColorRampStart(); i<=ColorRampEnd(); i++) {
    ColRamp(1.0*(i-ColorRampStart())/(ColorRampEnd()-ColorRampStart()),&red,&green,&blue);
    mymap [i].red   = 65535*red;
    mymap [i].green = 65535*green;
    mymap [i].blue  = 65535*blue;
  }
  for (i=0; i<ColorFieldX(); i++) {
    for (j=0; j<ColorFieldY(); j++) {
      ColField((1.0*i)/ColorFieldX(),(1.0*j)/ColorFieldY(),&red,&green,&blue);
      mymap [i*ColorFieldY()+j+ColorFieldStart()].red   = 65535*red;
      mymap [i*ColorFieldY()+j+ColorFieldStart()].green = 65535*green;
      mymap [i*ColorFieldY()+j+ColorFieldStart()].blue  = 65535*blue;
    }
  }
  setcolormap(mymap, myncolors);
}
