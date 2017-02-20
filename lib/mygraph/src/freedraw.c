#include <stdio.h>
#include <stdlib.h>
#include "freedraw.h"


static int NoFunct=0,MaxNoFunct=0;
static FrDr *FreeDrawArr=NULL;
static char **FreeDrawNameArr=NULL;

void AddFreedraw(char *name,FrDr FreeDraw){
  int i;
  NoFunct++;
  if (NoFunct>MaxNoFunct){
    MaxNoFunct+=10;
    FreeDrawArr= (FrDr *) realloc(FreeDrawArr,MaxNoFunct*sizeof(FrDr));
    FreeDrawNameArr= (char **) realloc(FreeDrawNameArr,MaxNoFunct*sizeof(char **));
    for (i=-10;i<0;i++) FreeDrawNameArr[MaxNoFunct+i]=(char *)malloc(256*sizeof(char));
  }
  sprintf(FreeDrawNameArr[NoFunct-1],"%s",name);
  FreeDrawArr[NoFunct-1]=FreeDraw;
}

FrDr FreeDraw(int i){
  if (i>NoFunct-1){
    printf("Error in FreeDraw, no function number %i\n",i);
    return NULL;
  }
  return FreeDrawArr[i];
}

char * FreeDrawName(int i){
 if (i>NoFunct-1){
    printf("Error in FreeDrawName, no function number %i\n",i);
    return NULL;
  }
  return FreeDrawNameArr[i];
}

int FreeDrawNumber(){
  return NoFunct;
}
