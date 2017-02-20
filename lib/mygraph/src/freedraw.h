#ifndef _freedraw_H 
#define _freedraw_H

typedef void (*FrDr)(int xsize,int ysize);
extern void AddFreedraw(char *name,FrDr FreeDraw);
extern FrDr FreeDraw(int i);
char * FreeDrawName(int i);
int FreeDrawNumber();

#endif
