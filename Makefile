simple: 
	gcc -pg -Wall -I ~/include -L ~/lib src/LatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o output/simple2D

toroidsimple: 
	gcc -pg -Wall -I ~/include -L ~/lib src/ToroidLatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o output/toroidsimple2D

memcpysimple: 
	gcc -O3 -Wall -I ~/include -L ~/lib src/MemcpyLatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o output/memcpysimple2D

fourcharges: 
	gcc -O3  -I ~/include -L ~/lib src/FourCharges.c -std=gnu99 -lgraph -lm -lX11 -o output/fourcharges

plates2d: 
	gcc -O3 -I ~/include -L ~/lib src/Plates2D.c -std=gnu99 -lgraph -lm -lX11 -o output/plates2D	

tausweep: 
	gcc -g -Wall -I ~/include -L ~/lib src/MemcpyTauSweep.c -std=gnu99 -lgraph -lm -lX11 -o output/tausweep	

d3q7: 
	gcc -O3 -Wall -I ~/include -L ~/lib src/D3Q7.c -std=gnu99 -lgraph -lm -lX11 -o output/d3q7

loopsim:
	gcc -O3 -Wall -I ~/include -L ~/lib src/LoopSim.c -std=gnu99 -lgraph -lm -lX11 -o output/loopsim -lpthread	

coilsim: 
	gcc -O3 -Wall -I ~/include -L ~/lib src/CoilSim.c -std=gnu99 -lgraph -lm -lX11 -o output/coilsim -lpthread	

coilsimfloat: 
	gcc -O3 -Wall -I ~/include -L ~/lib src/CoilSimFloat.c -std=gnu99 -lgraph -lm -lX11 -o output/coilsimfloat -lpthread

solenoidsim:
	gcc -O3 -Wall -I ~/include -L ~/lib src/SolenoidSim.c -std=gnu99 -lgraph -lm -lX11 -o output/solenoidsim -lpthread					

.PHONY: all
all: generic
	echo "Done!"

.PHONY: clean 
clean:
	rm -f LatticeSimple
