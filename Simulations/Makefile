simple: LatticeSimple.c
	gcc -pg -Wall -I ~/include -L ~/lib LatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o simple2D

toroidsimple: ToroidLatticeSimple.c
	gcc -pg -Wall -I ~/include -L ~/lib ToroidLatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o toroidsimple2D

memcpysimple: MemcpyLatticeSimple.c
	gcc -O3 -Wall -I ~/include -L ~/lib MemcpyLatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o memcpysimple2D

fourcharges: FourCharges.c
	gcc -O3  -I ~/include -L ~/lib FourCharges.c -std=gnu99 -lgraph -lm -lX11 -o fourcharges

plates2d: Plates2D.c
	gcc -O3 -I ~/include -L ~/lib Plates2D.c -std=gnu99 -lgraph -lm -lX11 -o plates2D	

tausweep: MemcpyTauSweep.c
	gcc -g -Wall -I ~/include -L ~/lib MemcpyTauSweep.c -std=gnu99 -lgraph -lm -lX11 -o tausweep	

d3q7: D3Q7.c
	gcc -O3 -Wall -I ~/include -L ~/lib D3Q7.c -std=gnu99 -lgraph -lm -lX11 -o d3q7

loopsim: LoopSim.c
	gcc -O3 -Wall -I ~/include -L ~/lib LoopSim.c -std=gnu99 -lgraph -lm -lX11 -o loopsim -lpthread	

coilsim: CoilSim.c
	gcc -O3 -Wall -I ~/include -L ~/lib CoilSim.c -std=gnu99 -lgraph -lm -lX11 -o coilsim -lpthread	

coilsimfloat: CoilSimFloat.c
	gcc -O3 -Wall -I ~/include -L ~/lib CoilSimFloat.c -std=gnu99 -lgraph -lm -lX11 -o coilsimfloat -lpthread

solenoidsim: SolenoidSim.c
	gcc -O3 -Wall -I ~/include -L ~/lib SolenoidSim.c -std=gnu99 -lgraph -lm -lX11 -o solenoidsim -lpthread					

.PHONY: all
all: generic
	echo "Done!"

.PHONY: clean 
clean:
	rm -f LatticeSimple
