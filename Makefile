simple: 
	gcc -pg -Wall -I include -I lib/mygraph/include -L lib/libgraph src/LatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o output/simple2D

toroidsimple: 
	gcc -pg -Wall -I include -I lib/mygraph/include -L lib/libgraph src/ToroidLatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o output/toroidsimple2D

memcpysimple: 
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/MemcpyLatticeSimple.c -std=gnu99 -lgraph -lm -lX11 -o output/memcpysimple2D

fourcharges: 
	gcc -O3  -I include -I lib/mygraph/include -L lib/libgraph src/FourCharges.c -std=gnu99 -lgraph -lm -lX11 -o output/fourcharges

plates2d: 
	gcc -O3 -I include -I lib/mygraph/include -L lib/libgraph src/Plates2D.c -std=gnu99 -lgraph -lm -lX11 -o output/plates2D	

tausweep: 
	gcc -g -Wall -I include -I lib/mygraph/include -L lib/libgraph src/MemcpyTauSweep.c -std=gnu99 -lgraph -lm -lX11 -o output/tausweep	

d3q7: 
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/D3Q7.c -std=gnu99 -lgraph -lm -lX11 -o output/d3q7

loopsim:
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/LoopSim.c -std=gnu99 -lgraph -lm -lX11 -o output/loopsim -lpthread	

coilsim: 
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/CoilSim.c -std=gnu99 -lgraph -lm -lX11 -o output/coilsim -lpthread	

coilsimfloat: 
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/CoilSimFloat.c -std=gnu99 -lgraph -lm -lX11 -o output/coilsimfloat -lpthread

solenoidsim:
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/SolenoidSim.c -std=gnu99 -lgraph -lm -lX11 -o output/solenoidsim -lpthread

mygraph_install:
	cd lib/mygraph/src/ && $(MAKE) install INCDIR="../include -I lib/mygraph/include" LIBDIR="../../libgraph"

mygraph_clean:
	cd lib/mygraph/src/ && $(MAKE) clean

d3q7_sim:
	gcc -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/d3q7_main.c -std=gnu99 -lgraph -lm -lX11 -o output/d3q7_sim -lpthread

.PHONY: all
all: generic
	echo "Done!"

.PHONY: clean 
clean:
	rm -f output/*
