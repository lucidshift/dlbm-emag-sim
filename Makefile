
mygraph_install:
	cd lib/mygraph/src/ && $(MAKE) install INCDIR="../include -I lib/mygraph/include" LIBDIR="../../libgraph"

mygraph_clean:
	cd lib/mygraph/src/ && $(MAKE) clean

d3q7_sim:
	g++ -O3 -Wall -I include -I lib/mygraph/include -L lib/libgraph src/d3q7_main.cpp src/d3q7.cpp -std=c++11 -lgraph -lm -lX11 -o output/d3q7_sim -lpthread

ss_d3q7_sim:
	g++ -g -Wall -I include -I lib/mygraph/include -L lib/libgraph src/SolenoidSim_d3q7.cpp src/d3q7.cpp -std=c++11 -lgraph -lm -lX11 -o output/ss_d3q7_sim -lpthread	

general_sim:
	g++ -g -Wall -I include -I lib/mygraph/include -L lib/libgraph src/GeneralSim.cpp src/d3q7.cpp -std=c++11 -lgraph -lm -lX11 -o output/general_sim -lpthread	


.PHONY: all
all: generic
	echo "Done!"

.PHONY: clean 
clean:
	rm -f output/*
