
cplex: esso_cplex.cpp data_store.hpp 
	g++ -std=c++11 -L/opt/ibm/ILOG/CPLEX_Studio125/cplex/lib/x86-64_sles10_4.1/static_pic -L/opt/ibm/ILOG/CPLEX_Studio125/concert/lib/x86-64_sles10_4.1/static_pic -I/opt/ibm/ILOG/CPLEX_Studio125/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio125/concert/include esso_cplex.cpp -lilocplex -lconcert -lcplex -lm -lpthread -DIL_STD -o esso_cplex.o

cplexww: esso_cplex.cpp data_store.hpp
	g++ -g -std=c++11 -Wall -L/opt/ibm/ILOG/CPLEX_Studio125/cplex/lib/x86-64_sles10_4.1/static_pic -L/opt/ibm/ILOG/CPLEX_Studio125/concert/lib/x86-64_sles10_4.1/static_pic -I/opt/ibm/ILOG/CPLEX_Studio125/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio125/concert/include esso_cplex.cpp -lilocplex -lconcert -lcplex -lm -lpthread -DIL_STD -o esso_cplex.o

process_topology.o: process_topology.cpp stop_watch.hpp
	g++ -std=c++11 process_topology.cpp -o process_topology.o
	
runsim: runsim.cpp
	g++ -std=c++11 -Wall -g -O3 runsim.cpp -o runsim.o

.PHONY: clean
clean:
	rm -f *.o 

dbg_process_topology: process_topology.cpp stop_watch.hpp
	g++ -std=c++11 -g -Wall process_topology.cpp -o process_topology.o

