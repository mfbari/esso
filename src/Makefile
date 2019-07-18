#
# Makefile for esso
#

CC := g++
CPPFLAGS := -std=c++11
DBGFLAGS := -g -Wall
OPZFLAGS := -O3
CPLEXINC := -L/opt/ibm/ILOG/CPLEX_Studio125/cplex/lib/x86-64_sles10_4.1/static_pic -L/opt/ibm/ILOG/CPLEX_Studio125/concert/lib/x86-64_sles10_4.1/static_pic -I/opt/ibm/ILOG/CPLEX_Studio125/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio125/concert/include 
TKCPLEXINC:= -L/home/mfbari/ILOG/CPLEX_Studio125/cplex/lib/x86-64_sles10_4.1/static_pic -L/home/mfbari/ILOG/CPLEX_Studio125/concert/lib/x86-64_sles10_4.1/static_pic -I/home/mfbari/ILOG/CPLEX_Studio125/cplex/include -I/home/mfbari/ILOG/CPLEX_Studio125/concert/include 
CPLEXLIB := -lilocplex -lconcert -lcplex -lm -lpthread -DIL_STD


cplex: esso_cplex.cpp data_store.hpp 
	$(CC) $(OPZFLAGS) $(CPPFLAGS) $(CPLEXINC) $< $(CPLEXLIB) -o esso_cplex.o

heuristic: esso_heuristic.cpp
	$(CC) $(OPZFLAGS) $(CPPFLAGS) $< -o esso_heuristic.o 

firstfit: esso_firstfit.cpp
	$(CC) $(OPZFLAGS) $(CPPFLAGS) $< -o esso_firstfit.o

dbg_heuristic: esso_heuristic.cpp
	$(CC) $(DBGFLAGS) $(CPPFLAGS) $< -o esso_heuristic.o

dbg_cplex: esso_cplex.cpp data_store.hpp
	$(CC) $(DBGFLAGS) $(CPPFLAGS) $(CPLEXINC) $< $(CPLEXLIB) -o esso_cplex.o

dbg_firstfit: esso_firstfit.cpp
	$(CC) $(DBGFLAGS) $(CPPFLAGS) $< -o esso_firstfit.o


tk_cplex: esso_cplex.cpp data_store.hpp
	$(CC) $(OPZFLAGS) $(CPPFLAGS) $(TKCPLEXINC) $< $(CPLEXLIB) -o esso_cplex.o

tk_dbg_cplex: esso_cplex.cpp data_store.hpp
	$(CC) $(DBGFLAGS) $(CPPFLAGS) $(TKCPLEXINC) $< $(CPLEXLIB) -o esso_cplex.o

process_topology: process_topology.cpp problem_instance.hpp \
		  iz_topology.hpp stop_watch.hpp \
		  esso_topology.hpp iz_priority_queue.hpp
	g++ -std=c++11 process_topology.cpp -o process_topology.o
	
dbg_process_topology: process_topology.cpp problem_instance.hpp \
		  iz_topology.hpp stop_watch.hpp \
		  esso_topology.hpp iz_priority_queue.hpp
	g++ -std=c++11 -g -Wall process_topology.cpp -o process_topology.o

.PHONY: clean
clean:
	rm -f *.o 


