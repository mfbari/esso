runsim:
	g++ -std=c++11 -Wall -g -O3 runsim.cpp
cplex:
	g++ -std=c++11 -Wall -L/opt/ibm/ILOG/CPLEX_Studio125/cplex/lib/x86-64_sles10_4.1/static_pic -L/opt/ibm/ILOG/CPLEX_Studio125/concert/lib/x86-64_sles10_4.1/static_pic -I/opt/ibm/ILOG/CPLEX_Studio125/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio125/concert/include esso_cplex.cpp -lilocplex -lconcert -lcplex -lm -lpthread -DIL_STD -o esso_cplex
