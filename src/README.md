## Simulation 

* `run_simulation.py`: main entry point into a simulation
* `green_shuffle.py`: shuffles the line of `green_cap.dat` file to randomize renewable energy availability at COs.   
* `run_acceptance_ratio_experiment.sh`: executes the `run_simulation.py` file in a batch mode with different parameter settings.
* `runsim.cpp`: deprecated  

## CPLEX implementation

* `esso_cplex.cpp`

## Heuristic implementation

* `esso_heuristic.hpp`
* `esso_heuristic.cpp`

## Topology util

* `process_topology.cpp`: takes a `co_topology.dat` file as input and generates the entire topology (inter and intra-CO) of the network. Basically adds server and switches at each CO node.

## Common util
* `data_store.hpp`: data structures to store network toplogy, shortest paths, current co2, energy utilization, etc.      
* `esso_topology.hpp`: represents the entire topology of the network   
* `problem_instance.hpp`: a specific input to the CPLEX and heuristic   
* `stop_watch.hpp` 
* `iz_priority_queue.hpp`  
* `iz_timer.hpp` 
* `iz_topology.hpp`        
    




