## res_topology.dat

Residual resource and topology for the entire newtwork, i.e., inter-CO and intra-CO. This file is generated/updated by `run_simulation.py` after each timeslot.

````python
<co-count>
# for each co in (0, co-count]
... <reneable-energy in kW> ... # 24 data points
<node-count> <edge-count> # for the entire network
# for each node in (0, node-count]
<node-id> <co-id> <'c'/'s'> <sleep-energy> <base-energy> <per-cpu-energy> <cpu-count> 
# c: compute/server, s: switch, last two data are for servers
# for each edge in (0, edge-count]
<node-id> <node-id> <'b'/'i'> <link-cap> <link-delay> 
# b: backbone edge, i: intra-CO edge
````

## path_edge.dat

