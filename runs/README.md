## res_topology.dat

Residual resource and topology for the entire newtwork, i.e., inter-CO and intra-CO. This file is generated/updated by `run_simulation.py` after each timeslot.

```python
<co-count>
# for each co in (0, co-count]
... <renewable-energy in kW> ... # 24 data points
<node-count> <edge-count> # for the entire network
# for each node in (0, node-count]
<node-id> <co-id> <'c'/'s'> <sleep-energy> <base-energy> <per-cpu-energy> <cpu-count> 
# c: compute/server, s: switch, last two data are for servers
# for each edge in (0, edge-count]
<node-id> <node-id> <'b'/'i'> <link-cap> <link-delay> 
# b: backbone edge, i: intra-CO edge
```

## paths.dat

This file is used to list the nodes (and switches) on a path. 

> The first `server-count` entries of this file includes dummy paths to represent a path within each server that is used for inter-server chain embedding.

```python
<path-count>
# for each path in (0, path-count]
<path-id> <node-count> [...<node-id>...] <switch-count> [...<switch-id>...]
# node-ids include both switch and server ids. The switch-ids are included again to simplify data processing.
```

## n_sfc_t`x`

New SFCs for timeslot `x`, `x` $\in$ `(0, timeslot-count]`

```python
<sfc-count>
# for each sfc in (0, sfc-count]
<sfc-id> <ingress-co> <egress-co> <ttl> <vnf-count> [... <vnf-cpu-requirement> ...] <bandwidth> <max-delay>
```

## n_sfc_t`x`_srv_map

Server mapping for SFCs in n_sfc_t`x`

```python
<sfc-count>
# for each sfc in (0, sfc-count]
<sfc-id> <vnf-count> [... <node-id> ...] 
# node-id represents the servers on which the VNFs are mapped
```

## n_sfc_t`x`_path_map

Path mapping for SFCs in n_sfc_t`x`

```python
<sfc-count>
# for each sfc in (0, sfc-count]
<sfc-id> <path-id> <edge-count> [...<node-id>...]
# node-id represents both servers and switches on the path
```

## x_sfc_t`x`

Pre-existing SFCs for timeslot `x`, `x` $\in$ `(0, timeslot-count]`

```python
<sfc-count>
# for each sfc in (0, sfc-count]
<sfc-id> <ingress-co> <egress-co> <ttl> <vnf-count> [... <vnf-cpu-requirement> ...] <bandwidth> <max-delay>
```

## x_sfc_t`x`_srv_map

Server mapping for SFCs in x_sfc_t`x`

```python
<sfc-count>
# for each sfc in (0, sfc-count]
<sfc-id> <0/1> <vnf-count> [... <node-id> ...] 
# the second entry indicate whether the SFC was migrated or not. 1: migrated, 0: not migrated
# node-id represents the servers on which the VNFs are mapped
```

## x_sfc_t`x`_path_map

Path mapping for SFCs in x_sfc_t`x`

```python
<sfc-count>
# for each sfc in (0, sfc-count]
<sfc-id> <0/1> <path-id> <edge-count> [...<node-id>...]
# the second entry indicate whether the SFC was migrated or not. 1: migrated, 0: not migrated
# node-id represents both servers and switches on the path
```