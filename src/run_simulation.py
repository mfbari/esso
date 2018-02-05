#!/usr/bin/python

import sys
import os
import argparse
import subprocess
import shutil

vnf_flavor_to_cpu = {}
sfcs = []
sfc_in = [set()]
sfc_out = [set()]
sfc_count = 0
timeslot_count = 0

def read_vnf_types_file(dataset_path):
    with open(os.path.join(dataset_path, "vnf_types.dat")) as f:
        type_count = int(f.readline())
        # skip over the types
        for i in range(type_count): f.readline()
        # read the flavors
        flavor_count = int(f.readline())
        for line in f.readlines():
            fid, tid, cpu, pd = [int(x) for x in line.split()]
            vnf_flavor_to_cpu[fid] = {'type_id':tid, 
                                      'cpu_count':cpu, 
                                      'proc_delay': pd}

def read_timeslots_file(dataset_path):
    global sfc_count
    global timeslot_count
    global sfc_in
    global sfc_out
    with open(os.path.join(dataset_path, "timeslots.dat")) as f:
        sfc_count = int(f.readline())
        timeslot_count = int(f.readline())
        sfc_in = [set() for t in range(timeslot_count)]
        sfc_out = [set() for t in range(timeslot_count)]
        sfc_id = 0
        for t in range(timeslot_count):
            t_sfc_count = int(f.readline())
            for s in range(t_sfc_count):
                sfc_in[t].add(sfc_id)
                values = [x for x in f.readline().split()]  
                ttl = int(values[2])
                sfc_out[t+ttl].add(sfc_id)
                sfcs.append(values)
                sfc_id += 1

if __name__ == '__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help = "path to dataset folder")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--cplex', action='store_true', 
                        help="run CPLEX code")
    group.add_argument('-t', '--tabusearch', action='store_true', 
                        help="run tabu search code")
    parser.add_argument('-i', '--id', required=True,
                        help="unique id for this run")
    parser.add_argument('-r', '--replace', action='store_true',
                        help="replace existing run dir")
    args = parser.parse_args()
    
    # path to the dataset and run folder
    dataset_path = args.path
    run_path = '../runs/'
    if args.cplex:
        run_path += 'c_'
    else:
        run_path += 't_'
    run_path += args.id

    # check whether dataset folder exists
    if not os.path.isdir(dataset_path):
        print "error: dataset folder " + dataset_path + " not found"
        exit()

    # check whether run folder exists
    if not args.replace and os.path.isdir(run_path):
        print "error: run folder " + run_path + " already exists"
        exit()

    # create run directory
    try:
        os.mkdir(run_path)
    except OSError:
        pass

    # copy data files to the run folder
    shutil.copy(os.path.join(dataset_path, 'path_switch.dat'), run_path)
    shutil.copy(os.path.join(dataset_path, 'path_edge.dat'), run_path)
    shutil.copyfile(os.path.join(dataset_path, 'init_topology.dat'), 
                    os.path.join(run_path,'res_topology.dat'))

    # build executable
    make_log = open('make_log.log', 'a')
    make_process = subprocess.Popen("make cplex", shell=True,
                    stdout=make_log, stderr=make_log)
    if make_process.wait() != 0:
        print "error: make failed"
        exit()

    # copy the executable to the run folder
    shutil.copy('esso_cplex.o', run_path)

    # read input files
    read_vnf_types_file(dataset_path)
    
    read_timeslots_file(dataset_path)

    # main simulation loop
    # x_sfcs: sfcs that can be considered for migration at a time instance
    x_sfcs = set()
    for t in range(timeslot_count):
        x_sfcs = x_sfcs.difference(sfc_out[t])
        # x_sfcs represent alive sfcs that arrived between [0,t)
        # --- --- start code for simulation
        print t, sfc_in[t], x_sfcs
        # write files for c++/CPLEX code
        with open(os.path.join(run_path, 'n_sfc_t' + str(t)), 'w') as f:
            f.write(str(len(sfc_in[t])) + "\n")
            for s in sfc_in[t]:
                f.write(" ".join(sfcs[s]) + "\n")
        with open(os.path.join(run_path, 'x_sfc_t' + str(t)), 'w') as f:
            f.write(str(len(x_sfcs)) + "\n")
            for s in x_sfcs:
                f.write(" ".join(sfcs[s]) + "\n")
        # invoke c++/CPLEX code 
        # how to represent mapping?
        # --- --- end code for simulation
        x_sfcs = x_sfcs.union(sfc_in[t])
