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
                # convert from flavor-id to cpu-count
                values[5:-2] = [str(vnf_flavor_to_cpu[int(x)]['cpu_count'])
                                    for x in values[5:-2]]
                ttl = int(values[3])
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
    parser.add_argument('-d', '--dryrun', action='store_true',
                        help="only generate data")
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
    shutil.copy(os.path.join(dataset_path, 'paths.dat'), run_path)
    shutil.copyfile(os.path.join(dataset_path, 'init_topology.dat'),
                    os.path.join(run_path,'res_topology.dat'))

    # copy the executable to the run folder
    if args.cplex:
        shutil.copy('esso_cplex.o', run_path)
    else:
        shutil.copy('esso_heuristic.o', run_path)

    # read input files
    read_vnf_types_file(dataset_path)
    #print vnf_flavor_to_cpu

    read_timeslots_file(dataset_path)

    # main simulation loop
    # x_sfcs: sfcs that can be considered for migration at a time instance
    x_sfcs = set()
    # change path to the run folder
    cdir = os.getcwd()
    os.chdir(run_path)
    for t in range(timeslot_count):
        x_sfcs = x_sfcs.difference(sfc_out[t])
        # x_sfcs represent alive sfcs that arrived between [0,t)
        # --- --- start code for simulation
        print "timeslot:", t, "n_sfc:", sfc_in[t], "x_sfc:", x_sfcs
        # write files for cplex/heuristic code
        n_sfc_filename = 'n_sfc_t' + str(t)
        with open(n_sfc_filename, 'w') as f:
            f.write(str(len(sfc_in[t])) + "\n")
            for s in sfc_in[t]:
                f.write(" ".join(sfcs[s]) + "\n")
        x_sfc_filename = 'x_sfc_t' + str(t)
        with open(x_sfc_filename, 'w') as f:
            f.write(str(len(x_sfcs)) + "\n")
            for s in x_sfcs:
                f.write(" ".join(sfcs[s]) + "\n")
        # invoke cplex/heuristic code if not dryrun
        if not args.dryrun:
            with open('run.log', 'w') as exe_log:
                exe_path = './esso_cplex.o res_topology.dat paths.dat ' + \
                            n_sfc_filename + ' ' + x_sfc_filename
                exe_proc = subprocess.Popen(exe_path, shell=True,
                            stdout=exe_log, stderr=exe_log)
                if exe_proc.wait() != 0:
                    print 'failed to execute code'
                    exit()
        # end of cplex/heuristic code execution
        x_sfcs = x_sfcs.union(sfc_in[t])
        # --- --- end code for simulation

    # change pwd to esso/src
    os.chdir(cdir)
